#!/usr/bin/env python3
'''Match binding sites to scaffolds
using standard Rosetta matcher.
'''
import os
import subprocess
import shutil
import json

import numpy as np

import pyrosetta
from pyrosetta import rosetta

def angle(x1, x2, x3):
    '''Return the angle x1-x2-x3 in degrees.
    x1, x2 and x3 are xyzVectors.
    '''
    v1 = (x1 - x2).normalized()
    v2 = (x3 - x2).normalized()

    return np.arccos(v1.dot(v2)) * 180 / np.pi

def dihedral(x1, x2, x3, x4):
    '''Return the dihedral x1-x2-x3-x4 in degrees.
    x1, x2, x3 and x4 are xyzVectors.
    '''
    v1 = (x1 - x2).normalized()
    v2 = (x3 - x2).normalized()
    v3 = (x4 - x3).normalized()

    v1_parallel = v2.normalized()
    v1_parallel *= v1.dot(v2)
    v3_parallel = v2.normalized()
    v3_parallel *= v3.dot(v2)
    
    v1_n = (v1 - v1_parallel).normalized()
    v3_n = (v3 - v3_parallel).normalized()

    cos = v1_n.dot(v3_n)
    sin = v1_n.cross(v3_n).dot(v2)

    return np.arctan2(sin, cos) * 180 / np.pi


def atom_residue_distance(pose, res, atom_id, res2):
    '''Return the minimum heavy atom distance
    between an atom and a residue.
    '''
    min_distance = float('inf')

    for a2 in range(1, pose.residue(res2).nheavyatoms() + 1):
        dist = pose.residue(res).xyz(atom_id).distance(pose.residue(res2).xyz(a2))
        resname = pose.residue(res).name3()
        atom_i_name = pose.residue(res).atom_name(atom_id)
        res2name = pose.residue(res2).name3()

        if dist < min_distance:
            a2_name = pose.residue(res2).atom_name(a2)
            min_distance = dist
    return min_distance

def find_contact_heavy_atoms(pose, res1, res2):
    '''Find contact heavy atoms for a residue. If only one heavy atom is available, return it. 
       Otherwise, find the heavy atom closest to the opposing residue and the two heavy atoms 
       closest to this atom within the residue.
    '''
    # Get the atoms that can be considered as contact atoms for residue 1
    potential_contact_atom_ids = []

    for a in range(1, pose.residue(res1).nheavyatoms() + 1):
        # Exclude protein backbone atoms except for N, CA, and C
        # If backbone variants are included, for example OXT, the matcher application would crash
        if pose.residue(res1).is_protein():
            if pose.residue(res1).atom_is_backbone(a):
                if pose.residue(res1).atom_name(a).strip() in ['N', 'CA', 'C','O']:
                    potential_contact_atom_ids.append(a)
            elif not pose.residue(res1).is_virtual(a):  # Ignore virtual atoms for protein residues
                potential_contact_atom_ids.append(a)
        else:
            #if not pose.residue(res1).is_virtual(a):  # Ignore virtual atoms for non-protein residues
            potential_contact_atom_ids.append(a)

    # If only one potential contact atom is found, return it immediately
    if len(potential_contact_atom_ids) == 1:
        contact_atoms = [pose.residue(res1).atom_name(potential_contact_atom_ids[0])]
        return contact_atoms

    # Calculate distances for all heavy atoms to the opposing residue
    distances = []
    for a in potential_contact_atom_ids:
        distances.append((a, atom_residue_distance(pose, res1, a, res2)))

    # Get the names of the closest atom
    distances_sorted = sorted(distances, key=lambda x: x[1])
    contact_atoms = [pose.residue(res1).atom_name(distances_sorted[0][0])]
        
    
    # Get the two closest intra-residue atoms
    intra_distances = []
    for a in potential_contact_atom_ids:
        if a == distances_sorted[0][0]:
            continue
        intra_distances.append((a, pose.residue(res1).xyz(distances_sorted[0][0]).distance(
            pose.residue(res1).xyz(a))))

    # If fewer than two atoms are available for intra-residue distances, handle accordingly
    if len(intra_distances) < 2:
        # If only one additional atom is available, add it
        if len(intra_distances) == 1:
            contact_atoms.append(pose.residue(res1).atom_name(intra_distances[0][0]))
        # If no additional atoms, just return the closest atom
        return contact_atoms
    else:
        intra_distances_sorted = sorted(intra_distances, key=lambda x: x[1])
        for i in range(2):
            contact_atoms.append(pose.residue(res1).atom_name(intra_distances_sorted[i][0]))

    return contact_atoms


def generate_cst_file(site_pose, output_file):
    '''Generate the constraint file for Rosetta matching.'''
    # Find the ligand and protein residues
    ligand_residue = 0
    protein_residues = []
    for i in range(1, site_pose.size() + 1):
        if site_pose.residue(i).is_ligand():
            ligand_residue = i
        elif site_pose.residue(i).is_protein():
            protein_residues.append(i)
    assert ligand_residue != 0
    
    # Define two templates: full constraints and distance-only
    full_template = \
'''CST::BEGIN
  TEMPLATE::   ATOM_MAP: 1 atom_name: {a1_1} {a1_2} {a1_3}
  TEMPLATE::   ATOM_MAP: 1 residue3: {r1}

  TEMPLATE::   ATOM_MAP: 2 atom_name: {a2_1} {a2_2} {a2_3}
  TEMPLATE::   ATOM_MAP: 2 residue3: {r2}

  CONSTRAINT:: distanceAB: {dist:7.2f}   0.2  100.00  1        1
  CONSTRAINT::    angle_A: {angA:7.2f}  20.00 100.00  360.00   2
  CONSTRAINT::    angle_B: {angB:7.2f}  20.00 100.00  360.00   2
  CONSTRAINT::  torsion_A: {torA:7.2f}  20.00 100.00  360.00   2
  CONSTRAINT::  torsion_B: {torB:7.2f}  20.00 100.00  360.00   2
  CONSTRAINT:: torsion_AB: {toAB:7.2f}  40.00 100.00  360.00   2
CST::END'''

 
    constraint_sections = []

    for pro_res in protein_residues:
        ligand_contact_atoms = find_contact_heavy_atoms(site_pose, ligand_residue, pro_res)
        protein_contact_atoms = find_contact_heavy_atoms(site_pose, pro_res, ligand_residue)


        # Use full template with all constraints
        template = full_template
        # Get positions of all atoms
        #print(f"res1 :{site_pose.residue(pro_res).name3()} ,res2 : {site_pose.residue(ligand_residue).name3()}")

        v_r1a1 = site_pose.residue(ligand_residue).xyz(ligand_contact_atoms[0])
        v_r1a2 = site_pose.residue(ligand_residue).xyz(ligand_contact_atoms[1])
        v_r1a3 = site_pose.residue(ligand_residue).xyz(ligand_contact_atoms[2])
        v_r2a1 = site_pose.residue(pro_res).xyz(protein_contact_atoms[0])
        v_r2a2 = site_pose.residue(pro_res).xyz(protein_contact_atoms[1])
        v_r2a3 = site_pose.residue(pro_res).xyz(protein_contact_atoms[2])
        
        # Calculate all geometric constraints
        dAB = v_r1a1.distance(v_r2a1)
        aA = angle(v_r1a2, v_r1a1, v_r2a1)
        aB = angle(v_r1a1, v_r2a1, v_r2a2)
        tA = dihedral(v_r1a3, v_r1a2, v_r1a1, v_r2a1)
        tB = dihedral(v_r1a1, v_r2a1, v_r2a2, v_r2a3)
        tAB = dihedral(v_r1a2, v_r1a1, v_r2a1, v_r2a2)
        
        # Format the constraint section with all parameters
        constraint_section = template.format(
            a1_1=ligand_contact_atoms[0], a1_2=ligand_contact_atoms[1], a1_3=ligand_contact_atoms[2], 
            r1=site_pose.residue(ligand_residue).name3(),
            a2_1=protein_contact_atoms[0], a2_2=protein_contact_atoms[1], a2_3=protein_contact_atoms[2], 
            r2=site_pose.residue(pro_res).name3(),
            dist=dAB, angA=aA, angB=aB, torA=tA, torB=tB, toAB=tAB)
    

        constraint_sections.append(constraint_section)

    # Write the constraint file
    with open(output_file, 'w') as f:
        f.write('\n\n'.join(constraint_sections))
    

def translate_binding_site(pose,output_path):
    generate_cst_file(pose, output_path)
