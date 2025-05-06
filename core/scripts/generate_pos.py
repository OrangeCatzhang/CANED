import pyrosetta
import numpy as np
import argparse
import os

# Define CustomResidue and CustomPose classes (consistent with what you provided earlier)
class CustomResidue:
    def __init__(self, residue):
        self._residue = residue
    
    def is_ligand(self):
        STANDARD_AMINO_ACIDS = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
        }
        STANDARD_NUCLEOTIDES = {"A", "T", "C", "G", "U"}
        res_name = self._residue.name().split(':')[0]
        if res_name in STANDARD_AMINO_ACIDS:
            return False
        if res_name in STANDARD_NUCLEOTIDES:
            return False
        if res_name == "HOH":
            return False
        if self._residue.natoms() <= 1:
            return False
        return True
    
    def __getattr__(self, name):
        return getattr(self._residue, name)

class CustomPose:
    def __init__(self, pose):
        self._pose = pose
    
    def residue(self, index):
        original_residue = self._pose.residue(index)
        return CustomResidue(original_residue)
    
    def get_pose(self):
        return self._pose
    
    def __getattr__(self, name):
        return getattr(self._pose, name)

# Main program
def main():
    # Set up command-line arguments
    parser = argparse.ArgumentParser(description="Generate pos files for residues around ligands")
    parser.add_argument("pdb_file", help="Input PDB file path")
    parser.add_argument("-d", "--distance", type=float, default=5.0, 
                        help="Distance threshold (in Å, default is 5.0)")
    parser.add_argument("-o", "--output_dir", default=".", 
                        help="Output directory (default is current directory)")
    args = parser.parse_args()

    # Initialize PyRosetta
    pyrosetta.init()

    # Load PDB file and create CustomPose object
    if not os.path.exists(args.pdb_file):
        raise FileNotFoundError(f"PDB file {args.pdb_file} does not exist")
    pose = pyrosetta.pose_from_pdb(args.pdb_file)
    custom_pose = CustomPose(pose)

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Identify all ligand residues
    ligand_residues = []
    for i in range(1, pose.total_residue() + 1):
        residue = custom_pose.residue(i)
        if residue.is_ligand():
            ligand_residues.append(i)
            print(f"Residue {i} ({residue.name()}) identified as ligand")

    if not ligand_residues:
        raise ValueError("No ligands found")

    # Generate pos file for each ligand
    for ligand_index in ligand_residues:
        ligand = custom_pose.residue(ligand_index)
        
        # Calculate ligand centroid
        atoms = [ligand.atom(i) for i in range(1, ligand.natoms() + 1)]
        coords = np.array([atom.xyz() for atom in atoms])
        centroid = np.mean(coords, axis=0)
        print(f"\nLigand {ligand_index} ({ligand.name()}) centroid: {centroid}")

        # Identify nearby residues
        nearby_residues = []
        for i in range(1, pose.total_residue() + 1):
            if i == ligand_index:
                continue  # Exclude the ligand itself
            residue = custom_pose.residue(i)
            if residue.is_protein():
                atom = residue.atom("CA")
            else:
                atom = residue.atom(1)
            distance = np.linalg.norm(np.array(atom.xyz()) - centroid)
            if distance <= args.distance:
                nearby_residues.append(i)
                print(f"Residue {i} ({residue.name()}) is {distance:.2f} Å from ligand {ligand_index}")

        # Write to pos file
        output_file = os.path.join(args.output_dir, f"active_site_ligand_{ligand_index}.pos")
        with open(output_file, "w") as f:
            for res in nearby_residues:
                f.write(f"{res}\n")
        print(f"Active site residues for ligand {ligand_index} written to {output_file}")

if __name__ == "__main__":
    main()