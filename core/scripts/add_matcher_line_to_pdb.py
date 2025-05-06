#!/usr/bin/env python3

import argparse
import pyrosetta

def extract_header_lines(template_pdb_path):
    """
    Extract HEADER, MODEL, and REMARK 666 lines from the template PDB file for the file beginning.
    
    Args:
        template_pdb_path (str): Path to the template PDB file
    
    Returns:
        list: List containing HEADER, MODEL, and REMARK 666 lines
    """
    header_lines = []
    try:
        with open(template_pdb_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(('HEADER', 'MODEL', 'REMARK 666')):
                    header_lines.append(line)
        if not header_lines:
            print("No HEADER, MODEL, or REMARK 666 lines found in template PDB, using empty header")
        return header_lines
    except Exception as e:
        print(f"Error reading template PDB file: {e}, using empty header")
        return []

def extract_ligand_header_lines(template_pdb_path):
    """
    Extract only HEADER and MODEL lines from the template PDB file for insertion before ligand HETATM records.
    
    Args:
        template_pdb_path (str): Path to the template PDB file
    
    Returns:
        list: List containing HEADER and MODEL lines
    """
    ligand_header_lines = []
    try:
        with open(template_pdb_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith(('HEADER', 'MODEL')):
                    ligand_header_lines.append(line)
        if not ligand_header_lines:
            print("No HEADER or MODEL lines found in template PDB for ligand section, using empty header")
        return ligand_header_lines
    except Exception as e:
        print(f"Error reading template PDB file for ligand header: {e}, using empty header")
        return []

def add_template_header_to_pdb(target_pdb_path, template_pdb_path, debug=False, params_file=None, output_pdb_path=None):
    """
    Add template PDB header (HEADER, MODEL, REMARK 666) to the beginning of the target PDB file,
    and insert only HEADER and MODEL lines before the ligand HETATM records.
    
    Args:
        target_pdb_path (str): Path to the target PDB file
        template_pdb_path (str): Path to the template PDB file
        debug (bool): Whether to print debug information, default is False
        params_file (str): Path to ligand params file (e.g., DGN.params), optional
        output_pdb_path (str): Path to the output PDB file, if None, overwrite the target file
    """
    # Initialize PyRosetta
    extra_options = "-run:preserve_header -beta "
    if params_file:
        extra_options += f"-extra_res_fa {params_file} "
    pyrosetta.init(extra_options=extra_options)
    
    # Extract header lines
    header_lines = extract_header_lines(template_pdb_path)
    ligand_header_lines = extract_ligand_header_lines(template_pdb_path)
    
    # Read target PDB file content
    try:
        with open(target_pdb_path, 'r') as f:
            target_lines = [line.rstrip() for line in f.readlines()]  # Strip trailing whitespace
    except Exception as e:
        print(f"Error loading target PDB file: {e}")
        return
    
    # Remove any existing HEADER, MODEL, or REMARK lines from target to avoid duplication
    target_lines_cleaned = [line for line in target_lines if not line.startswith(('HEADER', 'MODEL', 'REMARK'))]
    
    # Find the start of ligand (HETATM) records
    ligand_start_index = None
    for i, line in enumerate(target_lines_cleaned):
        if line.startswith('HETATM'):
            ligand_start_index = i
            break
    
    # Construct new PDB content
    new_pdb = []
    # Step 1: Add full header lines to the beginning
    new_pdb.extend([line + '\n' for line in header_lines])
    
    # Step 2: Add target content with ligand header before HETATM
    if ligand_start_index is not None:
        new_pdb.extend([line + '\n' for line in target_lines_cleaned[:ligand_start_index]])
        new_pdb.extend([line + '\n' for line in ligand_header_lines])
        new_pdb.extend([line + '\n' for line in target_lines_cleaned[ligand_start_index:]])
    else:
        print("No HETATM records found in target PDB, adding header only at the beginning")
        new_pdb.extend([line + '\n' for line in target_lines_cleaned])
    
    # Debug information
    if debug:
        print("Template header lines for file beginning:", header_lines)
        print("Template header lines for ligand section:", ligand_header_lines)
        print("Number of lines in new PDB:", len(new_pdb))
        remark_lines = [line.strip() for line in new_pdb if 'REMARK 666' in line]
        print("Added REMARK 666 lines:", remark_lines)
    
    # Determine output path
    if output_pdb_path is None:
        output_pdb_path = target_pdb_path
    
    # Save the modified PDB file
    with open(output_pdb_path, 'w') as f:
        f.write("".join(new_pdb))
    print(f"Saved PDB with template header to: {output_pdb_path}")

def main():
    parser = argparse.ArgumentParser(description="Add template PDB header (HEADER, MODEL, REMARK 666) to the beginning of a target PDB file, "
                                                 "and HEADER, MODEL before ligand HETATM records.")
    parser.add_argument("-i", "--target-pdb", type=str, required=True, help="Path to target PDB file")
    parser.add_argument("-t", "--template-pdb", type=str, required=True, help="Path to template PDB file")
    parser.add_argument("--debug", action="store_true", help="Print debug information (default: False)")
    parser.add_argument("--params", type=str, help="Path to ligand params file (e.g., DGN.params)")
    parser.add_argument("-o", "--output", type=str, help="Path to output PDB file (default: overwrite target)")
    
    args = parser.parse_args()

    add_template_header_to_pdb(
        target_pdb_path=args.target_pdb,
        template_pdb_path=args.template_pdb,
        debug=args.debug,
        params_file=args.params,
        output_pdb_path=args.output
    )

if __name__ == "__main__":
    main()