#!/usr/local/miniconda/bin/python
import argparse
import extract_ligand_binding_sites
import rosetta_standard_match
import pyrosetta

import os

class CustomResidue:
    """
    自定义残基类，封装 PyRosetta 的 Residue 对象并重写 is_ligand 方法
    """


    def __init__(self, residue):
        self._residue = residue  # 原始的 PyRosetta Residue 对象
    
    def is_ligand(self):
        """
        重写的 is_ligand 方法，判断残基是否是配体。
        逻辑：
        - 不是标准氨基酸或核苷酸。
        - 不是水分子（HOH）。
        - 原子数大于 1（排除单原子离子）。
        - 不参与聚合物连接（如肽键）。
        
        返回:
            bool: True 表示配体，False 表示非配体
        """
        
        STANDARD_AMINO_ACIDS = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
        }

        # 标准核苷酸的单字母代码（DNA/RNA）
        STANDARD_NUCLEOTIDES = {"A", "T", "C", "G", "U"}
        res_name = self._residue.name().split(':')[0]
        #print(res_name)
        # 如果是标准氨基酸，直接返回 False（无论是否在链端）
        
        if res_name in STANDARD_AMINO_ACIDS:
            return False
        
        # 如果是标准核苷酸，直接返回 False
        if res_name in STANDARD_NUCLEOTIDES:
            return False
        
        # 排除水分子
        if res_name == "HOH":
            print(f"{res_name} is HOH")
            return False
        
        # 排除单原子离子（如 Zn2+, Na+）
        if self._residue.natoms() <= 1:
            print(f"{res_name} is 单离子")
            return False
        
        # 非标准、非水、非单原子的残基视为配体
        return True
    
    def __getattr__(self, name):
        return getattr(self._residue, name)

class CustomPose:
    """
    自定义 Pose 类，封装 PyRosetta 的 Pose 对象并使用 CustomResidue
    """
    def __init__(self, pose):
        self._pose = pose  # 原始的 PyRosetta Pose 对象
    
    def residue(self, index):
        """
        返回封装后的 CustomResidue 对象
        """
        original_residue = self._pose.residue(index)
        return CustomResidue(original_residue)
    def get_pose(self):
        """返回原始的 Pose 对象"""
        return self._pose
    

    def __getattr__(self,name):
        return getattr(self._pose,name)
def remove_cryst1_anisou(input_pdb, output_pdb):
    """
    Remove CRYST1 and ANISOU lines from a PDB file.
    
    Args:
        input_pdb (str): Path to the input PDB file
        output_pdb (str): Path to the output processed PDB file
    Returns:
        str: Path to the output processed PDB file
    """
    with open(input_pdb, 'r') as f_in:
        lines = f_in.readlines()
    
    # Filter out CRYST1 and ANISOU lines
    filtered_lines = [line for line in lines if not (line.startswith('CRYST1') or line.startswith('ANISOU'))]
    
    with open(output_pdb, 'w') as f_out:
        f_out.writelines(filtered_lines)
    
    print(f"CRYST1 and ANISOU lines removed, saved to: {output_pdb}")
    return output_pdb

def main():
    parser = argparse.ArgumentParser(description='Generate cst from PDB files')
    
    # Existing arguments
    parser.add_argument('-i', '--input', required=True, help='Input PDB file path')
    parser.add_argument('-o', '--output', required=True, help='Output binding site directory')
    # New argument for extra_res_fa parameter
    parser.add_argument('-pa', '--params', 
                       required=True,
                       help='Path to params file for extra_res_fa')
    
    args = parser.parse_args()
    
    # Initialize PyRosetta with the command-line parameter
    pyrosetta.init(f"-extra_res_fa {args.params} -ignore_unrecognized_res false  -ignore_zero_occupancy false -load_PDB_components false ")
    
    # Process the input PDB file to remove CRYST1 and ANISOU lines before loading
    #cleaned_pdb = os.path.join(os.path.dirname(args.output), "cleaned_" + os.path.basename(args.output))
    #remove_cryst1_anisou(args.input, cleaned_pdb)
    
    # Load the cleaned PDB file
    pose = pyrosetta.pose_from_file(args.input)
    custompose = CustomPose(pose)
    bindingsite_path = os.path.join(args.output, "binding_site")  
    extract_ligand_binding_sites.get_cst(custompose, bindingsite_path)

    for file in os.listdir(bindingsite_path):
        print(file)
        if file.endswith(".pdb"):
            pdb_file_path = os.path.join(bindingsite_path, file)
            # Optionally clean generated PDB files as well
            cleaned_sub_pdb = os.path.join(bindingsite_path, "cleaned_" + file)
            remove_cryst1_anisou(pdb_file_path, cleaned_sub_pdb)
            
            pose = pyrosetta.pose_from_file(cleaned_sub_pdb)
            custompose = CustomPose(pose)
            output_file = os.path.join(args.output, file.split(".")[0] + ".cst")
            rosetta_standard_match.translate_binding_site(custompose, output_file)    
            print(f"Processing file: {pdb_file_path}")
 
if __name__ == '__main__':
    main()