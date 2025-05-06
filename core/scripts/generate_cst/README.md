Here is the translated documentation in English for the generate_cst.py script:

Usage of generate_cst.py:

This script automatically identifies the ligand and the catalytic conformation of the protein active site from a PDB file, and then automatically generates a constraint file in CST format.

How to Use:

-i Specify the input PDB file (must contain a ligand)
-o Specify the output result

Example:

python generate_cst.py -i example/7cj5.pdb -o out
Please ensure that the PDB file you are using contains the necessary ligand information for the script to function correctly.