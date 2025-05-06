import os
import shutil
import subprocess
import requests
import json
from pathlib import Path

# Load configuration from config.json
def load_config(config_path="config.json"):
    """Load tool paths from config.json."""
    config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), config_path)
    try:
        with open(config_file, "r") as f:
            config = json.load(f)
        return config
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file {config_path} not found. Run configure_caned.py to generate it.")
    except json.JSONDecodeError:
        raise ValueError(f"Invalid JSON in {config_path}.")

plugin_dir = os.path.dirname(os.path.abspath(__file__))
# Load paths from config.json
config = load_config()
ROSETTA_BIN = config["ROSETTA_BIN"]
OPENBABEL_BIN = config["OPENBABEL_BIN"]
CANED_PYTHON = config["CANED_PYTHON"]
CANED_DIR = plugin_dir


# Define script paths using config
SCRIPT_PATHS = {
    "molfile_to_params_script": os.path.join(CANED_DIR, "scripts/prepare_script/molfile_to_params.py"),
    "clean_pdb_script": os.path.join(CANED_DIR, "scripts/prepare_script/clean_pdb.py"),
    "generate_cst_script": os.path.join(CANED_DIR, "scripts/generate_cst/generate_cst.py"),  
    "gen_lig_grids_script": os.path.join(ROSETTA_BIN, "gen_lig_grids.static.linuxgccrelease"),  
    "ECNum_relate_PDBNum_script": os.path.join(CANED_DIR, "scripts/prepare_script/ECNum_relate_PDBNum.py"), 
    "index_txt": os.path.join(CANED_DIR, "scripts/prepare_script/index.txt")  
}

def download_pdb(pdb_id):
    """Download PDB file."""
    original_dir = os.getcwd()
    output_path = os.path.join(original_dir, "./CADPD_tmp/prepare/get_pdb")
    os.makedirs(output_path, exist_ok=True)
    os.chdir(output_path)

    pdb_url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    try:
        response = requests.get(pdb_url)
        output_file = f"{pdb_id}.pdb"
        with open(output_file, 'wb') as f:
            f.write(response.content)
        print(f"PDB file {pdb_id}.pdb downloaded to {output_path}")
    except Exception as e:
        print(f"Unable to download PDB file {pdb_id}.pdb: {e}")
        return None
    finally:
        os.chdir(original_dir)
    return output_path

def generate_cst(input_file, params_file):
    """Generate CST file."""
    output_path = os.path.join(os.getcwd(), "./CADPD_tmp/prepare/generate_cst")
    os.makedirs(output_path, exist_ok=True)
    try:
        script_path = SCRIPT_PATHS["generate_cst_script"]
        subprocess.run([CANED_PYTHON, script_path, "-i", input_file, "-pa", params_file, "-o", output_path], check=True)
        print(f"CST file generated to: {output_path}")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Generate CST file failed: {e}")
        return None
    return output_path

def clean_pdb(input_file):
    """Clean PDB file."""
    original_dir = os.getcwd()
    output_path = os.path.join(original_dir, "./CADPD_tmp/prepare/clean")
    os.makedirs(output_path, exist_ok=True)
    os.chdir(output_path)
    try:
        script_path = SCRIPT_PATHS["clean_pdb_script"]
        subprocess.run([CANED_PYTHON, script_path, input_file, "A"], check=True)
        print(f"Cleaned PDB file generated to: {output_path}")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Clean PDB failed: {e}")
        return None
    finally:
        os.chdir(original_dir)
    return output_path

def generate_params(input_file):
    """Generate parameter file."""
    original_dir = os.getcwd()
    output_path = os.path.join(original_dir, "./CADPD_tmp/prepare/generate_params")
    os.makedirs(output_path, exist_ok=True)
    os.chdir(output_path)
    try:
        script_path = SCRIPT_PATHS["molfile_to_params_script"]
        ligand_name = os.path.splitext(os.path.basename(input_file))[0]
        sdf_output = os.path.join(output_path, f"{ligand_name}_conformers.sdf")
        obabel_cmd = [os.path.join(OPENBABEL_BIN, "obabel"), input_file, "-O", sdf_output, "-h", "--conformer", "--nconf", "100", "--systematic", "--writeconformers"]
        subprocess.run(obabel_cmd, check=True)
        print(f"SDF file generated: {sdf_output}")

        params_output = os.path.join(output_path, ligand_name)
        command = [CANED_PYTHON, script_path, sdf_output, "--clobber", "--keep-names", "--long-names", "--kinemage=file", "-n", ligand_name, "--conformers-in-one-file"]
        subprocess.run(command, check=True)
        print(f"Parameter files generated in: {output_path}")

        params_file = os.path.join(output_path, f"{ligand_name}.params")
        conformers_pdb = os.path.join(output_path, f"{ligand_name}_conformers.pdb")
        if os.path.exists(params_file) and os.path.exists(conformers_pdb):
            with open(params_file, "r") as f:
                lines = f.readlines()
            if not any("PDB_ROTAMERS" in line for line in lines):
                with open(params_file, "a") as f:
                    f.write(f"PDB_ROTAMERS {conformers_pdb}\n")
                print(f"Added PDB_ROTAMERS line to {params_file}")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Generate params failed: {e}")
        return None
    finally:
        os.chdir(original_dir)
    return output_path

def generate_posfile(input_file, input_file2):
    """Generate POS file."""
    original_dir = os.getcwd()
    output_path = os.path.join(original_dir, "./CADPD_tmp/prepare/generate_posfile")
    os.makedirs(output_path, exist_ok=True)
    os.chdir(output_path)
    try:
        script_path = SCRIPT_PATHS["gen_lig_grids_script"]
        shutil.copy(input_file, output_path)
        shutil.copy(input_file2, output_path)
        base_name = os.path.basename(input_file)
        base_name2 = os.path.basename(input_file2)
        subprocess.run([script_path, "-s", base_name, base_name2, "-grid_active_res_cutoff", "5.0"], check=True)
        os.remove(base_name)
        os.remove(base_name2)
        print(f"POS file generated to: {output_path}")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Generate POS file failed: {e}")
        return None
    finally:
        os.chdir(original_dir)
    return output_path

def get_cstdatabasefile(target_name):
    """Get CST database file."""
    original_dir = os.getcwd()
    output_path = os.path.join(original_dir, "./CADPD_tmp/prepare/getcst")
    os.makedirs(output_path, exist_ok=True)
    try:
        script_path = SCRIPT_PATHS["ECNum_relate_PDBNum_script"]
        import re
        is_ec = re.match(r'^\d+\.\d+\.\d+\.\d+$', target_name) or target_name.startswith('EC')
        additional_arguments = ["-se", target_name] if is_ec else ["-sp", target_name]
        command = [CANED_PYTHON, script_path, "-i", SCRIPT_PATHS["index_txt"], *additional_arguments]
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        ec_paths = result.stdout.strip().splitlines()[1:]
        if not ec_paths:
            print("No EC paths found.")
            return None
        for path in ec_paths:
            dir_name = os.path.basename(path)
            shutil.copytree(path, os.path.join(output_path, dir_name), dirs_exist_ok=True)
        return output_path
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"Get CST database file failed: {e}")
        return None
    finally:
        os.chdir(original_dir)

def clean_temp():
    """Clean temporary directory."""
    temp_dir = './CADPD_tmp'
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
            print(f"Successfully deleted temporary directory: {temp_dir}")
        except Exception as e:
            print(f"Error deleting temporary directory: {e}")
    
