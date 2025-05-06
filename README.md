# CANED Plugin Environment Configuration Guide

CANED is a novel enzyme design plugin built on the PyMOL platform, integrating tools such as Rosetta, PyRosetta, and LigandMPNN to provide a complete workflow from file preparation to result evaluation. This document details the installation steps for the required software, libraries, and tools, as well as the configuration of the CANED environment. Follow the steps below to set up the environment.

---

## 1. System Requirements

CANED supports the following operating systems:

- **Linux** (recommended: Ubuntu 20.04 or higher)

Hardware requirements:

Rosetta and LigandMPNN computations require significant memory, and high-throughput results need substantial disk space. The following are recommended specifications:

- **Memory**: At least 16 GB, recommended 32 GB or more
- **Storage**: At least 50 GB of free disk space (Rosetta and dependencies require significant space)

---

## 2. Installing Rosetta Software

Rosetta is a powerful molecular modeling and design software suite. CANED relies on various Rosetta applications (e.g., Rosetta Match, Enzyme Design, FastRelax). Below are detailed instructions for installing Rosetta.

### 2.1 Obtaining a Rosetta License

Rosetta requires an academic or commercial license:

1. Visit the RosettaCommons website.
2. Register for an academic user account (free) or purchase a commercial license.
3. Download the Rosetta source code or precompiled binaries (binaries are recommended for easier installation).

### 2.2 Installing Rosetta

1. **Download Rosetta**:
    - Obtain the latest Rosetta binary from the RosettaCommons download link (e.g., `rosetta_bin_linux_2024.XX.tar.gz`).
    - Extract the file:
        
        ```bash
        tar -xvzf rosetta_bin_linux_2024.XX.tar.gz
        mv rosetta_bin_linux_2024.XX /path/to/rosetta
        ```
        
2. **Verify Installation**:
    - Run the following command:
        
        ```bash
        /path/to/rosetta/main/source/bin/rosetta_scripts -h
        ```
        
        If the Rosetta help information is displayed, the installation is successful.
        

### 2.3 Installing PyRosetta

PyRosetta is the Python interface to Rosetta, used by CANED’s automatic constraint generation tool to access Rosetta functionality.

1. **Download PyRosetta**:
    - Visit the PyRosetta download page and log in with your RosettaCommons account.
    - Download the PyRosetta wheel file matching your Python version and OS (e.g., `PyRosetta4.Release.python39.ubuntu.whl`).
2. **Install PyRosetta**:
    - In a virtual environment, run:
        
        ```bash
        pip install /path/to/PyRosetta4.Release.python39.ubuntu.whl
        ```
        
3. **Verify PyRosetta**:
    - Run the Python interpreter and execute:
        
        ```python
        import pyrosetta
        pyrosetta.init()
        ```
        
        If no errors occur, PyRosetta is installed successfully.
        

---

## 3. Installing LigandMPNN

LigandMPNN is a deep learning model for structure-based sequence design, developed by Dauparas et al. (original repository: https://github.com/dauparas/LigandMPNN). CANED integrates a specific version of LigandMPNN (2024) in the `scripts/mpnn_script` directory, adapted for CANED’s workflow, eliminating the need to clone the original repository.

This section explains how to configure the integrated LigandMPNN environment.

### 3.1 Configuration Steps

1. **Verify LigandMPNN Path**:
    - Check the `scripts/mpnn_script` directory in the CANED repository to confirm it contains LigandMPNN core files (e.g., `run.py`, excluding model weights, which must be downloaded separately).
    - Example path: `/path/to/CANED/scripts/mpnn_script`.
2. **Create Python Environment**:
    - Create a dedicated Python environment for LigandMPNN (Python 3.8 or 3.9 recommended):
        
        ```bash
        python3 -m venv /path/to/ligandmpnn_env
        source /path/to/ligandmpnn_env/bin/activate
        ```
        
3. **Install Dependencies**:
    - In the activated virtual environment, install the required dependencies using the provided `requirements.txt` file in `scripts/mpnn_script`:
        
        ```bash
        pip install -r /path/to/CANED/scripts/mpnn_script/requirements.txt
        ```
        
    - Key dependencies include:
        - PyTorch (choose a version compatible with your GPU/CPU, refer to the PyTorch website).
        - Other libraries: `numpy`, `biopython`, etc. (see `requirements.txt` for details).
4. **Verify Installation**:
    - Test if LigandMPNN is operational:
        
        ```bash
        /path/to/ligandmpnn_env/bin/python3 /path/to/CANED/scripts/mpnn_script/run.py --help
        ```
        
        If help information is displayed, LigandMPNN is configured successfully.
        

### 3.2 Installing the Latest LigandMPNN (Optional)

To use the latest LigandMPNN version (e.g., for new features or fixes), install from the original repository:

1. **Clone the LigandMPNN Repository**:
    
    ```bash
    git clone https://github.com/dauparas/LigandMPNN.git
    cd LigandMPNN
    ```
    
2. **Create Python Environment**:
    
    ```bash
    python3 -m venv /path/to/ligandmpnn_env
    source /path/to/ligandmpnn_env/bin/activate
    ```
    
3. **Install Dependencies**:
    
    ```bash
    pip install torch numpy biopython
    ```
    
    Install the appropriate PyTorch version for your hardware (GPU/CPU) as per the PyTorch website.
    
4. **Verify Installation**:
    
    ```bash
    /path/to/ligandmpnn_env/bin/python3 /path/to/LigandMPNN/run.py --help
    ```
    

**Note**: The latest LigandMPNN version may not be fully compatible with CANED’s workflow. Test thoroughly or contact the CANED developers for compatibility support.

### 3.3 Licensing and Acknowledgments

The LigandMPNN version integrated into CANED adheres to the original project’s MIT License (see the LICENSE file in https://github.com/dauparas/LigandMPNN). We express gratitude to Dauparas et al. for their outstanding contribution to the structural design community. For issues or updates, visit the original repository or contact the CANED developers.

---

## 4. Installing Additional Dependency Tools

CANED integrates the following tools and libraries, which require separate installation.

### 4.1 OpenBabel

OpenBabel is used to generate substrate parameter files.

```bash
sudo apt install openbabel
```

### 4.2 Rosetta Helper Scripts

Ensure the following Rosetta scripts are available (typically provided with the Rosetta installation):

- `molfile_to_params.py`
- `clean_pdb.py`
- `gen_lig_grids`

These scripts are usually located in Rosetta’s `main/source/scripts` directory. Ensure this directory is accessible.

### 4.3 Automatic Constraint Generation Tool

CANED’s automatic constraint generation tool generates ligand binding site constraint files for Rosetta, based on the `match_ligand_binding_sites` tool developed by Kortemme-Lab (original repository: https://github.com/Kortemme-Lab/match_ligand_binding_sites). CANED integrates a modified version of this tool in the `scripts/generate_cst` directory, tailored for CANED’s workflow to support constraint file generation compatible with Rosetta Match and Enzyme Design.

### 4.3.1 CANED’s Integrated Constraint Generation Tool

CANED’s `scripts/generate_cst` directory contains a customized version of `match_ligand_binding_sites`, based on the original repository (https://github.com/Kortemme-Lab/match_ligand_binding_sites). We have made the following modifications to adapt it for CANED:

- **Optimized Input Format**: Adjusted the input file parsing logic to support CANED’s protein and ligand complex structures, simplifying input preparation.
- **Optimized Output Format**: Improved constraint file generation to ensure seamless compatibility with Rosetta’s `rosetta_match` and `enzyme_design` tools, reducing manual post-processing.
- **Integrated Processing**: Consolidated the generation of binding sites and constraint files into an automated workflow, executed via a script (e.g., `generate_cst.py`), enhancing efficiency and minimizing errors.

This version preserves the core binding site generation algorithm from the original tool while adding enhancements for CANED’s needs. It adheres to the original project’s license ([to be confirmed, e.g., MIT or other], see the LICENSE file in the original repository). We thank Kortemme-Lab for their open-source contribution to the Rosetta community.

### 4.3.2 Configuration Steps

1. **Verify Tool Path**:
    - Check the `scripts/generate_cst` directory in the CANED repository to confirm it contains the core files for the constraint generation tool (e.g., `generate_cst.py`).
    - Example path: `/path/to/CANED/scripts/generate_cst`.
2. **Install Dependencies**:
    - The constraint generation tool relies on Python libraries and PyRosetta.
    - Key dependencies include:
        - Python 3.8+ .
        - PyRosetta (see **2.3**).
        - Other libraries: `numpy`, `biopython`, etc.
3. **Verify Installation**:

Test if the constraint generation tool is operational (assuming the main script is `generate_cst.py`):

```bash
/home/user/caned_env/bin/python3 /path/to/CANED/core/scripts/generate_cst/generate_cst.py --help
```

### 4.4.3 Why Integrate the Constraint Generation Tool?

CANED integrates a modified version of `match_ligand_binding_sites` to provide a streamlined solution for generating ligand binding site constraints. The optimized input/output formats and integrated workflow make the tool user-friendly, particularly for CANED’s Rosetta Match and Enzyme Design tasks. We acknowledge Kortemme-Lab’s original algorithm for binding site generation as the foundation for our enhancements.

### 4.4.4 Installing the Original match_ligand_binding_sites (Optional)

To use the latest version of `match_ligand_binding_sites` (e.g., for new features or fixes), install from the original repository:

1. **Clone the Repository**:
    
    ```bash
    git clone https://github.com/Kortemme-Lab/match_ligand_binding_sites.git
    cd match_ligand_binding_sites
    ```
    
2. **Install Dependencies**:
    - Follow the original repository’s README to install required dependencies (typically Python libraries and Rosetta).
    - Example:
        
        ```bash
        pip install numpy biopython
        ```
        
    - Ensure Rosetta is installed (with `ROSETTA_BIN` accessible).
3. **Verify Installation**:
    - Run the original tool’s main script (e.g., `generate_constraints.py`, adjust based on actual filename):
        
        ```bash
        /home/user/caned_env/bin/python3 /path/to/match_ligand_binding_sites/generate_constraints.py --help
        ```
        

**Note**: The original version may not be fully compatible with CANED’s workflow (e.g., different input/output formats). Test thoroughly or contact the CANED developers for compatibility support.

---

## 5. CANED Environment Configuration

CANED executes Rosetta tools (e.g., Rosetta Match, Enzyme Design, FastRelax) and Python scripts (e.g., constraint generation, LigandMPNN) via shell commands. To ensure proper invocation, users must configure the paths to these tools and Python interpreters. This section provides a configuration method using a script to collect tool paths and generate a CANED configuration file.

### 5.1 Why Configure Paths?

CANED relies on the following external tools:

- **Executables**: Rosetta’s `rosetta_match`, `enzyme_design`, etc., which are compiled binaries invoked via shell commands.
- **Python Scripts**: Rosetta’s `molfile_to_params.py`, LigandMPNN’s `run.py`, etc., requiring a Python interpreter.
- **Other Tools**: OpenBabel’s `obabel` for generating substrate parameter files.

Since these tools may be installed in different paths on user systems, CANED requires their exact locations. Using a configuration file (`config.json`) to specify paths avoids hardcoding, enhancing compatibility and portability.

### 5.2 Configuration Keys and Descriptions

CANED uses a `config.json` file to store tool paths. Below are the required configuration keys and their descriptions:

| Configuration Key | Description | Example Path |
| --- | --- | --- |
| `ROSETTA_BIN` | Directory for Rosetta executables (e.g., rosetta_match) | `/path/to/rosetta/main/source/bin` |
| `OPENBABEL_BIN` | Directory for OpenBabel executables (e.g., obabel) | `/usr/bin` or `/path/to/openbabel/bin` |
| `CANED_PYTHON` | Python interpreter path for non-LigandMPNN scripts | `/path/to/caned_env/bin/python3` |
| `LIDMPNN_PYTHON` | Python interpreter path for LigandMPNN scripts | `/path/to/ligandmpnn/bin/python3` |
| `CONSTRAINT_GEN_DIR` | Directory for the constraint generation tool | `/path/to/CANED/scripts/generate_cst` |

### 5.3 Using the Configuration Script

To simplify configuration, CANED provides a Python script, `configure_caned.py`, which collects user-specified tool paths and generates the `config.json` file for internal use.

### **Firstly , you need to clone the CANED Repository**:

```bash
git clone https://github.com/OrangeCatzhang/CANED.git
cd CANED
```

The `configure_caned.py` script is located in the CANED repository root directory.

### How to Use the Configuration Script

1. **Run the Script**:
    - In the terminal, execute:
        
        ```bash
        python configure_caned.py
        ```
        
2. **Enter Paths**:
    - Follow the prompts to input the paths for Rosetta, OpenBabel, CANED Python interpreter, LigandMPNN Python interpreter, and the constraint generation tool. Example:
        
        ```
        Enter the Rosetta executable directory (e.g., /path/to/rosetta/main/source/bin): /opt/rosetta/main/source/bin
        Enter the OpenBabel executable directory (e.g., /usr/bin): /usr/bin
        Enter the CANED Python interpreter path (e.g., /path/to/caned_env/bin/python3): /home/user/caned_env/bin/python3
        Enter the LigandMPNN Python interpreter path (e.g., /path/to/ligandmpnn/bin/python3): /home/user/ligandmpnn/bin/python3
        Enter the constraint generation tool directory (e.g., /path/to/CANED/scripts/generate_cst): /path/to/CANED/scripts/generate_cst
        ```
        
3. **Script Functionality**:
    - Validates that input paths exist and meet requirements (e.g., directories or executable files).
    - Generates the `config.json` file with all paths.
    - Replaces path placeholders in the CANED main script (`main.py`) with user-specified paths (e.g., `{{ ROSETTA_BIN }}`).

**Execution Context**: Run the script in a PyMOL environment with necessary dependencies installed.

**Script Location**: CANED repository root (e.g., `/path/to/CANED/configure_caned.py`).

**Generated File**: The `config.json` file is saved in the CANED repository root.

### Generated Configuration File

After running the script, a `config.json` file is created in the CANED repository root. Example content:

```json
{    "ROSETTA_BIN": "/opt/rosetta/main/source/bin",    "OPENBABEL_BIN": "/usr/bin",    "CANED_PYTHON": "/home/user/caned_env/bin/python3",    "LIDMPNN_PYTHON": "/home/user/ligandmpnn/bin/python3",    "CONSTRAINT_GEN_DIR": "/path/to/CANED/scripts/generate_cst"}
```

CANED reads this file at runtime to access tool paths.

### 5.4 Verifying Configuration

1. **Check Configuration File**:
    - Verify that `config.json` exists and contains correct paths:
        
        ```bash
        cat config.json
        ```
        
2. **Test Tool Invocation**:
    - Use the paths from `config.json` to run the following commands, replacing paths with those in `config.json`:
        
        ```bash
        /opt/rosetta/main/source/bin/rosetta_match -h/usr/bin/obabel -h/home/user/caned_env/bin/python3 --version/home/user/ligandmpnn/bin/python3 --version/home/user/caned_env/bin/python3 /path/to/CANED/scripts/generate_cst/generate_cst.py --help
        ```
        
        If the commands execute successfully (e.g., display help or version information), the paths are valid.
        
3. **Run CANED Test**:
    - Launch PyMOL, load the CANED plugin, and execute a simple task (e.g., file preparation module) to ensure no path-related errors occur.

---

## 6. Installing the CANED Plugin

1. **Access to CANED's repository path**
    
    ```bash
    cd CANED
    ```
    
2. **Add CANED to PyMOL**:
    - Copy the CANED plugin directory to PyMOL’s plugin directory:
        
        ```bash
        cp -r /path/to/CANED ~/.pymol/startup/
        ```
        
    - Alternatively, in PyMOL, use “Plugin -> Manage Plugins -> Install” to manually install the CANED directory.
3. **Verify CANED Installation**:
    - Start PyMOL:
        
        ```bash
        pymol
        ```
        
    - Check if the “CANED” plugin option appears in the PyMOL menu bar.

---

## 7. Next Steps

After completing the above steps, the CANED plugin environment should be fully configured. You can test the plugin using example cases provided in the CANED documentation (e.g., maleate isomerase design). For issues, refer to the CANED repository’s issue page or contact the development team.

---

## 8. Reference Resources

- Rosetta Website: https://www.rosettacommons.org/
- PyRosetta Documentation: http://www.pyrosetta.org/
- LigandMPNN Repository: https://github.com/dauparas/LigandMPNN
- match_ligand_binding_sites Repository: https://github.com/Kortemme-Lab/match_ligand_binding_sites
- PyMOL Website: https://pymol.org/
- OpenBabel Documentation: http://openbabel.org/
- CANED Project: https://github.com/OrangeCatzhang/CANED