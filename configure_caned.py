import os
import json
from pathlib import Path

def validate_path(path, is_dir=True, is_executable=False):
    """Validate if the path exists and meets requirements"""
    path = Path(path).resolve()
    if not path.exists():
        raise ValueError(f"Path does not exist: {path}")
    if is_dir and not path.is_dir():
        raise ValueError(f"Directory path required, but a file was provided: {path}")
    if is_executable and not os.access(path, os.X_OK):
        raise ValueError(f"Path is not executable: {path}")
    return str(path)

def configure_caned():
    """Collect tool paths from user input and generate configuration file"""
    print("Welcome to the CANED environment configuration tool! Please provide the paths for the following tools.")

    config = {}

    # Collect paths
    try:
        config["ROSETTA_BIN"] = validate_path(
            input("Enter the Rosetta executable directory (e.g., /path/to/rosetta/main/source/bin): "),
            is_dir=True
        )
        config["OPENBABEL_BIN"] = validate_path(
            input("Enter the OpenBabel executable directory (e.g., /usr/bin): "),
            is_dir=True
        )
        config["CANED_PYTHON"] = validate_path(
            input("Enter the CANED Python interpreter path (e.g., /path/to/caned_env/bin/python3): "),
            is_dir=False, is_executable=True
        )
        config["LIDMPNN_PYTHON"] = validate_path(
            input("Enter the LigandMPNN Python interpreter path (e.g., /path/to/ligandmpnn/bin/python3): "),
            is_dir=False, is_executable=True
        )
        config["ROSETTA_DATABASE"] = validate_path(
            input("Enter the Rosetta database directory (e.g., /path/to/rosetta/main/database): "),
            is_dir=True
        )
    except ValueError as e:
        print(f"Error: {e}")
        return

    # Save configuration file
    config_file = Path("core/config.json")
    try:
        with open(config_file, "w") as f:
            json.dump(config, f, indent=4)
        print(f"Configuration file saved to: {config_file}")
    except Exception as e:
        print(f"Failed to save configuration file: {e}")
        return

if __name__ == "__main__":
    configure_caned()