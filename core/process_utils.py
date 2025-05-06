import os
import subprocess
import threading
from PyQt5 import QtWidgets, QtCore
from glob import glob
import pandas as pd
import json

from pymol import cmd
from .ui_utils import initialize_download_button
from PyQt5.QtWidgets import QFileSystemModel

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
config = load_config()
ROSETTA_BIN = config["ROSETTA_BIN"]
CANED_PYTHON = config["CANED_PYTHON"]
LIDMPNN_PYTHON = config["LIDMPNN_PYTHON"]
ROSETTA_DATABASE = config["ROSETTA_DATABASE"]
CANED_DIR = plugin_dir

# Define script paths using config
SCRIPT_PATHS = {
    "match_script": os.path.join(ROSETTA_BIN, "match.static.linuxgccrelease"),
    "design_script": os.path.join(ROSETTA_BIN, "enzyme_design.static.linuxgccrelease"),
    "mpnn_script": os.path.join(CANED_DIR, "scripts/mpnn_script/run.py"),
    "add_matcher_line_to_pdb_script": os.path.join(CANED_DIR, "scripts/add_matcher_line_to_pdb.py"),
    "FastRelax_script": os.path.join(CANED_DIR, "scripts/FastRelax.py")
}

class ProgressDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Task Progress")
        self.setFixedSize(300, 150)
        layout = QtWidgets.QVBoxLayout()
        self.status_label = QtWidgets.QLabel("Starting task...")
        layout.addWidget(self.status_label)
        self.cancel_button = QtWidgets.QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.cancel_task)
        layout.addWidget(self.cancel_button)
        self.setLayout(layout)
        self.is_canceled = False
        self.processes = []

    def update_status(self, message):
        self.status_label.setText(message)
        QtWidgets.QApplication.processEvents()

    def cancel_task(self):
        self.is_canceled = True
        self.update_status("Canceling task...")
        for process in self.processes:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    process.kill()
        self.close()

    def closeEvent(self, event):
        self.is_canceled = True
        event.accept()

class TaskSignalEmitter(QtCore.QObject):
    task_finished = QtCore.pyqtSignal(bool)

def run_command(form):
    """Run Rosetta Match command with progress dialog and threading."""
    dialog = ProgressDialog()
    dialog.show()
    signal_emitter = TaskSignalEmitter()

    def run_task():
        try:
            dialog.update_status("Starting Rosetta Match task...")
            ligand_params = form.cst_params_entry.text()
            cst_file = form.cst_file_entry.text()
            pos_file = form.pos_file_entry.text()
            pdb_file = form.pdb_file_entry.text()
            script_path = SCRIPT_PATHS["match_script"]
            ligand_name = os.path.splitext(os.path.basename(ligand_params))[0]
            command = [
                script_path,
                "-extra_res_fa", ligand_params,
                "-match:geometric_constraint_file", cst_file,
                "-match:scaffold_active_site_residues", pos_file,
                "-s", pdb_file,
                "-match:lig_name", ligand_name,
                "-in:ignore_unrecognized_res",
                "-ex1", "-ex2",
                "-match:consolidate_matches", "True",
                "-match:output_matches_per_group", "1",
                "-ignore_zero_occupancy", "false"
            ]
            os.makedirs("./CADPD_tmp/match_results", exist_ok=True)
            os.chdir("./CADPD_tmp/match_results")
            process = subprocess.Popen(command)
            dialog.processes.append(process)
            process.wait()
            os.chdir("../../")
            if dialog.is_canceled:
                return
            dialog.update_status("Rosetta Match task completed successfully")
            signal_emitter.task_finished.emit(True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"Match error: {e}")
            dialog.update_status(f"Error: {e}")
            signal_emitter.task_finished.emit(False)

    def on_task_finished(success):
        if success and not dialog.is_canceled:
            browserpath = "CADPD_tmp/match_results/"
            update_file_browser(form, browserpath)
            QtCore.QTimer.singleShot(1000, dialog.close)

    signal_emitter.task_finished.connect(on_task_finished)
    task_thread = threading.Thread(target=run_task)
    task_thread.start()

def run_design_command(form):
    """Run Rosetta Design command with progress dialog and threading."""
    dialog = ProgressDialog()
    dialog.show()
    signal_emitter = TaskSignalEmitter()

    def run_task():
        try:
            dialog.update_status("Starting Rosetta Design task...")
            ligand_params = form.cst_params_entry_2.text()
            cst_file = form.cst_file_entry_2.text()
            design_pdb_file = form.design_show_pdbpath.text()
            choosdesign_choosenum = form.design_choosenum.text() or "15"
            script_path = SCRIPT_PATHS["design_script"]
            protocol_file = os.path.join("scripts/", "enzdes_new.xml")
            command = [
                script_path,
                "-s", design_pdb_file,
                "-enzdes::cstfile", cst_file,
                "-extra_res_fa", ligand_params,
                "-run::preserve_header",
                "-parser:protocol", protocol_file,
                "-database", ROSETTA_DATABASE,
                "-enzdes::detect_design_interface",
                "-enzdes::cut1", "6.0",
                "-enzdes::cut2", "8.0",
                "-enzdes::cut3", "10.0",
                "-enzdes::cut4", "12.0",
                "-mute", "core.io.database",
                "-jd2::enzdes_out",
                "-nstruct", choosdesign_choosenum,
                "-jd2:ntrials", "1",
                "-out:file:o", "score.sc",
                "-ignore_zero_occupancy", "false"
            ]
            os.makedirs("./CADPD_tmp/design_results", exist_ok=True)
            os.chdir("./CADPD_tmp/design_results")
            process = subprocess.Popen(command)
            dialog.processes.append(process)
            process.wait()
            os.chdir("../../")
            if dialog.is_canceled:
                return
            dialog.update_status("Rosetta Design task completed successfully")
            signal_emitter.task_finished.emit(True)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"Design error: {e}")
            dialog.update_status(f"Error: {e}")
            signal_emitter.task_finished.emit(False)

    def on_task_finished(success):
        if success and not dialog.is_canceled:
            browserpath = "CADPD_tmp/design_results/"
            update_file_browser(form, browserpath)
            QtCore.QTimer.singleShot(1000, dialog.close)

    signal_emitter.task_finished.connect(on_task_finished)
    task_thread = threading.Thread(target=run_task)
    task_thread.start()

def parse_pdb_constraints(pdb_file):
    """Parse constraints from PDB file."""
    remark_residues = []
    link_residues = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("REMARK 666 MATCH TEMPLATE"):
                parts = line.split()
                chain = parts[-5]
                residue_num = parts[-3]
                remark_residues.append(f"{chain}{residue_num}")
            elif line.startswith("LINK"):
                parts = line.split()
                chain = parts[3]
                residue_num = parts[4]
                link_residues.append(f"{chain}{residue_num}")
    if not remark_residues and not link_residues:
        raise ValueError("No REMARK 666 or LINK constraints found in PDB file!")
    return f"'{' '.join(link_residues or remark_residues)}'"

def run_mpnn_command(form):
    """Run MPNN command."""
    checkpoint_file = form.checkpoint_file_entry.text()
    pdb_infile = form.pdb_file_entry_2.text()
    model_type = form.select_model_type.currentText()
    chain = form.mpnn_show_chain.text()
    cst_file = form.show_cstpath.text()
    params_file = form.show_paramsfilepath.text()
    mpnn_choosenum = form.mpnn_choosenum.text() or "1"
    choose_cstcutoff = form.choose_cstcutoff.text() or "0"

    dialog = ProgressDialog()
    dialog.show()
    signal_emitter = TaskSignalEmitter()

    def run_task():
        try:
            dialog.update_status("Starting MPNN task...")
            remain_cs = parse_pdb_constraints(pdb_infile)
            mpnn_script = SCRIPT_PATHS["mpnn_script"]
            command = [
                LIDMPNN_PYTHON,
                mpnn_script,
                "--pdb_path", pdb_infile,
                "--checkpoint_ligand_mpnn", checkpoint_file,
                "--model_type", model_type,
                "--seed", "37",
                "--batch_size", mpnn_choosenum,
                "--temperature", "0.05",
                "--save_stats", "1",
                "--chains_to_design", chain,
                "--out_folder", "./CADPD_tmp/mpnn_results",
                "--fixed_residues", remain_cs
            ]
            process = subprocess.Popen(command)
            process.wait()

            mpnn_output_dir = "./CADPD_tmp/mpnn_results"
            tmp_dir = os.path.join(mpnn_output_dir, "tmp")
            fastrelax_output_dir = os.path.join(mpnn_output_dir, "fastrelax_results")
            os.makedirs(tmp_dir, exist_ok=True)
            os.makedirs(fastrelax_output_dir, exist_ok=True)

            pdb_files = glob(os.path.join(mpnn_output_dir, "./backbones/*.pdb"))
            for pdb_path in pdb_files:
                pdb_basename = os.path.splitext(os.path.basename(pdb_path))[0]
                tmp_pdb_path = os.path.join(tmp_dir, f"{pdb_basename}_with_header.pdb")
                matcher_script = SCRIPT_PATHS["add_matcher_line_to_pdb_script"]
                matcher_command = [
                    CANED_PYTHON,
                    matcher_script,
                    "-i", pdb_path,
                    "-t", pdb_infile,
                    "--params", params_file,
                    "--debug",
                    "-o", tmp_pdb_path
                ]
                subprocess.Popen(matcher_command).wait()

                output_path = os.path.join(fastrelax_output_dir, f"{pdb_basename}_relaxed")
                fastrelax_script = SCRIPT_PATHS["FastRelax_script"]
                fastrelax_command = [
                    CANED_PYTHON,
                    fastrelax_script,
                    "-i", tmp_pdb_path,
                    "-c", cst_file,
                    "-o", output_path,
                    "-pa", params_file,
                    "-cut", choose_cstcutoff
                ]
                process = subprocess.Popen(fastrelax_command)
                dialog.processes.append(process)
                process.wait()
                if dialog.is_canceled:
                    return

            score_files = glob(os.path.join(fastrelax_output_dir, "*.sc"))
            if score_files and not dialog.is_canceled:
                combined_scores_file = os.path.join(fastrelax_output_dir, "combined_scores.sc")
                combined_df = pd.concat([pd.read_csv(f, sep='\s+', engine='python') for f in score_files], ignore_index=True)
                combined_df.to_csv(combined_scores_file, sep='\t', index=False)
                print(f"Combined scores saved to {combined_scores_file}")

            dialog.update_status("Task completed successfully")
            signal_emitter.task_finished.emit(True)
        except (subprocess.CalledProcessError, FileNotFoundError, ValueError) as e:
            print(f"MPNN task error: {e}")
            dialog.update_status(f"Error: {e}")
            signal_emitter.task_finished.emit(False)

    def on_task_finished(success):
        if success and not dialog.is_canceled:
            browserpath = "CADPD_tmp/mpnn_results/"
            update_file_browser(form, browserpath)
            QtCore.QTimer.singleShot(1000, dialog.close)

    signal_emitter.task_finished.connect(on_task_finished)
    task_thread = threading.Thread(target=run_task)
    task_thread.start()



def update_file_browser(form, result_path):
    """Update file browser (placeholder function, implement as needed)."""
    print(f"Updating file browser with path: {result_path}")
    if result_path:
        file_system_model = QFileSystemModel()
        file_system_model.setRootPath(result_path)
        form.file_tree_view.setModel(file_system_model)
        form.file_tree_view.setRootIndex(file_system_model.index(result_path))
        form.file_tree_view.selectionModel().selectionChanged.connect(lambda selected, _: load_pdb_structure(file_system_model, selected))
        initialize_download_button(form, file_system_model)

def load_pdb_structure(model, selected):
    """加载选中的 PDB 结构"""
    indexes = selected.indexes()
    if indexes and not model.isDir(indexes[0]):
        cmd.load(model.filePath(indexes[0]))
        cmd.show('cartoon')
        cmd.color('auto')
        cmd.reset()

