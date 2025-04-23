import sys
import os
import uuid
from PyQt5 import QtWidgets, uic
from pymol import cmd
from pymol.plugins import addmenuitemqt
from PyQt5.QtWidgets import QFileSystemModel, QMainWindow
from core.ui_utils import find_ui_file, select_file, switch_to_page, initialize_download_button
from core.file_utils import download_pdb, generate_cst, clean_pdb, generate_params, generate_posfile, get_cstdatabasefile, clean_temp
from core.process_utils import run_command, run_design_command, run_mpnn_command
from core.analysis_utils import analysis_design_view, analysis_mpnn_view

    
def run_plugin_gui():


    
    dialog = QtWidgets.QDialog()
    uifile = find_ui_file('CANED.ui')
    form = uic.loadUi(uifile, dialog)
    form.result_windows = {}

    dialog.rejected.connect(clean_temp)
    form.stackedWidget.setCurrentIndex(0)
    
    dialog.show()
    # 页面切换事件
    form.select_match.clicked.connect(lambda: switch_to_page(3, form))
    form.select_design.clicked.connect(lambda: switch_to_page(1, form))
    form.select_prepare.clicked.connect(lambda: switch_to_page(0, form))
    form.select_mpnn.clicked.connect(lambda: switch_to_page(2, form))

    # Prepare 功能按钮
    form.prepare_upload_pdb_2.clicked.connect(lambda: select_file('prepare_upload_pdb_2', 'prepare_show_cstpath', form))
    form.prepare_upload_pdb_3.clicked.connect(lambda: select_file('prepare_upload_pdb_3', 'prepare_show_cleanpdb', form))
    form.prepare_upload_mol2.clicked.connect(lambda: select_file('prepare_upload_mol2', 'prepare_show_mol2path', form))
    form.prepare_upload_pdbfile.clicked.connect(lambda: select_file('prepare_upload_pdbfile', 'prepare_show_pdbpath', form))
    form.select_paramsfile_2.clicked.connect(lambda: select_file('select_paramsfile_2', 'show_paramsfile', form))
    form.prepare_upload_ligandfile.clicked.connect(lambda: select_file('prepare_upload_ligandfile', 'prepare_show_ligandpath', form))

    form.prepare_get_pdbfile.clicked.connect(lambda: open_result_window(form, form.prepare_show_pdbname.text(), "get_pdb"))
    form.prepare_genrate_cst.clicked.connect(lambda: open_result_window(form, form.prepare_show_cstpath.text(), "generate_cst", form.show_paramsfile.text()))
    form.prepare_generate_cleanpdb.clicked.connect(lambda: open_result_window(form, form.prepare_show_cleanpdb.text(), "clean_pdb"))
    form.prepare_generate_params.clicked.connect(lambda: open_result_window(form, form.prepare_show_mol2path.text(), "generate_params"))
    form.prepare_generate_posfile.clicked.connect(lambda: open_result_window(form, form.prepare_show_pdbpath.text(), "generate_posfile", form.prepare_show_ligandpath.text()))
    form.prepare_get_cstfile.clicked.connect(lambda: open_result_window(form, form.prepare_show_pdbname_2.text(), "get_cstdatabasefile"))

    # Match 功能按钮
    form.select_cst_params.clicked.connect(lambda: select_file('select_cst_params', 'cst_params_entry', form))
    form.select_cst_file.clicked.connect(lambda: select_file('select_cst_file', 'cst_file_entry', form))
    form.select_pos_file.clicked.connect(lambda: select_file('select_pos_file', 'pos_file_entry', form))
    form.select_pdb_file.clicked.connect(lambda: select_file('select_pdb_file', 'pdb_file_entry', form))
    form.run_button.clicked.connect(lambda: update_file_browser(form, run_command(form)))

    # Design 功能按钮
    form.select_cst_params_2.clicked.connect(lambda: select_file('select_cst_params_2', 'cst_params_entry_2', form))
    form.select_cst_file_2.clicked.connect(lambda: select_file('select_cst_file_2', 'cst_file_entry_2', form))
    form.select_design_pdb_file.clicked.connect(lambda: select_file('select_design_pdb_file', 'design_show_pdbpath', form))
    form.design_upload_scfile.clicked.connect(lambda: select_file('design_upload_scfile', 'design_show_scfilepath_2', form))
    form.design_start_scanalysis.clicked.connect(lambda: analysis_design_view(form.design_show_scfilepath_2.text(), form))
    form.run_design_button.clicked.connect(lambda: update_file_browser(form, run_design_command(form)))

    # MPNN 功能按钮
    form.select_checkpoint_file.clicked.connect(lambda: select_file('select_checkpoint_file', 'checkpoint_file_entry', form))
    form.select_pdb_file_2.clicked.connect(lambda: select_file('select_pdb_file_2', 'pdb_file_entry_2', form))
    form.select_model_type.currentIndexChanged.connect(lambda: form.model_type_entry.setText(f"current select: {form.select_model_type.currentText()}"))
    form.model_type_entry.setText("default model type: ligand mpnn")
    form.select_cst_file_3.clicked.connect(lambda: select_file('select_cst_file_3', 'show_cstpath', form))
    form.select_paramsfile.clicked.connect(lambda: select_file('select_paramsfile', 'show_paramsfilepath', form))
    form.select_score_button.clicked.connect(lambda: select_file('select_score_button', 'mpnn_show_scpath', form))
    form.analysisscore.clicked.connect(lambda: analysis_mpnn_view(form.mpnn_show_scpath.text(), form))
    form.run_button_2.clicked.connect(lambda: run_mpnn_command(form))

    

def open_result_window(form, file_path, type, option_file2_path=None):
    """打开结果窗口"""
    window_id = uuid.uuid4().hex
    new_window = ResultWindow(file_path, type, option_file2_path, window_id, form)
    form.result_windows[window_id] = new_window
    new_window.show()
    return window_id

class ResultWindow(QMainWindow):
    def __init__(self, file_path, type, option_file2_path=None, window_id=None, form=None):
        super().__init__()
        self.file_path = file_path
        self.type = type
        self.option_file2_path = option_file2_path
        self.window_id = window_id
        self.form = form

        uifile = find_ui_file('result_ui.ui')
        self.ui = uic.loadUi(uifile, self)
        

        output_path = {
            "get_pdb": download_pdb,
            "generate_cst": lambda x: generate_cst(x, self.option_file2_path),
            "clean_pdb": clean_pdb,
            "generate_params": generate_params,
            "generate_posfile": lambda x: generate_posfile(x, self.option_file2_path),
            "get_cstdatabasefile": get_cstdatabasefile
        }[type](file_path)

        self.file_list = [f for root, _, files in os.walk(output_path) for f in files]
        for file in self.file_list:
            self.ui.listWidget.addItem(file)

        self.ui.pushButton.clicked.connect(lambda: self.download_file(output_path))


    def download_file(self, output_path):
        selected_file = self.ui.listWidget.currentItem().text()
        current_file = os.path.join(output_path, selected_file)
        save_path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save file", selected_file)
        if save_path and os.path.exists(current_file):
            with open(current_file, 'rb') as f_in, open(save_path, 'wb') as f_out:
                f_out.write(f_in.read())
            print(f"File saved to: {save_path}")

    def closeEvent(self, event):
        
        if self.window_id in self.form.result_windows:
            del self.form.result_windows[self.window_id]
        event.accept()
        

def update_file_browser(form, result_path):
    """更新文件浏览器"""
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

if __name__ == "__main__":
    addmenuitemqt("CADPD Plugin", run_plugin_gui)