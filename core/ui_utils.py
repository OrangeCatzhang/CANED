import os
from PyQt5 import QtWidgets, uic
from PyQt5.QtWidgets import QFileDialog, QMessageBox

def find_ui_file(filename):
    """查找 UI 文件的路径"""
    plugin_dir = os.path.dirname(os.path.abspath(__file__))
    ui_path = os.path.join(plugin_dir,"resources",filename)
    print(ui_path)
    if os.path.exists(ui_path):
        return os.path.abspath(ui_path)
    raise FileNotFoundError(f"Unable to find the UI file: {filename}")


def select_file(button_name, line_edit_name, form):
    """选择文件并更新对应的 QLineEdit"""
    filename, _ = QFileDialog.getOpenFileName(form, "Select File")
    if filename:
        getattr(form, line_edit_name).setText(filename)

def switch_to_page(page_index, form):
    """切换 stackedWidget 的页面"""
    form.stackedWidget.setCurrentIndex(page_index)

def download_selected_file(file_system_model, file_tree_view):
    """下载选中的文件"""
    selected_indexes = file_tree_view.selectionModel().selectedIndexes()
    if not selected_indexes:
        QMessageBox.warning(None, "No file selected!", "Please select a file to download.", QMessageBox.Ok)
        return

    file_path = file_system_model.filePath(selected_indexes[0])
    if not os.path.exists(file_path):
        QMessageBox.warning(None, "Error", f"Invalid file path: {file_path}!", QMessageBox.Ok)
        return

    download_directory = QFileDialog.getExistingDirectory(None, "Select download directory", os.getcwd())
    if not download_directory:
        return

    import shutil
    file_name = os.path.basename(file_path)
    target_path = os.path.join(download_directory, file_name)
    try:
        shutil.copy(file_path, target_path)
        QMessageBox.information(None, "Success", f"File successfully downloaded to {target_path}", QMessageBox.Ok)
    except Exception as e:
        QMessageBox.critical(None, "Error", f"File download failed: {e}", QMessageBox.Ok)

def initialize_download_button(form, file_system_model):
    """初始化下载按钮"""
    download_button = QtWidgets.QPushButton("Download File", form)
    download_button.clicked.connect(lambda: download_selected_file(file_system_model, form.file_tree_view))
    form.layout().addWidget(download_button)