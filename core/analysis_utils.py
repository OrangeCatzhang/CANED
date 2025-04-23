import pandas as pd
from PyQt5.QtCore import QAbstractTableModel, Qt
from io import StringIO

class PandasModel(QAbstractTableModel):
    def __init__(self, data):
        super().__init__()
        self._data = data

    def rowCount(self, parent=None):
        return self._data.shape[0]

    def columnCount(self, parent=None):
        return self._data.shape[1]

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid() and role == Qt.DisplayRole:
            return str(self._data.iloc[index.row(), index.column()])
        return None

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._data.columns[section])
            else:
                return str(section)
        return None

    def sort(self, column, order):
        self.layoutAboutToBeChanged.emit()
        self._data = self._data.sort_values(self._data.columns[column], ascending=(order == Qt.AscendingOrder))
        self.layoutChanged.emit()

def read_and_display_data(file_path, view_widget):
    """读取并显示数据"""
    with open(file_path, 'r') as f:
        lines = f.readlines()

    is_design_sc = False
    is_mpnn = False
    for line in lines:
        if line.strip().startswith("total_score") or line.strip().startswith("cst_score"):
            is_design_sc = True
            break
        elif line.startswith('>'):
            is_mpnn = True
            break

    if is_design_sc:
        header_index = next((i for i, line in enumerate(lines) if line.strip().startswith("total_score") or line.strip().startswith("cst_score")), None)
        if header_index is not None:
            table_data = lines[header_index:]
            df = pd.read_csv(StringIO(''.join(table_data)), sep='\s+')
        else:
            print("Header row not found in Design SC format")
            return
    elif is_mpnn:
        data = []
        for line in lines:
            if line.startswith('>'):
                parts = line[1:].strip().split(',')
                entry = {}
                for part in parts:
                    if '=' in part:
                        key, value = part.split('=', 1)
                        entry[key.strip()] = value.strip()
                    else:
                        entry['pdb name'] = part.strip()
                data.append(entry)
        df = pd.DataFrame(data)
    else:
        print("Unknown file format")
        return

    model = PandasModel(df)
    view_widget.setModel(model)
    view_widget.setSortingEnabled(True)

def analysis_design_view(ligand_score_path, form):
    """分析并显示 Design 数据"""
    read_and_display_data(ligand_score_path, form.design_show_scanalyresults)

def analysis_mpnn_view(ligand_score_path, form):
    """分析并显示 MPNN 数据"""
    read_and_display_data(ligand_score_path, form.analysis_view)
