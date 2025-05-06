# __init_plugin__.py
import os,sys
from pymol import cmd
from . import main
from pymol.plugins import addmenuitemqt



# Gets the Python version.
python_version = sys.version_info.major
python_minor_version = "%s.%s" % (sys.version_info.major, sys.version_info.minor)
python_micro_version = "%s.%s.%s" % (sys.version_info.major, sys.version_info.minor, sys.version_info.micro)


# Checks for dependencies.
try:
    # Checks for some Qt bindings in PyMOL.
    from pymol.Qt import QtWidgets

    def showerror(title, message):
        QtWidgets.QMessageBox.critical(None, title, message)

    has_gui = "qt"
    pyqt_found = True

except ImportError:
    pyqt_found = False





# Sets the version of the CADED plugin.
__CADED_version__ = "0.9"
__revision__ = ".0"
__version__ = float(__CADED_version__ + __revision__.replace(".", ""))
CADED_plugin_name = "CADED " + __CADED_version__


plugin_dir = os.path.dirname(__file__)
sys.path.append(plugin_dir)



def __init_plugin__(app=None):
    if has_gui is None:
        print("\n# No GUI library (either Tkinter or Qt bindings) was found. CADED"
                " can not be launched.")
        return None

    # Check if a CADED main window is already open.
    if pyqt_found:
        try:
            for widget in QtWidgets.QApplication.instance().topLevelWidgets():
                if hasattr(widget, "is_CADED_main_window") and widget.isVisible():
                    title = "CADED Error"
                    message = ("CADED is already running. Please close its main"
                               " window or restart PyMOL in order to launch it again.")
                    showerror(title, message)
                    return None
        except Exception as e:
            pass

    # Checks if Python 3 is available.
    if python_version != 3:
        title = "Python Version Error"
        message = "CADED %s requires Python 3. Your current Python version is %s." % (__CADED_version__, python_micro_version)
        showerror(title, message)
        return None

    # Checks the PyMOL version.
    pymol_version = float(".".join(cmd.get_version()[0].split(".")[0:2]))
    if pymol_version < 2.3:
        title = "PyMOL Version Error"
        message = "CADED %s requires a PyMOL version of 3.1.1 or higher. Your current PyMOL version is %s." % (__CADED_version__, pymol_version)
        showerror(title, message)
        return None

    # Checks for PyQt.
    if not pyqt_found:
        title = "Import Error"
        message = "PyQt5 is not installed on your system. Please install it in order to use CADED."
        showerror(title, message)
        return None


    addmenuitemqt("CANED",lambda: main.run_plugin_gui())

    

    

