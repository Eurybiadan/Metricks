# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import sys
import tkinter.filedialog

import scipy
from PySide6 import QtCore, QtWidgets, QtGui

class PygmyFeeder(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()

        self.button = QtWidgets.QPushButton("Click me damnit")
        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.button)


class PygmyMetricks():
    def __init__(self):
        super().__init__()

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    app = QtWidgets.QApplication([])

    widget = PygmyFeeder()
    widget.resize(800, 600)
    widget.show()

    sys.exit(app.exec())
