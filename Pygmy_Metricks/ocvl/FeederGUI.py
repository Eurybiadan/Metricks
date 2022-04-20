
#  Copyright (c) 2021. Robert F Cooper
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import sys, time

from PySide6 import QtGui, QtWidgets, QtCore
from PySide6.QtWidgets import QApplication, QWizard, QWizardPage, QVBoxLayout, QGridLayout, QPushButton, QLabel, \
    QLineEdit, QFileDialog, QListWidget, QAbstractItemView, QProgressBar, QWidget
from PySide6.QtCore import Qt, Signal, Slot, Property

class PygmyFeeder(QWizard):
    def __init__(self, parent=None):
        QWizard.__init__(self,parent)
        self.setWizardStyle(QWizard.ModernStyle)  # added to get the title to be formatted correctly
        self.setWindowTitle("Welcome to OCVL's Metricks Master (Pygmy Python editon)")
        self.addPage(WelcomePage(self))
        self.addPage(MetricSelectPage(self))  # added the second page in the wizard
        self.addPage(UnitSelect(self))  # added the third page in the wizard
        self.addPage(ResultsSaveLocation(self))  # added the fourth page in the wizard
        self.addPage(Calculate(self))  # added the fifth page in the wizard


class WelcomePage(QWizardPage):
    butt_signal = Signal(str)  # Make a signal, pass it

    def __init__(self, parent=None):
        QWizardPage.__init__(self, parent)
        self.setTitle("Select the data source(s) to analyze:")

        self.label = QLabel("<b>Images are not required.</b><br> However, if they"+
                            " are not present, then bounded metrics will be calculated from the bounding rectangle"
                            " of the provided coordinates.")
        self.label.setWordWrap(True)
        self.label.setTextFormat(Qt.RichText)
        self.label.setAlignment(Qt.AlignCenter)

        self.file_dir_form = QGridLayout()

        self.coord_label = QLineEdit()
        self.coord_label.setReadOnly(True)
        self.coord_label.text()
        self.coord_butt = QPushButton("Select...")

        self._coord_path = ""

        self.image_label = QLineEdit()
        self.image_label.setReadOnly(True)
        #self.image_label.setMaximumWidth()
        self.image_butt = QPushButton("Select...")
        self.image_path = ""

        self.file_dir_form.addWidget(self.coord_label, 0, 0)
        self.file_dir_form.addWidget(self.coord_butt, 0, 1)
        self.file_dir_form.addWidget(self.image_label, 1, 0)
        self.file_dir_form.addWidget(self.image_butt, 1, 1)

        self.v_layout = QVBoxLayout()
        self.v_layout.setSpacing(32)
        self.v_layout.addWidget(self.label)
        self.v_layout.addLayout(self.file_dir_form)

        self.file_dir_form.setSpacing(4)

        self.setLayout(self.v_layout)

        self.coord_butt.clicked.connect(self.select_coord_path)
        self.image_butt.clicked.connect(self.select_image_path)



    def readCoordPath(self):
        return self._coord_path

    def setCoordPath(self, val):
        self._coord_path = val

    coord_path = Property(str, readCoordPath, setCoordPath, notify=butt_signal)

    @Slot()
    def select_coord_path(self):
        self._coord_path = QFileDialog.getExistingDirectory(parent=self,
                                                           caption="Select the folder containing the coordinates of interest.",
                                                           options=QFileDialog.ShowDirsOnly)
        # added to display the chosen path
        if self._coord_path:
            self.coord_label.setText(self._coord_path)


    @Slot()
    def select_image_path(self):

        self.image_path = QFileDialog.getExistingDirectory(parent=self,
                                                           caption="Select the folder containing the images of interest.",
                                                           options=QFileDialog.ShowDirsOnly)
        # added to display the chosen path
        if self.image_path:
            self.image_label.setText(self.image_path)

class MetricSelectPage(QWizardPage):
    butt_signal = Signal(str)  # Make a signal, pass it

    def __init__(self, parent=None):
        QWizardPage.__init__(self, parent)
        self.setTitle("Select metrics to run:")
        # https://stackoverflow.com/questions/4008649/qlistwidget-and-multiple-selectionllllll
        # https://www.geeksforgeeks.org/pyqt5-qlistwidget-setting-selection-mode/
        self.layout = QtWidgets.QVBoxLayout()
        self.listWidget = QtWidgets.QListWidget()
        self.listWidget.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
        self.listWidget.setGeometry(QtCore.QRect(10, 10, 211, 291))

        numUnboundCells = QtWidgets.QListWidgetItem("Number of Unbound Cells")
        numBoundCells = QtWidgets.QListWidgetItem("Number of Bound Cells")
        totalArea = QtWidgets.QListWidgetItem("Total Area")
        totalBoundedArea = QtWidgets.QListWidgetItem("Total Bounded Area")
        meanVorArea = QtWidgets.QListWidgetItem("Mean Voronoi Area")
        perc6SidesVor = QtWidgets.QListWidgetItem("Percent Six-Sided Voronoi")
        densityUn = QtWidgets.QListWidgetItem("Density (Uncorrected)")
        densityCor = QtWidgets.QListWidgetItem("Density (Corrected)")
        nndUn = QtWidgets.QListWidgetItem("Nearest Neighbor Distance (Uncorrected)")
        nndCor = QtWidgets.QListWidgetItem("Nearest Neighbor Distance (Corrected)")
        icdUn = QtWidgets.QListWidgetItem("Inter-Cell Distance (Uncorrected)")
        icdCor = QtWidgets.QListWidgetItem("Inter-Cell Distance (Corrected)")
        fndUn = QtWidgets.QListWidgetItem("Furthest Neighbor Distance (Uncorrected)")
        fndCor = QtWidgets.QListWidgetItem("Furthest Neighbor Distance (Corrected)")
        densityRecProDist = QtWidgets.QListWidgetItem("Density Recovery Profile Distance")
        vorAreaRegIndex = QtWidgets.QListWidgetItem("Voronoi Area Regularity Index")
        vorNumSidesRegIndex = QtWidgets.QListWidgetItem("Voronoi Number of Sides Regularity Index")
        nnRegIndex = QtWidgets.QListWidgetItem("Nearest Neighbor Regularity Index")
        icrIndex = QtWidgets.QListWidgetItem("Inter-Cell Regularity Index")

        self.listWidget.addItem(numUnboundCells)
        self.listWidget.addItem(numBoundCells)
        self.listWidget.addItem(totalArea)
        self.listWidget.addItem(totalBoundedArea)
        self.listWidget.addItem(meanVorArea)
        self.listWidget.addItem(perc6SidesVor)
        self.listWidget.addItem(densityUn)
        self.listWidget.addItem(densityCor)
        self.listWidget.addItem(nndUn)
        self.listWidget.addItem(nndCor)
        self.listWidget.addItem(icdUn)
        self.listWidget.addItem(icdCor)
        self.listWidget.addItem(fndUn)
        self.listWidget.addItem(fndCor)
        self.listWidget.addItem(densityRecProDist)
        self.listWidget.addItem(vorAreaRegIndex)
        self.listWidget.addItem(vorNumSidesRegIndex)
        self.listWidget.addItem(nnRegIndex)
        self.listWidget.addItem(icrIndex)

        self.listWidget.itemClicked.connect(self.printItemText)
        self.layout.addWidget(self.listWidget)
        self.setLayout(self.layout)

    def printItemText(self):
        items = self.listWidget.selectedItems()
        x = []
        for i in range(len(items)):
            x.append(str(self.listWidget.selectedItems()[i].text()))

        print(x)

class UnitSelect(QWizardPage):
    butt_signal = Signal(str)  # Make a signal, pass it

    def __init__(self, parent=None):
        QWizardPage.__init__(self, parent)
        self.setTitle("Select units:")

        self.layout = QtWidgets.QVBoxLayout()
        self.listWidget = QtWidgets.QListWidget()
        self.listWidget.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.listWidget.setGeometry(QtCore.QRect(10, 10, 211, 291))

        microns = QtWidgets.QListWidgetItem("Microns (mm density)")
        degrees = QtWidgets.QListWidgetItem("Degrees")
        arcmin = QtWidgets.QListWidgetItem("Arcmin")

        self.listWidget.addItem(microns)
        self.listWidget.addItem(degrees)
        self.listWidget.addItem(arcmin)

        self.listWidget.itemClicked.connect(self.printItemText)
        self.layout.addWidget(self.listWidget)
        self.setLayout(self.layout)

    def printItemText(self):
        items = self.listWidget.selectedItems()
        x = []
        for i in range(len(items)):
            x.append(str(self.listWidget.selectedItems()[i].text()))

        print(x)

class ResultsSaveLocation(QWizardPage):
    butt_signal = Signal(str)  # Make a signal, pass it

    def __init__(self, parent=None):
        QWizardPage.__init__(self, parent)
        self.setTitle("Select file location to save results to:")

        self.save_path = ""

        self.file_dir_form = QGridLayout()

        self.save_path_label = QLineEdit()
        self.save_path_label.setReadOnly(True)
        self.save_path_label.text()
        self.save_path_butt = QPushButton("Select...")

        self.file_dir_form.addWidget(self.save_path_label, 0, 0)
        self.file_dir_form.addWidget(self.save_path_butt, 0, 1)

        self.v_layout = QVBoxLayout()
        self.v_layout.setSpacing(32)
        self.v_layout.addLayout(self.file_dir_form)

        self.file_dir_form.setSpacing(4)

        self.setLayout(self.v_layout)

        self.save_path_butt.clicked.connect(self.select_save_path)



    @Slot()
    def select_save_path(self):
        self.save_path = QFileDialog.getExistingDirectory(parent=self,
                                                           caption="Specify file location to save results to.",
                                                           options=QFileDialog.ShowDirsOnly)
        # added to display the chosen path
        if self.save_path:
            self.save_path_label.setText(self.save_path)


class Calculate(QWizardPage):
    butt_signal = Signal(str)  # Make a signal, pass it

    def __init__(self, parent=None):
        QWizardPage.__init__(self, parent)
        self.setTitle("Calculate Metrics:")

        self.file_dir_form = QGridLayout()

        # https://www.google.com/search?q=qprogressbar&rlz=1C1GCEU_enUS949US949&oq=Qprogress&aqs=chrome.0.0i67j69i57j0i67l5j0i512l3.2106j0j7&sourceid=chrome&ie=UTF-8#kpvalbx=_YBpgYvOqLpXV9AOB2bqgDQ11
        n = 500  #steps for progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(n)

        self.start_butt = QPushButton("Start")
        self.file_dir_form.addWidget(self.progress_bar, 0, 0)
        self.file_dir_form.addWidget(self.start_butt, 0, 1)

        self.v_layout = QVBoxLayout()
        self.v_layout.setSpacing(32)
        self.v_layout.addLayout(self.file_dir_form)

        self.file_dir_form.setSpacing(4)

        self.setLayout(self.v_layout)
        # https://stackoverflow.com/questions/6784084/how-to-pass-arguments-to-functions-by-the-click-of-button-in-pyqt
        self.start_butt.clicked.connect(lambda: self.start(n))

    @Slot()
    def start(self, n):
        for i in range(n):
            time.sleep(0.01)
            self.progress_bar.setValue(i+1)


if __name__ == '__main__':

    app = QApplication([])
    widget = PygmyFeeder()
    widget.resize(512, 384)
    widget.show()

    sys.exit(app.exec())