
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

import sys
from PySide6.QtWidgets import QApplication, QWizard, QWizardPage, QVBoxLayout, QGridLayout, QPushButton, QLabel, \
                              QLineEdit, QFileDialog
from PySide6.QtCore import Qt, Signal, Slot, Property

class PygmyFeeder(QWizard):
    def __init__(self, parent=None):
        QWizard.__init__(self,parent)
        self.setWizardStyle(QWizard.ModernStyle)  # added to get the title to be formatted correctly
        self.setWindowTitle("Welcome to OCVL's Metricks Master (Pygmy Python editon)")
        self.addPage(WelcomePage(self))


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

        self.butt_signal
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


if __name__ == '__main__':

    app = QApplication([])
    widget = PygmyFeeder()
    widget.resize(512, 384)
    widget.show()

    sys.exit(app.exec())