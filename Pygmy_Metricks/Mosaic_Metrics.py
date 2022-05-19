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
import os
import sys
import importlib.resources as pkg_resources
from pint import UnitRegistry


from PySide6.QtWidgets import QApplication

from ocvl.FeederGUI import PygmyFeeder


class PygmyMetricks():
    def __init__(self):
        super().__init__()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    ureg = UnitRegistry()
    with pkg_resources.path("ocvl", "ocvl-pint-defs.txt") as pint_path:
        deffile = open(pint_path, "r")
        ureg.load_definitions(deffile)

    app = QApplication([])
    widget = PygmyFeeder()
    widget.resize(800, 600)
    widget.show()

    sys.exit(app.exec())
