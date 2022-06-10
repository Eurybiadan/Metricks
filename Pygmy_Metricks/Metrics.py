#  Copyright (c) $originalComment.match("Copyright \(c\) (\d+)", 1, "-")2022. Robert F Cooper
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
from os.path import exists
from warnings import warn

import cv2
import numpy
from PIL import Image

import ocvl.FeederGUI
import csv
import pandas

class Metricks():
    def __init__(self):
        super().__init__()

    # adapted from read_folder_contents.m
    def readFolderContents(self, directoryName, extension):
        self.directoryName = directoryName
        x=1
        self.fileList = []
        for file in os.listdir(directoryName):
            if file.endswith(extension):
                self.fileList.append(os.path.join(directoryName, file))
        self.numOfFiles = len(self.fileList)
        ocvl.FeederGUI.Calculate.n = self.numOfFiles

    def loadLutFile(self, path):
        self.LUT = pandas.read_csv(path[0], header=None)
        self.findScale()

    def findScale(self):
        # Find which row index matches our files
        for i in self.fileList:
            i = i.split('\\')
            i = i[1].split('_')
            self.LUT["idMatch"] = self.LUT[0].str.contains(i[1], regex=True)
            index = self.LUT.all(axis=1)
            for j in range(len(index)):
                if self.LUT["idMatch"][j] == True:
                    self.LUT["eyeMatch"] = self.LUT[3].str.contains(i[3], regex=True)
                    if self.LUT["eyeMatch"][j] == True:
                        self.LUTindex = j
                        break
            break  # not sure if this can really go here for all scenarios

        # Calculate scale
        self.axialLength = self.LUT[1][self.LUTindex]
        self.pixelsPerDegree = self.LUT[2][self.LUTindex]
        self.micronsPerDegree = (291*self.axialLength)/24

    def scaleInput(self, scaleInput):
        self.scaleInput = scaleInput
        self.scaleval = self.scaleInput

    def selectUnit(self, unit):
        self.selectedUnit = unit
        if self.selectedUnit == 'Microns (mm density)':
            self.scaleval = 1/(self.pixelsPerDegree/self.micronsPerDegree)
        if self.selectedUnit == 'Degrees':
            self.scaleval = 1/self.pixelsPerDegree
        if self.selectedUnit == 'Arcmin':
            self.scaleval = 60/self.pixelsPerDegree

    def runMetricks(self, i):
        self.coords = pandas.read_csv(self.fileList[i])
        # https://note.nkmk.me/en/python-pandas-len-shape-size/
        if len(self.coords.columns) != 2:
            # https://docs.python.org/3/library/warnings.html
            warn("Coordinate list contains more than 2 columns! Skipping...")
            return
        self.matchingName = self.fileList[i].split("_coords.csv")
        self.imageFilePath = self.matchingName[0] + '.tif'
        if exists(self.imageFilePath):
            self.image = Image.open(self.imageFilePath, 'r')
            self.imageValues = numpy.array(self.image)
            self.imageValues = pandas.DataFrame(self.imageValues)
            # self.width = len(self.image, 2)
            # self.height = len(self.image, 1)
