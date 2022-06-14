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
from operator import xor
from os.path import exists
from warnings import warn

import cv2
import numpy
from PIL import Image
from numpy import mean, std, matlib

import ocvl.FeederGUI
import csv
import pandas
from scipy.spatial.distance import squareform, pdist


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

    def runMetricks(self, i, windowSize):
        self.coords = pandas.read_csv(self.fileList[i], header=None)
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
            # https://datatofish.com/numpy-array-to-pandas-dataframe/
            self.imageValues = pandas.DataFrame(self.imageValues)
            self.width = len(self.imageValues.columns)
            self.height = len(self.imageValues)
            print(len(windowSize))

            if len(windowSize) != 0:
                self.pixelWindowSze = windowSize/self.scaleval
                self.diffWidth = (self.width-self.pixelWindowSze)/2
                self.diffHeight = (self.height-self.pixelWindowSze)/2

                if self.diffWidth < 0:
                    self.diffWidth = 0
                if self.diffHeight < 0:
                    self.diffHeight = 0
            else:
                self.pixelWindowSze = [self.height, self.width]
                self.diffWidth = 0
                self.diffHeight = 0

            self.clippedCoords = self.coordClip(self.coords, [self.diffWidth, self.width-self.diffWidth], [self.diffHeight, self.height-self.diffHeight], 'i')
            self.clipStartEnd = [self.diffWidth, self.width-self.diffWidth, self.diffHeight, self.height-self.diffHeight]


        else:
            self.width = max(self.coords.iloc[:, 0]) - min(self.coords.iloc[:, 0])
            self.height = max(self.coords.iloc[:, 1]) - min(self.coords.iloc[:, 1])

            if len(windowSize) != 0:
                self.pixelWindowSze = windowSize/self.scaleval
                self.diffWidth = (self.width-self.pixelWindowSze)/2
                self.diffHeight = (self.height-self.pixelWindowSze)/2
            else:
                self.pixelWindowSze = [self.height, self.width]
                self.diffWidth = 0
                self.diffHeight = 0

            self.clippedCoords = self.coordClip(self.coords, [min(self.coords.iloc[:, 0]) - 0.01 + self.diffWidth,
                                                              max(self.coords.iloc[:, 0]) - self.diffWidth + 0.01],
                                                             [min(self.coords.iloc[:, 1]) + self.diffHeight - 0.01,
                                                              max(self.coords.iloc[:, 1]) - self.diffHeight + 0.01], 'i')
            # HERE JG 6/10/22
            self.clipStartEnd = [min(self.coords.iloc[:, 0]) + self.diffWidth - 0.01,
                                 max(self.coords.iloc[:, 0]) - self.diffWidth + 0.01,
                                 min(self.coords.iloc[:, 1]) + self.diffHeight - 0.01,
                                 max(self.coords.iloc[:, 1]) - self.diffHeight + 0.01]

        # NEED TO ADD IN THIS FUNCTION CALL AND FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        self.statistics = self.determineMosaicStats(self.clippedCoords, self.scaleval, self.selectedUnit, self.clipStartEnd, [self.pixelWindowSze, self.pixelWindowSze], 4)
        # HERE JG 6/10/22
    def coordClip(self, coords, thresholdx, thresholdy, inoutorxor):
        # determine max image size - assumes there are border coordinates
        imsizecol = max(coords.iloc[:, 0])
        imsizerow = max(coords.iloc[:, 1])

        # making threshold variables
        minYthresh = thresholdy[0]
        maxYthresh = thresholdy[1]
        minXthresh = thresholdx[0]
        maxXthresh = thresholdx[1]

        clippedCoords = coords

        if inoutorxor == 'i':
            # ensures that all coordinates are inside the box. - In accordance with notebook decision
            boolKey = ((coords.iloc[:, 1] > minYthresh) & (coords.iloc[:, 1] < maxYthresh) & (coords.iloc[:, 0] > minXthresh) & (coords.iloc[:, 0] < maxXthresh))

        elif inoutorxor == 'o':
            boolKey = ((coords.iloc[:, 1] > minYthresh) | (coords.iloc[:, 1] < maxYthresh) |
                       (coords.iloc[:, 0] > minXthresh) | (coords.iloc[:, 0] < maxXthresh))

        elif inoutorxor == 'xor':
            # Check rows coordinates for includable entries - in accordance with notebook decision
            boolKey = xor((coords.iloc[:, 1] > minYthresh) | (coords.iloc[:, 1] < maxYthresh) |
                          (coords.iloc[:, 0] > minXthresh) | (coords.iloc[:, 0] < maxXthresh))

        elif inoutorxor == 'and':
            # Check rows coordinates for includable entries - in accordance with notebook decision
            boolKey = (((coords.iloc[:, 1] > minYthresh) | (coords.iloc[:, 1] < maxYthresh)) &
                       ((coords.iloc[:, 0] > minXthresh) | (coords.iloc[:, 0] < maxXthresh)))

        else:
            return None

        for i in range(len(coords)):
            if not boolKey[i]:
                clippedCoords = clippedCoords.drop(index=i)  # delete the row from the coordinate list if it was false from the boolean key

        return clippedCoords

    def determineMosaicStats(self, coords, scale, unit, bounds, clippedRowCol, reliability):
        # This function takes in a list of coordinates in a m-2 matrix and calculates metrics
        # calculates the mean nearest neighbor, cell area created by the coordinates, and calculates the density of the coordinates
        # Coords are in X,Y

        # Determine Mean N-N
        # Measure the distance from each set of points to the other
        # https://stackoverflow.com/questions/32946241/scipy-pdist-on-a-pandas-dataframe
        distBetweenPts = pdist(coords, 'euclidean')
        squareform(distBetweenPts)
        distBetweenPts = pandas.DataFrame(squareform(distBetweenPts), index=coords.index, columns=coords.index)
        # distBetweenPts.to_excel("C:\\Users\\6794grieshj\\Documents\\help.xlsx")

        maxDist = max(distBetweenPts.max())
        maxIdent = numpy.eye(len(distBetweenPts)) * maxDist

        ##

        # [minval, minind] = min(distBetweenPts + maxIdent)
        # meanNNDist = mean(minval * scale)
        # regularityNNIndex = meanNNDist/std(minval*scale)




