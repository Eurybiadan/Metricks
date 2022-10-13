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
import math
import os
from operator import xor
from os.path import exists
from warnings import warn

import cv2
import numpy
import scipy.spatial
import self
from PIL import Image
from numpy import mean, std, matlib
from numpy.ma import size
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.spatial import voronoi_plot_2d

import ocvl.FeederGUI
import csv
import pandas
from scipy.spatial.distance import squareform, pdist, cdist


class Metricks():
    def __init__(self):
        super().__init__()

    # adapted from read_folder_contents.m
    def readFolderContents(self, directoryName, extension):
        self.directoryName = directoryName
        x = 1
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
        self.micronsPerDegree = (291 * self.axialLength) / 24

    def scaleInput(self, scaleInput):
        self.scaleInput = scaleInput
        self.scaleval = self.scaleInput

    def selectUnit(self, unit):
        self.selectedUnit = unit
        if self.selectedUnit == 'Microns (mm density)':
            self.scaleval = 1 / (self.pixelsPerDegree / self.micronsPerDegree)
        if self.selectedUnit == 'Degrees':
            self.scaleval = 1 / self.pixelsPerDegree
        if self.selectedUnit == 'Arcmin':
            self.scaleval = 60 / self.pixelsPerDegree

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
                self.pixelWindowSze = windowSize / self.scaleval
                self.diffWidth = (self.width - self.pixelWindowSze) / 2
                self.diffHeight = (self.height - self.pixelWindowSze) / 2

                if self.diffWidth < 0:
                    self.diffWidth = 0
                if self.diffHeight < 0:
                    self.diffHeight = 0
            else:
                self.pixelWindowSze = [self.height, self.width]
                self.diffWidth = 0
                self.diffHeight = 0

            self.clippedCoords = self.coordClip(self.coords, [self.diffWidth, self.width - self.diffWidth],
                                                [self.diffHeight, self.height - self.diffHeight], 'i')
            self.clipStartEnd = [self.diffWidth, self.width - self.diffWidth, self.diffHeight,
                                 self.height - self.diffHeight]


        else:
            self.width = max(self.coords.iloc[:, 0]) - min(self.coords.iloc[:, 0])
            self.height = max(self.coords.iloc[:, 1]) - min(self.coords.iloc[:, 1])

            if len(windowSize) != 0:
                self.pixelWindowSze = windowSize / self.scaleval
                self.diffWidth = (self.width - self.pixelWindowSze) / 2
                self.diffHeight = (self.height - self.pixelWindowSze) / 2
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
        self.statistics = self.determineMosaicStats(self.clippedCoords, self.scaleval, self.selectedUnit,
                                                    self.clipStartEnd, [self.pixelWindowSze, self.pixelWindowSze], 4)
        # HERE JG 6/10/22

    def determineMosaicStats(self, coords, scale, unit, bounds, clippedRowCol, reliability):
        # This function takes in a list of coordinates in a m-2 matrix and calculates metrics
        # calculates the mean nearest neighbor, cell area created by the coordinates, and calculates the density of the coordinates
        # Coords are in X,Y

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Determine Mean N-N
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # Measure the distance from each set of points to the other
        # https://stackoverflow.com/questions/32946241/scipy-pdist-on-a-pandas-dataframe
        distBetweenPts = pdist(coords, 'euclidean')
        squareform(distBetweenPts)
        distBetweenPts = pandas.DataFrame(squareform(distBetweenPts), index=coords.index,
                                          columns=coords.index)  # Measure the distance from each set of points to the other

        # Make diagonal not the minimum for any observation
        maxDist = max(distBetweenPts.max())
        maxIdent = numpy.eye(len(distBetweenPts)) * maxDist

        # Find the minimum distance from one set of obs to another
        combined = distBetweenPts.add(maxIdent, fill_value=0)
        minval = combined.min()
        meanNNDist = mean(minval * scale)  # Distance in units

        scaledMinVal = minval * scale
        stdScaledMinVal = std(scaledMinVal,
                              ddof=1)  # https://stackoverflow.com/questions/27600207/why-does-numpy-std-give-a-different-result-to-matlab-std
        regularityNNIndex = meanNNDist / stdScaledMinVal

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Determine Voronoi Cell Area
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        sixSided = 0
        bound = pandas.DataFrame(numpy.zeros((size(coords, 1), 1)))  #https://www.adamsmith.haus/python/answers/how-to-create-a-zero-filled-pandas-dataframe-in-python#:~:text=Use%20numpy.,size%20shape%20populated%20with%200.0%20.
        cellArea = pandas.DataFrame()
        numEdges = pandas.DataFrame(numpy.zeros((size(coords, 1), 1)))
        coordsBound = pandas.DataFrame()

        if size(coords, 1) > 2:
            points = scipy.spatial.Voronoi(coords, qhull_options='QJ')
            # fig = voronoi_plot_2d(points)
            # plt.show()
            V = pandas.DataFrame(points.vertices)
            # V.to_excel("C:\\Users\\6794grieshj\\Documents\\help.xlsx")
            C = pandas.DataFrame(points.regions)
            # C.to_excel("C:\\Users\\6794grieshj\\Documents\\help2.xlsx")

            for i in range(len(C)):
                cTemp = C.iloc[i]
                cTemp = [item for item in cTemp if not(math.isnan(item))]
                vertices = V.iloc[cTemp, :]
                x = all(cTemp)
                a = vertices.iloc[:, 0] < bounds[1]
                b = vertices.iloc[:, 1] < bounds[3]
                c = vertices.iloc[:, 0] > bounds[0]
                d = vertices.iloc[:, 1] > bounds[2]
                a1 = all(a)
                b1 = all(b)
                c1 = all(c)
                d1 = all(d)
                if x & a1 & b1 & c1 & d1:  # [xmin xmax ymin ymax]
                    cellArea.loc[i, 0] = self.PolyArea(V.iloc[cTemp, 0], V.iloc[cTemp, 1])
                    numEdges.iloc[i] = len(V.iloc[cTemp, 0])
                    # print(numEdges.iloc[i, 0])
                    if numEdges.iloc[i, 0] == 6:
                        sixSided = sixSided + 1
                    # print(coords.iloc[i])
                    coordsBound.loc[i, 0] = coords.iloc[i, 0]
                    coordsBound.loc[i, 1] = coords.iloc[i, 1]
                    bound.iloc[i] = 1

        if len(coordsBound) != 0:
            cellArea = cellArea * (scale ** 2)  # convert to square microns
            numEdges = numEdges.loc[(numEdges != 0).any(axis=1)]  # removes the zeros; https://stackoverflow.com/questions/22649693/drop-rows-with-all-zeros-in-pandas-data-frame

            meanCellArea = mean(cellArea)
            regularityVoroIndex = meanCellArea/std(cellArea, ddof=1)
            regularityVoroSides = mean(numEdges)/std(numEdges, ddof=1)
            print(len(coordsBound))
            # coordsBound.to_excel("C:\\Users\\6794grieshj\\Documents\\coordBound.xlsx")  # there is one extra bound cell in python
            percentSixSided = 100*sixSided/len(coordsBound)
        else:
            cellArea = 0
            meanCellArea = 0
            regularityVoroIndex = 0
            regularityVoroSides = 0
            percentSixSided = 0

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Determine number of cells, density direct count (D_dc)
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        numCells = len(coords)  # Total number of cells
        totalCellArea = cellArea.sum()

        if unit == "Microns (mm density)":
            self.totalCoordArea = ((clippedRowCol[0][0] * clippedRowCol[0][1]) * ((scale ** 2) / (1000 ** 2)))

        else:
            self.totalCoordArea = ((clippedRowCol[0][0] * clippedRowCol[0][1]) * (scale ** 2))

        pixelDensity = numCells / (clippedRowCol[0][0] * clippedRowCol[0][1])
        densityDc = numCells / self.totalCoordArea

        if len(coordsBound) != 0:  # check if not empty
            if unit == "Microns (mm density)":
                densityBound = (1000 ** 2) * len(coordsBound)/totalCellArea
            else:
                densityBound = len(coordsBound)/totalCellArea
        else:
            densityBound = 0


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Determine Inter-Cell Distance
        m = 1
        interCellDist = list()
        maxCellDist = list()

        correctInterCellDist = pandas.DataFrame(numpy.zeros((size(coords, 1), 1)))
        correctMaxCellDist = pandas.DataFrame(numpy.zeros((size(coords, 1), 1)))
        correctNNCellDist = pandas.DataFrame(numpy.zeros((size(coords, 1), 1)))

        if len(coords) > 2:
            dt = scipy.spatial.Delaunay(coords)

            # Find all instances of each coordinate point
            for k in range(0, len(coords)):
                i = numpy.where(dt.simplices == k)
                i = i[0]
                connInd = dt.simplices[i, :]
                # connInd = connInd[0]
                coordRow = numpy.unique(connInd[connInd != k])

                if size(i, 0) != 1:
                    coordRow = numpy.insert(coordRow, 0, k)
                else:
                    # https://stackoverflow.com/questions/39885495/what-is-the-meaning-of-single-quote-in-matlab-and-how-to-change-it-to-python
                    coordRow = numpy.insert(coordRow, 0, k)
                    coordRow = coordRow.T

                temp1 = coords.iloc[coordRow, 0]
                temp2 = coords.iloc[coordRow, 1]
                temp3 = numpy.empty((0, 2), int)
                for z in range(len(temp1)):
                    temp3 = numpy.append(temp3, numpy.array([[temp1.iloc[z], temp2.iloc[z]]]), axis=0)
                cellDist = squareform(pdist(temp3))

                if bound.iloc[k][0] == 1:
                    # if it is bound, then we've flagged it as such, and can use it in the triangulation
                    # only take the first row because that is the cell of interest's relative distance to its neighboring cells
                    correctInterCellDist[0][m] = scale * (sum(cellDist[0, :]) / (len(cellDist[0, :])-1))
                    correctMaxCellDist[0][m] = scale * max(cellDist[0, :])
                    correctNNCellDist[0][m] = scale * min(cellDist[0, 1:])
                    m = m + 1

                interCellDist.append(scale * (sum(cellDist[0, :]) / (len(cellDist[0, :])-1)))
                maxCellDist.append(scale*max(cellDist[0, :]))

            m = m-1

            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! All these distances need to be verified !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            meanInterCellDist = mean(interCellDist)
            meanMaxCellDist = mean(maxCellDist)

        else:
            meanInterCellDist = scale * pdist(coords)
            meanMaxCellDist = meanInterCellDist

        if len(coordsBound != 0):
            meanCorrectNNDist = mean(correctNNCellDist[0:m])
            meanCoorectICDist = mean(correctInterCellDist[0:m])
            regularityICIndex = mean(correctInterCellDist[0:m] / std(correctInterCellDist[0:m]))
            meanCorrectMaxCellDist = mean(correctMaxCellDist[0:m])
        else:
            regularityICIndex = 0
            meanCorrectNNDist = 0
            meanCoorectICDist = 0
            meanCorrectMaxCellDist = 0


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Determine Density Recovery Profile

        coordBounds = (bounds[0:1], bounds[2:3])
        densityDc = pixelDensity

        width = coordBounds[0][1] - coordBounds[0][0]
        height = coordBounds[1][1] - coordBounds[1][0]

        # Equation for determine min bin size Rodieck (1990)
        s = min([width, height])  # smallest length of the two sides
        k = reliability  # reliability factor
        D = densityDc  # Density in cones/pixels^2

        pixPerBin = round(k / (s * D * math.sqrt(math.pi)))
        # If it's at a subpixel level, cap it to be a single pixel and a single micron/pixel
        if pixPerBin < 0.5:
            pixPerBin = 1

        numOfBins = round((s/pixPerBin)/2)
        if numOfBins < 5:
            numOfBins = 5

        drpSizes = pandas.DataFrame(numpy.zeros((size(numOfBins, 1), 1)))
        drpSizes[0] = 1  # to prevent dividing by zero later

        for i in range(2,numOfBins-1):
            drpSizes[i] = pixPerBin * (i-1)

        # Make sure we don't have duplicate sizes
        drpSizes = numpy.unique(drpSizes)
        scaledDrpSizes = drpSizes * scale

        numCellsInAnnulus = pandas.DataFrame(numpy.zeros((size(coords, 1), len(drpSizes))))
        densityPerAnnulus = pandas.DataFrame(numpy.zeros((len(drpSizes-1), 1)))
        edgeFactor = 1-1/math.pi  # Edge factor integrated

        for i in range(0,len(drpSizes), 1):
            xMin = coordBounds[0][0] + drpSizes[i]
            xMax = coordBounds[0][1] - drpSizes[i]

            yMin = coordBounds[1][0] + drpSizes[i]
            yMax = coordBounds[1][1] - drpSizes[i]

            xBounds = [xMin, xMax]
            yBounds = [yMin, yMax]

            # Establish which cells are inside and outside hte region we'll use

            cellCoordsXOR = self.coordClip(coords, xBounds, yBounds, 'xor')
            cellCoordsInside = self.coordClip(coords, xBounds, yBounds, 'i')
            cellCoordsAND = self.coordClip(coords, xBounds, yBounds, 'and')

            numCellsInside = size(cellCoordsInside, 0)
            numCellsXOR = size(cellCoordsXOR, 0)
            numCellsAND = size(cellCoordsAND, 0)

            coordsReordered = [cellCoordsInside, cellCoordsXOR, cellCoordsAND]

            # https://stackoverflow.com/questions/43650931/python-alternative-for-calculating-pairwise-distance-between-two-sets-of-2d-poin
            # Take reordered coordinates and find the distance between them all - each coord is along the row
            umDistBetweenPts = cdist(coordsReordered, coordsReordered) * scale
            unadjustedArea = math.pi * (drpSizes[i] * drpSizes[i] - drpSizes[i-1] * drpSizes[i-1])  # In pixels

            # factor calculated for each
            edgeArea = (unadjustedArea * edgeFactor * (scale * scale) / (1000 * 1000))

            cornerFactor = 1 - (((2 * drpSizes[i]) / (math.pi * width * height)) * (width + height)) + ((drpSizes[i] * drpSizes[i]) / (math.pi * width * height))
            cornerArea = (unadjustedArea * cornerFactor * (scale * scale) / (1000 * 1000))

            unadjustedArea = unadjustedArea * (scale * scale) / (1000 * 1000)

            numCellsInAnnulus[:, i-1] = sum((scaledDrpSizes[i-1] < umDistBetweenPts) & (umDistBetweenPts <= scaledDrpSizes[i]), 2)
            centerDens = numCellsInAnnulus[1:numCellsInside, i-1] / unadjustedArea

            edgeDens = numCellsInAnnulus[numCellsInside+1: numCellsInside + numCellsXOR, i-1] / edgeArea

            cornerDens = numCellsInAnnulus[numCellsInside+numCellsXOR+1: numCellsInside+numCellsXOR+numCellsAND, i-1] / cornerArea

            densityPerAnnulus[i-1] = mean([centerDens, edgeDens, cornerDens])

        scaledDrpSizes = scaledDrpSizes[1:]  # add transpose here? - jg

        # Auto-peak finding approach
        change = numpy.diff(scaledDrpSizes)
        resampleX = change/10
        minSample = min(scaledDrpSizes)
        maxSample = max(scaledDrpSizes)
        interPDrpX = range(minSample, resampleX, maxSample)

        # If we don't have enough to fit splines, then just pick the 2nd value
        if (len(scaledDrpSizes) >= 2) and (len(densityPerAnnulus) >= 2) and (len(interPDrpX) >= 2):
            # https://stackoverflow.com/questions/20694809/matlab-spline-function-in-python
            splined = UnivariateSpline(scaledDrpSizes, densityPerAnnulus, interPDrpX)

            # Start looking for points AFTER the known "0's"
            lastZero = 1
            for i in range(0, len(splined)):
                if splined[i] == 0:
                    lastZero = i
                else:
                    break

            localExtremeaBins = interPDrpX[lastZero:]
            localExtrema = splined[lastZero:]

            if len(localExtrema) >= 3:
                self.extrema(localExtrema)

                # Only use the first peak that is higher than the mean height (designed to kill off the ER-level peaks)
                # https://stackoverflow.com/questions/28512237/python-equivalent-to-matlab-a-b-sorty
                self.maxesInds = self.maxesInds.sort(axis=1)
                sortId = self.maxesInds.argsort(axis=1)
                self.localMaxY = self.localMaxY[sortId]
                localMaxX = localExtremeaBins[self.maxesInds]

                meanDensity = mean(densityPerAnnulus)
                for i in range(0,len(self.localMaxY)):
                    if self.localMaxY[i] >= meanDensity:
                        localMaxX = localMaxX[i]
                        self.localMaxY = self.localMaxY[i]
                        break
                if len(localMaxX) != 0:
                    estSpacing = localMaxX
                else:
                    # If there aren't any peaks, pick the first peak that is above the mean density
                    for i in range(0, len(localExtrema)):
                        if localExtrema[i] >= meanDensity:
                            estSpacing = localExtremeaBins[i]
                            self.localMaxY = localExtrema[i]
                            break

            else:
                estSpacing = localExtrema[0]
        else:
            estSpacing = scaledDrpSizes[1]

        #PLOT something

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Output List Formatting
    # Make the returned struct
    if self.selectedUnit == 'microns(mm density)':
        self.totalCoordArea = self.totalCoordArea * 1000 ^ 2

    # mosaicStats = ...


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
            boolKey = ((coords.iloc[:, 1] > minYthresh) & (coords.iloc[:, 1] < maxYthresh) & (
                        coords.iloc[:, 0] > minXthresh) & (coords.iloc[:, 0] < maxXthresh))

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
                clippedCoords = clippedCoords.drop(
                    index=i)  # delete the row from the coordinate list if it was false from the boolean key

        return clippedCoords
    def extrema(self, x):
        """
        :param x:
        :return:
        """
        xMax = []
        iMax = []
        xMin = []
        iMin = []

        # Vector Input?
        Nt = numpy.size(x)
        if Nt != len(x):
            print("Entry must be a vector")

        # NaNs
        # https://www.mathworks.com/matlabcentral/answers/466615-what-is-an-equivalent-of-find-in-python
        iNan = numpy.argwhere(numpy.isnan(x))
        indX = range(0, Nt)
        if len(iNan) != 0:
            indX[iNan] = []
            x[iNan] = []
            Nt = len(x)

        # Difference between subsequent elements:
        dx = numpy.diff(x)

        # Is a horizontal line?
        if ~numpy.any(dx):
            return

        # Flat peaks? Put the middle element:
        a = numpy.argwhere(dx != 0)                     # Indexdes where x changes
        lm = numpy.argwhere(numpy.diff(a) != 1) + 1     # Indexes where a do not changes
        d = a[lm] - a[lm-1]                             # Number of elements in the flat peak
        a[lm] = a[lm] - math.floor(d/2)                 # Save middle element
        a[-1+1] = Nt

        # Peaks?
        xa = x[a]                               # Serie without flat peaks
        b = (numpy.diff(xa) > 0)                # 1 => positive slopes (mimima begin)
                                                # 0  =>  negative slopes (maxima begin)
        xb = numpy.diff(b)                      # -1 =>  maxima indexes (but one)
                                                # +1 =>  minima indexes (but one)

        iMax = numpy.argwhere(xb == -1) + 1     # Maxima indexes
        iMin = numpy.argwhere(xb == +1) + 1     # Minima indexes
        iMax = a[iMax]
        iMin = a[iMin]

        nMaxi = len(iMax)
        nMini = len(iMin)

        # Maximum or minimum on a flat peak at the end?
        if (nMaxi == 0) and (nMini == 0):
            if x[0] > x[Nt]:
                xMax = x[0]
                iMax = indX[0]
                xMin = x[Nt]
                iMin = indX[Nt]
            elif x[0] < x[Nt]:
                xMax = x[Nt]
                iMax = indX[Nt]
                xMin = x[0]
                iMin = indX[0]
            return

        # Maximum or minimum at the ends?
        if (nMaxi == 0):
            iMax[0, 1] = [0, Nt]
        elif (nMini == 0):
            iMin[0, 2] = [0, Nt]
        else:
            if iMax[0] < iMin[0]:
                iMin[1, nMini+1] = iMin
                iMin[0] = 1
            else:
                iMax[1, nMaxi+1] = iMin
                iMax[0] = 1
            if iMax[-1] > iMin[-1]:
                iMin[-1+1] = Nt
            else:
                iMax[-1+1] = Nt

        xMax = x[iMax]
        xMin = x[iMin]

        # NaN's
        if len(iNan) != 0:
            iMax = indX[iMax]
            iMin = indX[iMin]

        # Same size as x:
        # https://stackoverflow.com/questions/11892358/matlab-vs-python-reshape
        iMax = numpy.reshape(iMax, len(xMax))
        iMin = numpy.reshape(iMin, len(xMin))

        # Descending order:
        temp = -iMax.sort(axis=1)
        inMax = -iMax.argsort(axis=1)
        xMax = xMax[inMax]
        iMax = iMax[inMax]
        xMin = xMin.sort(axis=1)
        inMin = xMin.argsort(axis=1)
        iMin = iMin[inMin]


        # self.localMaxY = None
        # self.maxesInds = None


    def PolyArea(self, x, y):
        # https: // stackoverflow.com / questions / 24467972 / calculate - area - of - polygon - given - x - y - coordinates
        return 0.5 * numpy.abs(numpy.dot(x, numpy.roll(y, 1)) - numpy.dot(y, numpy.roll(x, 1)))
