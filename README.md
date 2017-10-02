# Metricks - A MATLAB package for analyzing the cone photoreceptor mosaic.

Metricks is designed to take away the typically "ad-hoc" and unvalidated nature of much of the metrics extracted from Adaptive Optics Opthalmoscopy images.

## Information about running Coordinate_Mosaic_Metrics.m:

When run, the script will prompt the user to select a folder with image/coordinate pairs.

**At present, images must be 8-bit grayscale tifs, coordinates must be formatted as a 2 column matrix (x,y), and must be named using the following convention, where [imagename] can be any valid filename:**
* Image File: [imagename].tif
* Coordinate File: [imagename]\_coords.csv

It will then prompt the user to select what the output unit should be. At present, the options are:
* Microns (using millimeters^2 for density)
* Degrees
* Arcminutes

Once the output unit is select, it will give the user the option to pick a lookup table. The lookup table allows the software to analyze a folder of images from different subjects/timepoints/conditions. The lookup table itself **must** be a 3 column 'csv' file, where the **first column** is a common identifier for image/coordinate pairs, the **second column** is the axial length (or '24' if the axial length is unknown) of the image/coordinate pairs, and the **third column** is the pixels per degree of the image/coordinate pairs. Each row must contain a different identifier/axial length/pixels per degree tuple.

An example common identifier could be a subject number, e.g, when working with the files
- 1235_dateoftheyear_OD_0004.tif
- 1235_dateoftheyear_OD_0005.tif

Common identifiers could be "1235", "1235_dateoftheyear", "1235_dateoftheyear_OD". If all three were placed in a LUT, then the one that matches the most (as determined via levenshtein distance) will be used. In this example, we would use "1235_dateoftheyear_OD".

If we had another date, say: 1235_differentdateoftheyear_OD_0005.tif, then _only_ the identifier "1235" would match between all images. However, say the two dates have different scales, then you would want to create two rows in the look up table for each date, with identifiers like: "1235_dateoftheyear" and "1235_differentdateoftheyear".

**If you do not wish to use a lookup table, then press "cancel", and the software will allow you put in your own scale in UNITS/pixel**

The software will then run, and calculate every metric currently validated.

At present, it calculates the following metrics from each image and coordinate pair:

- Number of Unbound Cells
- Number of Bound Cells
- Total Area
- Total Bounded Area
- Mean Voronoi Area
- Percent Six-Sided Voronoi
- Density (uncorrected/corrected)
- Nearest Neighbor Distance (uncorrected/corrected)
- Inter-Cell Distance (uncorrected/corrected)
- Furthest Neighbor Distance (uncorrected/corrected)
- Density Recovery Profile Distance
- Voronoi Area Regularity Index
- Voronoi Number of Sides Regularity Index
- Nearest Neighbor Regularity Index
- Inter-Cell Regularity Index

The results will then be placed in to a datestamped file within a "Results" folder as a subfolder of the one selected for analysis.

Every metric that is run via the main "Coordinate_Mosaic_Metrics.m" script has been validated and used in the following manuscript: **Cooper RF, Wilk MA, Tarima S, Dubra A, Carroll J. “Evaluating descriptive metrics of the human cone mosaic.” Invest Ophthalmol Vis Sci. 2016 57(7):2993.** You can also find formal definitions of each metric calculated here in that paper.

**This package is free for use under GPL v3, but I ask that you please cite the above paper if you use this package.**


