# Metricks - A MATLAB package for analyzing the cone photoreceptor mosaic.

Metricks is designed to take away the typically "ad-hoc" and unvalidated nature of much of the metrics extracted from Adaptive Optics Opthalmoscopy images.

### It performs three main functions:
* [Calculates metrics from a set of images, coordinates, and a scale file.](#Metrics) (Citation: Cooper RF, Wilk MA, Tarima S, Dubra A, Carroll J. “Evaluating descriptive metrics of the human cone mosaic.” Invest Ophthalmol Vis Sci. 2016 57(7):2993)
* [Creates maps of any metric we can calcuate in the above paper.](#MetricsMap)
* [Autmoatically estimates spacing and density in images with periodic structures.](#DFTMetrics) (Citation: Cooper RF, Aguirre GK, Morgan JIWM. "Fully-Automated Estimation of Cone Spacing and Density for Retinal Montages." Submitted.)
* [Calculates the local anisotropy or orientation of retinal structures.](#Orientation) (Citation: Cooper RF, Lombardo M, Carroll J, Sloan KR, Lombardo G. 2015 “Methods for investigating the local spatial anisotropy and the preferred orientation of cones in adaptive optics retinal images.” Visual Neurosci. 2016 33:E005)

## Information about how to calculate a set of metrics from images, coordinates and scale files: <a name="Metrics"></a>

### Run Coordinate_Mosiac_Metrics.m:

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

**If you do not wish to use a lookup table, then press "cancel", and the software will allow you put in your own scale in UNITS/pixel.**

**This software has the ability to pre-crop the input data (if, for example, you have 80 pixels of coordinates and you only want to analyze the middle 50).**

To specify a cropping window, input the size (in the units you are going to use) in to the brackets on line 10 of Coordinate_Mosaic_Metricks.m.

**Cropping is governed by the following rules:**

- If the tif is present and windowsize is not specified, the analysis will be done on everything within the dimensions of the image.
- If the tif is present and windowsize is specified, the assumed center of the image is calculated according to the borders of the tif. **In either case, it doesn’t “care” how many (or even if there are any) cells in the image.**
- If the tif is not present and windowsize is not specified, the analysis will be done on everything within the min and max coordinates in both x and y directions. So if you have an image in which there is an absence of cells on one side, for example, you might end up with a clipped area that is not a square.
- If the tif is not present and windowsize is specified, the assumed center of the image is calculated according to the min and max coordinates in both x and y directions. So if you have an image in which there is an absence of cells on one side, the center will shift towards the other side of the image.

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

## How to create a map of the above metrics: <a name="MetricsMap"></a>

### Run Coordinate_Mosiac_Metrics_MAP.m:

This script creates a map of an image/coordinate set across a set of image/coordinates. It shares a lot of the previous rules, so I will simply note the differences.

As before, to run this script, the script will prompt the user to select an folder containing image/coordinate pair. It will then create a map of each image/coordinate pair in a folder.

**However, unlike the above script, This software will automatically adjust its window size to encompass up 100 coordinates.**

If you wish to specify the sliding window size, input the size (in the units you are going to in to the brackets of the variable "WINDOW_SIZE" on line ~92 of Coordinate_Mosaic_Metrics_MAP.m. The window is governed by the same rules outlined previously.

## How to automatically assess the spacing/density of objects in an image or montage: <a name="DFTMetrics"></a>

### Run fit_fourier_spacing.m:

Running this function will ask the you to select the image that will be used for analysis. It will then calculate the DFT-derived spacing across the entire image, and return the average pixel spacing of the input image.

This script can also be run with arguments, and has the form:

`[avg_pixel_spac, interped_spac_map, interped_conf_map, sum_map, imbox ] = fit_fourier_spacing(test_image, roi_size, supersampling, row_or_cell)`

#### Inputs:

- **test_image**: The image that will be analyzed. The only requirement is that it is a 2d, grayscale (1 channel) image.
- **supersampling**: If "true", then each roi will be super-sampled in accordance with: [Bernstein et al.](https://arxiv.org/pdf/1401.2636.pdf) before calculating the DFT-derived spacing.
- **roi_size**: The side length (in pixels) of a sliding roi window- The roi will march along the image you've provided at a rate of 1/4 the size of the ROI, creating a "map" of spacing of the image.
- **row_or_cell**: The range of angles from the polar DFT that will be used to calculate the DFT-derived spacing. If "row", then it will be the upper and lower 90 degrees of the DFT. If "cell", it will be the left and right 90 degrees.

#### Outputs:

- **avg_pixel_spac**: The average spacing of the image.
- **interped_spac_map**: The spacing map of the input image (in pixel spacing).
- **interped_conf_map**: The confidence map of the input image.
- **sum_map**: The map corresponding to the amount of ROI overlap across the output map.
- **imbox**: The bounding region of valid (nonzero, NaN, or Inf) pixels.

### Run Montage_DFT_Analysis.m:

Running this function will ask you to select the set of montage images that you wish to analyze. Montage images will have a single layer within a montage placed within an image the size of the montage canvas.

**If you only have a PSD of a montage and do not have data in this format, please use [my PSD Layer Exporter](https://github.com/Eurybiadan/PSD_Layer_Export/releases/tag/v1.0). To dump them to disk in this format.

**This function uses parfor loops, which requires the Parallel Toolbox from MATLAB.** If you don't have that toolbox, change the "parfor" on line 93 to a "for" loop.

As above, it will then prompt the user to select what the output unit should be. At present, the options are:
* Microns (using millimeters^2 for density)
* Degrees
* Arcminutes

Once the output unit is select, it will give the user the option to pick a lookup table. The lookup table allows the software to analyze a folder of images from different subjects/timepoints/conditions. The lookup table itself **must** be a 3 column 'csv' file, where the **first column** is a common identifier for image/coordinate pairs, the **second column** is the axial length (or '24' if the axial length is unknown) of the image/coordinate pairs, and the **third column** is the pixels per degree of the image/coordinate pairs. Each row must contain a different identifier/axial length/pixels per degree tuple.

An example common identifier could be a subject number, e.g, when working with the files
- 1235_dateoftheyear_OD_0004.tif
- 1235_dateoftheyear_OD_0005.tif

Common identifiers could be "1235", "1235_dateoftheyear", "1235_dateoftheyear_OD". If all three were placed in a LUT, then the one that matches the most (as determined via levenshtein distance) will be used. In this example, we would use "1235_dateoftheyear_OD".

If we had another date, say: 1235_differentdateoftheyear_OD_0005.tif, then _only_ the identifier "1235" would match between all images. However, say the two dates have different scales, then you would want to create two rows in the look up table for each date, with identifiers like: "1235_dateoftheyear" and "1235_differentdateoftheyear".

**If you do not wish to use a lookup table, then press "cancel", and the software will allow you put in your own scale in UNITS/pixel.**

The software will then run, showing the spacing montage, the density montage, the confidence montage, and the sum map. It will also save them to disk in the same folder you ran from alongside a mat file that contains the results.

## How to calculate the orientation of cells in an image: <a name="Orientation"></a>

Ask me via email or via issue- this code is older and less user-friendly.

# Don't thank me; cite me:
Every metric that is run via the main "Coordinate_Mosaic_Metrics.m" script has been validated and used in the following manuscript: **Cooper RF, Wilk MA, Tarima S, Dubra A, Carroll J. “Evaluating descriptive metrics of the human cone mosaic.” Invest Ophthalmol Vis Sci. 2016 57(7):2993.** You can also find formal definitions of each metric calculated here in that paper.

**This package is free for use under GPL v3, but I ask that you please cite the above paper if you use this package.**

