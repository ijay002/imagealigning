README

Step 1: unzip lib.zip and copy the lib folder to the main IDL directory
step 2:	download and unzip align_sparkmass_github directory
	- unzip 2019_05_10_H_.zip into the same directory (this is the rendered DNA-PAINT image file required for image alignment
	- open script align.pro with idlde.exe
	- update string defined as dir in line 1 to the current path containing the files
	- update string defined as sr in line 2 to the file name of DNA-PAINT image (2D tiff format)
	- update string defined as ca in line 3 to the file name of averaged Ca2+ image (2D tiff format)
step 4:	(optional) change the integer number of variable shrink to alter the factor by which the image data will be downsampled for display
step 5: Iteratively run the script 
	- each time, update x & y manual shift vectors (in pixels) in variable defined as svec on line 11. Default values optimal for supplied example data. Note that the shift values are sensitive to the shrink factor (defined in step 4)
	- (optional) uncomment the command 'end' on line 31 to shortcut the script for this iterative loop

step 6: close idlde.exe and open correlative_spark.pro with a new idlde.exe session
	- update string defined as dir in line 1 to the current path containing the files. Make sure that this is the same path as the outputs of the align.pro script
	- update string defined as ca in line 12 to the file name of averaged Ca2+ image (2D tiff format)
	- update string defined as dots in line 38 to the .txt file output from the xyspark localisation of the Ca2+ sparks (it is essential that the x-y dimensions of the Ca2+ time series is identical to those of the array ca)
	- update string defined as dotc in line 57 to the .txt file containing the x-y coordinates of the centroids of punctate RyRs detected from the DNA-PAINT image used in the previous script
step 7: Compile and run the script
	- during the spark-by-spark analysis, the display window will show an overlay of the averaged Ca2+ image (red), The RyR coordinates (in blue circles), and the local Ca2+ spark area that is being analysed (green)
	
Outputs - Window 3 will display a scattergram of the spark mass vs the local RyR count. This will have excluded sparks which were filtered out with an exclusion criteria (detailed in Hurley et al manuscript)
	- sparkmass.csv containing a table of RyR counts and fitting parameter for each spark detected in the Ca2+ time series.



# Copyright Izzy Jayasinghe, 2020
# i.jayasinghe@leeds.ac.uk
#
# This code is free software: you can redistribute it and/or modify
# it. It is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  


