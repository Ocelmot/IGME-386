# Main script for Watershed analysis project
# Author: Chrys Adams

# This module imports the other two modules,
# instantiates and configures the Watershed object
# and runs the analysis


#import library modules
import sys
import os
import time
import importlib

import arcpy

# setup script metadata for project imports for use in arcpy
project = arcpy.mp.ArcGISProject("CURRENT")
project_directory = os.path.dirname(project.filePath)
scripts_directory = os.path.join(project_directory, "scripts")
data_directory = os.path.join(project_directory, "data")
# add the scripts directory to the search path for python modules
# to make importing them easier. This assumes that the scripts are
# located in a directory called scripts in the project directory.
if sys.path[0] != scripts_directory: # fresh run
    sys.path.insert(1, scripts_directory)

    # Import project modules
    import utils
    import watershed
# This was an attempt to make the arcgis python console reload
# the modules. It does not work. The jupyter notebooks may be
# better at this.
else: # rerun requires reimports
    utils = importlib.reload(utils)
    watershed = importlib.reload(watershed)


# To allow overwriting the outputs change the overwrite option to true.
arcpy.env.overwriteOutput = True

# This function was used to develop the merge functionality
def test_operations():
    # Test operations!
    raster1 = arcpy.Raster('USGS_13_n44w078.tif')
    raster2 = arcpy.Raster('USGS_13_n43w078.tif')
    raster3 = arcpy.Raster('USGS_13_n43w077.tif')

    # slow_start = time.time()
    # result_a = utils.merge_slow(raster1, raster2)
    # slow_end = time.time()

    print("using merge")
    fast_start = time.time()
    result_b = utils.merge_raster(raster1, raster2)
    # result_b1 = utils.merge_raster(result_b, raster3)
    fast_end = time.time()

    print("using multi merge")
    multi_start = time.time()
    result_c = utils.merge_rasters([raster1, raster2])
    multi_end = time.time()

    # print("Slow: ", slow_end-slow_start)
    print("Fast: ", fast_end-fast_start)
    print("Multi: ", multi_end-multi_start)

# This function manages and runs the main functionality of this script
def find_watershed():
    # Get raster objects from the current arcgis project
    # High res
    raster1 = arcpy.Raster('USGS_13_n44w078.tif')
    raster2 = arcpy.Raster('USGS_13_n43w078.tif')

    # Low res
    raster1_LR = arcpy.Raster('USGS_1_n44w078.tif')

    # Instantiate the watershed object and include the desired raster datasets
    ws = watershed.Watershed()
    ws.IncludeRaster(raster1)
    # ws.IncludeRaster(raster2)
    # ws.IncludeRaster(raster1_LR)

    # Create the starting point for the watershed analysis.
    # This point is used to determine which watershed to analyse.
    point = arcpy.Point(-77.6989588, 43.0640209)

    # Run and return the results of the analysis with the selected parameters.
    # return ws.find_watershed_simple(point, slope_tolerance=0.05)
    return ws.find_watershed_simple(point, slope_tolerance=0.00, neighborhood_radius=3)
    # return ws.find_watershed_rainfill(point)


# If this is the main file, execute the desired functionality
if __name__ == "__main__":
    # test_operations()
    (ws, inflows) = find_watershed()
