# This module defines various functions to facilitate the main functionality.
# Author: Chrys Adams

#import library modules
import functools
import numpy as np
import arcpy

# A list of offsets to calculate each direct neighbor around a cell
Offsets = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]

# A function to calculate a set of offsets for a
# neighborhood of a given radius around a cell.
# Radius defaults to 1 which would be the same as Offsets list above.
def offset_generator(radius=1):
    offsets = []
    for i in range (-radius, radius+1): # x offset
        for j in range (-radius, radius+1): # y_offset
            if i == 0 and j == 0:
                continue
            offsets.append((i, j))
    return offsets

# Converts a latitude, longitude point to a
# pair of x, y indices into the provided raster.
def coord_to_index(point, raster):
    if not raster.extent.contains(point):
        return None
    
    x = int((point.X - raster.extent.XMin) / raster.meanCellWidth)
    y = int((raster.extent.YMax - point.Y) / raster.meanCellHeight)
    return (x, y)


# Converts the x, y indices into a raster into a latitude, longitude point
def index_to_coord(x, y, raster):

    x = raster.extent.XMin + ((float(x)+0.5)*raster.meanCellWidth)
    y = raster.extent.YMin + ((float(y)+0.5)*raster.meanCellHeight)
    return arcpy.Point(x, y)




# Merges two rasters into a single raster using python loops.
# This function is very slow. Do not use it.
def merge_slow(raster1, raster2):
    spatialReference = arcpy.Describe(raster1).spatialReference
    # Calculate extent of resulting raster.
    new_extent = arcpy.Extent(
        XMin=min(raster1.extent.XMin, raster2.extent.XMin),
        XMax=max(raster1.extent.XMax, raster2.extent.XMax),
        YMin=min(raster1.extent.YMin, raster2.extent.YMin),
        YMax=max(raster1.extent.YMax, raster2.extent.YMax),
    )
    # Convert inputs into numpy arrays
    arr1 = arcpy.RasterToNumPyArray(raster1)
    arr2 = arcpy.RasterToNumPyArray(raster2)

    # Calculate the size and dimensions of the resulting raster
    x_increment = raster1.meanCellWidth
    y_increment = raster1.meanCellHeight
    new_dimensions = {
        "height": int(new_extent.height / y_increment),
        "width": int(new_extent.width / x_increment)
    }
    new_arr = np.full((new_dimensions['height'], new_dimensions['width']), None, float)

    # Copy array 1 into the result array
    arr1_offset = {
        "x": int((raster1.extent.XMin - new_extent.XMin)/x_increment),
        "y": int((raster1.extent.YMax - new_extent.YMax)/y_increment)
    }
    for index, value in np.ndenumerate(arr1):
        if value < 0:
            value = None
        new_arr[index[0]-arr1_offset['y'], index[1]+arr1_offset['x']] = value

    # Copy array 2 into the result array
    arr2_offset = {
        "x": int((raster2.extent.XMin - new_extent.XMin)/x_increment),
        "y": int((raster2.extent.YMax - new_extent.YMax)/y_increment)
    }
    for index, value in np.ndenumerate(arr2):
        if value < 0:
            value = None
        new_arr[index[0]-arr2_offset['y'], index[1]+arr2_offset['x']] = value

    # Convert the result array into a raster. Use extent position, resolution and spatial reference calculated earlier
    result_a = arcpy.NumPyArrayToRaster(new_arr, new_extent.lowerLeft, x_increment, y_increment)
    arcpy.DefineProjection_management(result_a, spatialReference)
    return result_a

# This function merges two rasters into a single raster using numpy assignment.
# This function has less functionality than merge_rasters(). Prefer that function.
def merge_raster(raster1, raster2):
    spatialReference = arcpy.Describe(raster1).spatialReference

    # Calculate new extent
    new_extent = arcpy.Extent(
        XMin=min(raster1.extent.XMin, raster2.extent.XMin),
        XMax=max(raster1.extent.XMax, raster2.extent.XMax),
        YMin=min(raster1.extent.YMin, raster2.extent.YMin),
        YMax=max(raster1.extent.YMax, raster2.extent.YMax),
    )

    # Convert rasters to arrays
    arr1 = arcpy.RasterToNumPyArray(raster1, nodata_to_value=np.nan)
    arr2 = arcpy.RasterToNumPyArray(raster2, nodata_to_value=np.nan)

    # Calculate resolution and dimensions of result array
    x_increment = raster1.meanCellWidth
    y_increment = raster1.meanCellHeight
    new_dimensions = {
        "height": int(new_extent.height / y_increment),
        "width": int(new_extent.width / x_increment)
    }
    new_arr = np.full((new_dimensions['height'], new_dimensions['width']), np.nan, float)

    # Copy arr1 into result array
    arr1_offset = {
        "x": int((raster1.extent.XMin - new_extent.XMin)/x_increment),
        "y": -int((raster1.extent.YMax - new_extent.YMax)/y_increment)
    }
    # Calculate index offsets for numpy slice assignment for arr1.
    tgt_x_begin = arr1_offset['x']
    tgt_x_end = arr1_offset['x'] + arr1.shape[1]
    tgt_y_begin = arr1_offset['y']
    tgt_y_end = arr1_offset['y'] + arr1.shape[0]
    src_x_begin = 0
    src_x_end = arr1.shape[1]
    src_y_begin = 0
    src_y_end = arr1.shape[0]
    # Assign the source array to the result array at the calculated location
    new_arr[
        tgt_y_begin: tgt_y_end,
        tgt_x_begin: tgt_x_end
    ] = arr1[
        src_y_begin: src_y_end,
        src_x_begin: src_x_end
    ]

    # Copy arr2 into result array
    arr2_offset = {
        "x": int((raster2.extent.XMin - new_extent.XMin)/x_increment),
        "y": -int((raster2.extent.YMax - new_extent.YMax)/y_increment)
    }
    # Calculate index offsets for numpy slice assignment for arr2.
    tgt_x_begin = arr2_offset['x']
    tgt_x_end = arr2_offset['x'] + arr2.shape[1]
    tgt_y_begin = arr2_offset['y']
    tgt_y_end = arr2_offset['y'] + arr2.shape[0]
    src_x_begin = 0
    src_x_end = arr2.shape[1]
    src_y_begin = 0
    src_y_end = arr2.shape[0]
    # Assign the source array to the result array at the calculated location
    new_arr[
        tgt_y_begin: tgt_y_end,
        tgt_x_begin: tgt_x_end
    ] = arr2[
        src_y_begin: src_y_end,
        src_x_begin: src_x_end
    ]

    # Convert the result array into a raster. Use extent position, resolution and spatial reference calculated earlier
    result_a = arcpy.NumPyArrayToRaster(new_arr, new_extent.lowerLeft, x_increment, y_increment)
    arcpy.DefineProjection_management(result_a, spatialReference)
    return result_a

# Helper function for merge_rasters()
# Takes two extents or rasters, and returns a larger extent that encloses both.
# The second extent will be projected into the spatial reference of the first extent.
# The resulting extent will have the same spatial reference as the first parameter.
def combine_extents(extent1: arcpy.Extent, extent2: arcpy.Extent):
    # If either parameter is a raster, get that raster's extent
    if isinstance(extent1, arcpy.Raster):
        extent1 = extent1.extent
    if isinstance(extent2, arcpy.Raster):
        extent2 = extent2.extent
    # Project the second parameter into the spatial reference of the first
    extent2 = extent2.projectAs(extent1.spatialReference)
    # Return the extent that encloses both
    return arcpy.Extent(
        XMin=min(extent1.XMin, extent2.XMin),
        XMax=max(extent1.XMax, extent2.XMax),
        YMin=min(extent1.YMin, extent2.YMin),
        YMax=max(extent1.YMax, extent2.YMax),
        spatial_reference=extent1.spatialReference,
    )


# Merges all rasters in a list of rasters into a single raster.
# The resulting raster will have the spatial reference as the first raster in the list.
def merge_rasters(rasters):
    if len(rasters) == 0:
        return None
    
    # assuming each raster is the same resolution,
    # use these parameters from the first.
    x_increment = rasters[0].meanCellWidth
    y_increment = rasters[0].meanCellHeight
    
    # Use the helper funtion to combine all extents in the list into a single enclosing extent.
    new_extent = functools.reduce(combine_extents, rasters)
    spatialReference = new_extent.spatialReference

    # Calculate array size from the extent
    new_dimensions = {
        "height": int(new_extent.height / y_increment),
        "width": int(new_extent.width / x_increment)
    }
    new_arr = np.full((new_dimensions['height'], new_dimensions['width']), np.nan, float)
    
    # Copy each raster into the new array
    for raster in rasters:
        # Convert raster to array
        arr = arcpy.RasterToNumPyArray(raster, nodata_to_value=np.nan)
        # Calculate indices
        arr_offset = {
            "x": int((raster.extent.XMin - new_extent.XMin)/x_increment),
            "y": -int((raster.extent.YMax - new_extent.YMax)/y_increment)
        }
        
        # Calculate numpy slice assignment
        tgt_x_begin = arr_offset['x']
        tgt_x_end = arr_offset['x'] + arr.shape[1]
        tgt_y_begin = arr_offset['y']
        tgt_y_end = arr_offset['y'] + arr.shape[0]
        src_x_begin = 0
        src_x_end = arr.shape[1]
        src_y_begin = 0
        src_y_end = arr.shape[0]
        # Assign to result array
        new_arr[
            tgt_y_begin: tgt_y_end,
            tgt_x_begin: tgt_x_end
        ] = arr[
            src_y_begin: src_y_end,
            src_x_begin: src_x_end
        ]

    # Convert array into raster with calculated resolution and spatial reference
    result = arcpy.NumPyArrayToRaster(new_arr, new_extent.lowerLeft, x_increment, y_increment)
    arcpy.DefineProjection_management(result, spatialReference)
    return result
