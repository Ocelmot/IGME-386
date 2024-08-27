# This module defines the bulk of how the watershed is calculated
# Author: Chrys Adams

from collections import deque

import numpy as np
import arcpy

import utils

# The Watershed class stores the rasters and starting point for any watershed analyses.
class Watershed:
    # Create a new Watershed class
    def __init__(self) -> None:
        self.rasters = []
        self.complete_raster = None
        self.point = None

    # Add a raster to the set of rasters to use for analysis
    def IncludeRaster(self, raster):
        self.rasters.append(raster)
        self.complete_raster = None
        # If there is only one raster, it must be the complete raster.
        if len(self.rasters) == 1:
            self.complete_raster = self.rasters[0]

    # Combine all available rasters into a single complete raster for analysis.
    def compile_rasters(self):
        if self.complete_raster is None:
            print("Watershed: merging rasters...")
            self.complete_raster = utils.merge_rasters(self.rasters)
            print("Watershed: merging rasters: Done")

    # Set the point to be used in the analysis
    def set_point(self, point):
        self.point = point

    # Find the minimum point in all rasters.
    def find_min(self):
        self.compile_rasters()
        if self.complete_raster is None:
            return None
        print("Watershed: finding minimum index...")
        # Convert complete raster into array
        arr = arcpy.RasterToNumPyArray(self.complete_raster, nodata_to_value=np.nan)
        min_val = float("inf")
        min_index = None
        # iterate through array searching for the minimum
        for index, value in np.ndenumerate(arr):
            if value is np.nan:
                continue
            if value < min_val:
                min_val = value
                min_index = (index[1], index[0])
        print("Watershed: found minimum index:", min_index)
        # Return the index of the minimum value
        return min_index
    
    # Find the watershed of a point using the simple algorithm to find all uphill cells
    # The point argument can be used to override the point set with set_point()
    # The slope_tolerance defines how much lower a cell can be from the current cell and still be considered.
    # The neighborhood radius defines how many cells around the current cell will be considered.
    # Returns: (Watershed, inflow) tuple
    def find_watershed_simple(self, point=None, slope_tolerance = 0.0, neighborhood_radius=1):
        # Prepare the raster variable
        self.compile_rasters()
        if self.complete_raster is None:
            return None
        
        # Prepare the point variable
        if point is None:
            point = self.point
        if point is None:
            point = self.find_min()
        if point is None:
            return None

        # Convert point to indices into the complete raster
        if isinstance(point, arcpy.Point):
            point = utils.coord_to_index(point, self.complete_raster)

        # Create arrays to operate on
        data = arcpy.RasterToNumPyArray(self.complete_raster)
        watershed = np.full_like(data, None)
        inflows = np.full_like(data, None)

        # Define a queue of candidate points to process
        process_queue = deque()
        process_queue.append(point)

        # Calculate offsets for the neighborhood using the radius
        neighborhood_offsets = utils.offset_generator(neighborhood_radius)

        # Process candidate points
        while True:
            try:
                current_point = process_queue.popleft()
            except:
                break

            # if this point is already included, skip as it was already processed
            if watershed[current_point[1], current_point[0]] == 1:
                continue
            # include this point, then determine which of its neighbors should also be considered
            watershed[current_point[1], current_point[0]] = 1
            if current_point[0] == 0 or current_point[0] >= data.shape[1] - 1:
                inflows[current_point[1], current_point[0]] = 1
            if current_point[1] == 0 or current_point[1] >= data.shape[0] - 1:
                inflows[current_point[1], current_point[0]] = 1

            current_value = data[current_point[1], current_point[0]]
            for offset in neighborhood_offsets:
                offset_x = current_point[0]+offset[0]
                if offset_x < 0 or offset_x >= data.shape[1]:
                    continue
                offset_y = current_point[1]+offset[1]
                if offset_y < 0 or offset_y >= data.shape[0]:
                    continue
                if watershed[offset_y, offset_x] == 1:
                    continue
                offset_value = data[offset_y, offset_x]
                if offset_value >= (current_value - slope_tolerance):
                    process_queue.append((offset_x, offset_y))

        # Calculate and assign extent, resolution, and spatial reference to output arrays.
        lower_left = self.complete_raster.extent.lowerLeft
        x_increment = self.complete_raster.meanCellWidth
        y_increment = self.complete_raster.meanCellHeight
        watershed = arcpy.NumPyArrayToRaster(watershed, lower_left, x_increment, y_increment)
        arcpy.DefineProjection_management(watershed, self.complete_raster.spatialReference)
        inflows = arcpy.NumPyArrayToRaster(inflows, lower_left, x_increment, y_increment)
        arcpy.DefineProjection_management(inflows, self.complete_raster.spatialReference)
        return (watershed, inflows)

    # Find watershed by simulating trickles of water flowing down the slope until they meet the starting point.
    # The point parameter can override the starting point defined by the class.
    # The slope_tolerance governs how fast concavities are filled in.
    # Returns a raster defining the watershed.
    #
    # The overall strategy can be described as follows:
    # have point flow downhill until it pools (all cells within radius are within tolerance)
    # have each cell flow downhill until it pools, if its path intersects the point trickle, include it
    # if the cell successfully pools without intersecting the point trickle, exclude it
    # if the cell reaches a concavity, but is not pooled, then elevate this cell to the height of its lowest neighbor
    def find_watershed_rainfill(self, point=None, slope_tolerance = 0.0):
        print("Watershed: finding watershed (rainfill method) ...")
        # Prepare the raster variable
        self.compile_rasters()
        if self.complete_raster is None:
            return None
        
        # Prepare the point variable
        if point is None:
            point = self.point
        if point is None:
            point = self.find_min()
        if point is None:
            return None

        # Convert point to indices
        if isinstance(point, arcpy.Point):
            point = utils.coord_to_index(point, self.complete_raster)

        # Create arrays for processing
        data = arcpy.RasterToNumPyArray(self.complete_raster)
        watershed = np.full_like(data, None)

        print("Watershed: generating start trickle ...")
        # trickle start point until it pools
        (path, pool) = trickle_fill(data, point)
        set_indices(watershed, path)
        set_indices(watershed, pool)
        print("Watershed: generated start trickle")

        # This algorithm as it stands is unlikely to get even this far in a short amount of time.
        # Uncomment the following to only calculate the starting trickle.
        # lower_left = self.complete_raster.extent.lowerLeft
        # x_increment = self.complete_raster.meanCellWidth
        # y_increment = self.complete_raster.meanCellHeight
        # ret = arcpy.NumPyArrayToRaster(watershed, lower_left, x_increment, y_increment)
        # arcpy.DefineProjection_management(ret, self.complete_raster.spatialReference)
        # print("Watershed: found watershed (rainfill method)")
        # return ret

        print("Watershed: filling candidate queue ...")
        candidate_queue = deque()
        for index, _ in np.ndenumerate(data): # all points are candidates at first
            candidate_queue.append(index)
        print("Watershed: filled candidate queue")

        print("Watershed: trickling candidate cells ...")
        print("Watershed: initial candidates:", len(candidate_queue))
        candidates_trickled = 0
        while True:
            try:
                current_point = candidate_queue.popleft()
            except:
                break

            candidates_trickled += 1
            if candidates_trickled % 10000 == 0:
                print("trickled", candidates_trickled, "candidates...")
                print(len(candidate_queue), "candidates remain.")
            
            # for each candidate index, find its path and pool
            (path, pool) = trickle(data, current_point)
            # if it pools, check if the final step in the path has intersected the existing watershed
            if pool is not None:
                intersects_watershed = False
                for cell in path:
                    if watershed[cell[1], cell[0]] == 1:
                        intersects_watershed = True
                        break
                # if it has, add the path to the watershed
                if intersects_watershed:
                    set_indices(watershed, path)
                # since this cell has pooled, it is finished and not re-added to the candidate set
            # if it failed to pool, increase the height of where it ended
            else:
                last_indices = path[-1]
                current_height = data[last_indices[1], last_indices[0]]
                lowest_neighbor_height = float('inf')
                for offset in utils.Offsets:
                    offset_x = last_indices[0]+offset[0]
                    if offset_x < 0 or offset_x >= data.shape[1]:
                        continue
                    offset_y = last_indices[1]+offset[1]
                    if offset_y < 0 or offset_y >= data.shape[0]:
                        continue
                    neighbor_height = data[offset_y, offset_x]
                    
                    if neighbor_height < lowest_neighbor_height:
                        lowest_neighbor_height = neighbor_height
                # fill to lowest neighbor
                if current_height < lowest_neighbor_height:
                    data[last_indices[1], last_indices[0]] = lowest_neighbor_height
                else:
                    data[last_indices[1], last_indices[0]] += slope_tolerance / 2
        print("Watershed: all candidates pooled")
        

        # Prepare the return raster
        lower_left = self.complete_raster.extent.lowerLeft
        x_increment = self.complete_raster.meanCellWidth
        y_increment = self.complete_raster.meanCellHeight
        ret = arcpy.NumPyArrayToRaster(watershed, lower_left, x_increment, y_increment)
        arcpy.DefineProjection_management(ret, self.complete_raster.spatialReference)
        print("Watershed: found watershed (rainfill method)")
        return ret


# Watershed specific helper functions
    
# Calculate the downhill flow from a single cell.
# start_indices is the starting point (x,y)
# pool_radius is the required size of the pool to determine the stopping point
# min_flow is how level the stopping pool needs to be
# Returns: (path, pool) sets of coordinates to represent the path and resulting pool.
def trickle(arr, start_indices, pool_radius=5, min_flow=0.05):
    # print("starting trickle at", start_indices)
    path = []
    pool = None

    current_indices = start_indices
    while True:
        path.append(current_indices)
        # try to flow to lowest neighbor
        current_height = arr[current_indices[1], current_indices[0]]
        flow_cell = None
        flow_magnitude = min_flow
        for offset in utils.Offsets:
            offset_x = current_indices[0]+offset[0]
            if offset_x < 0 or offset_x >= arr.shape[1]:
                continue
            offset_y = current_indices[1]+offset[1]
            if offset_y < 0 or offset_y >= arr.shape[0]:
                continue
            offset_height = arr[offset_y, offset_x]
            flow = current_height - offset_height
            if flow > flow_magnitude:
                flow_magnitude = flow
                flow_cell = (offset_x, offset_y)
        # if flowed, continue and try to flow some more
        if flow_cell is not None:
            current_indices = flow_cell
            continue   
        # failed to flow, check if pooled
        hit_edge = False
        pool_flow = False
        pool_coords = []
        # For each cell within the pool radius, check if all cells are within the min_flow range
        for i in range (-pool_radius, pool_radius+1): # x offset
            for j in range (-pool_radius, pool_radius+1): # y_offset
                # remember to not count 0,0
                offset_x = current_indices[0] + i
                if offset_x < 0 or offset_x >= arr.shape[1]:
                    hit_edge = True
                    continue
                offset_y = current_indices[1] + j
                if offset_y < 0 or offset_y >= arr.shape[0]:
                    hit_edge = True
                    continue
                pool_coords.append((offset_x, offset_y))
                offset_height = arr[offset_y, offset_x]
                flow = current_height - offset_height
                if abs(flow) > min_flow:
                    pool_flow = True
        # Hitting the edge of the data does not require a full pool
        if hit_edge or not pool_flow:
            pool = pool_coords
        # Even if the trickle did not pool, we return since it cannot trickle further
        break
    return (path, pool)

# Trickle as in the trickle() function, but if the function failed to pool,
# attempt to fill in the concavity until pooling has occured or flow has been restored.
# start_indices is the starting point (x,y)
# pool_radius is the required size of the pool to determine the stopping point
# min_flow is how level the stopping pool needs to be
# Returns: (path, pool) sets of coordinates to represent the path and resulting pool.
def trickle_fill(arr, start_indices, pool_radius=5, min_flow=0.05):
    while True:
        (path, pool) = trickle(arr, start_indices, pool_radius, min_flow)
        if pool is not None:
            return (path, pool)
        # could not pool, must adjust height of last cell
        last_indices = path[-1]
        current_height = arr[last_indices[1], last_indices[0]]
        # Find height of nearest neighbor.
        lowest_neighbor = None
        lowest_neighbor_height = float('inf')
        for offset in utils.Offsets:
            offset_x = last_indices[0]+offset[0]
            if offset_x < 0 or offset_x >= arr.shape[1]:
                continue
            offset_y = last_indices[1]+offset[1]
            if offset_y < 0 or offset_y >= arr.shape[0]:
                continue
            neighbor_height = arr[offset_y, offset_x]
            
            if neighbor_height < lowest_neighbor_height:
                lowest_neighbor = (offset_x, offset_y)
                lowest_neighbor_height = neighbor_height
        # fill to lowest neighbor, or by some minimum if the neighbor was not lower.
        # this is to avoid flat, non-pool areas from getting stuck as the height of
        # the neighbor will be the same as the current cell.
        if current_height < lowest_neighbor_height:
            arr[last_indices[1], last_indices[0]] = lowest_neighbor_height + min_flow
        else:
            arr[last_indices[1], last_indices[0]] += min_flow

# Sets cells in a numpy array to a specified value.
# The indices of the cells to by set are specified by the arr parameter.
# Useful with the trickle functions which return arrays of indices.
def set_indices(arr, indices, value=1):
    for index in indices:
        arr[index[1], index[0]] = value
