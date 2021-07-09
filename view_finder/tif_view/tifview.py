##############
#Title: tifview
#Author: Liam Monninger SID 004748349
#Purpose: The tifview module
##############
import math
import requests
import numpy as np
from osgeo import gdal
from osgeo import osr
from osgeo import ogr

class Bounds:

    '''
    Initializes a set of bounds

    @ param south is the southern limit for the bounded area
    @ param west is the western limit for the bounded area
    @ param north is the northern limit for the bounded area
    @ param east is the eastern limit for the bounded area
    '''
    def __init__(self, south, west, north, east):
        self.north = north
        self.west = west
        self.south = south
        self.east = east

class DEM:
    # projection and pixel size are constant
    EPSG = 4326
    EPSG_STRING = '4326'

    '''
    Creates a DEM using requests and saves said DEM to a given location

    @ param bounds are the lng, lat bounds for the DEM
    @ param filepath is the output filepath
    @ res is an optional parameter for the pixel width and height in DD

    '''
    def __init__(self, bounds, filepath, *, res=0.0008333333333333334):

        # store those variables for future access 
        self.bounds = bounds
        self.filepath = filepath

        # interpret the resolution variable as follows
        self.pxlXsize = res
        self.pxlYsize = res # the pixels will be squares; I use this notation for ease of understanding

        # get the geotransform tuple
        self.GeoTransform = (self.bounds.west, self.pxlXsize, 0, self.bounds.north, 0, -self.pxlYsize)

        # we'll start by making two parallel matrices each containing the lat on lng values which satisfy the bounds
        # that is...
        lngs = np.arange(self.bounds.west, self.bounds.east, self.pxlXsize)
        lats = np.arange(self.bounds.south, self.bounds.north, self.pxlYsize)
        lats = np.flip(lats)
        self.rasterXSize = len(lngs) # no need for abs if done this way
        self.rasterYSize = len(lats) # no need for abs if done this way
        
        # target projection
        target = osr.SpatialReference()
        target.ImportFromEPSG(self.EPSG)

        # set driver
        driver = gdal.GetDriverByName('GTiff')
        filepath = self.filepath + ".elevation.tif"
        output = driver.Create(filepath, self.rasterXSize, self.rasterYSize, 1, gdal.GDT_Float64)
        
        # we'll then make a mesh of th lngs and lats matrices
        LNGS, LATS = np.meshgrid(lngs, lats)

        # and initialize a matrix of 0's with the same dimension
        self.elev_matrix = np.zeros_like(LNGS) # Note: this could just as well be LATS
        
        with np.nditer(self.elev_matrix, flags=['multi_index'], op_flags=['writeonly']) as pixel: # with a numpy iterator that allows for indexing and writing
            while not pixel.finished: # for all elements in the elev_matrix
                lng = LNGS[pixel.multi_index[0], pixel.multi_index[1]] 
                lat = LATS[pixel.multi_index[0], pixel.multi_index[1]]

                # make requests for elevation values
                r = requests.get('https://apps.gis.ucla.edu/elevation/api/v1/lookup?locations=%s,%s'%(lat, lng))
                result = r.json()
                elevation = result['results'][0]['elevation'] # get elevation
                elevation*=3.28084 # convert to feet
                pixel[0]=elevation
                pixel.iternext()

        # write to tiff
        output.SetGeoTransform(self.GeoTransform)
        output.SetProjection(target.ExportToWkt())
        output.GetRasterBand(1).WriteArray(self.elev_matrix)
        output.GetRasterBand(1).SetNoDataValue(-9999)
        output = None


'''
This is used to store lng, lat coordinates and their 'neighbors'

Note: I do 
'''
class Point:

    POINT_STEP = 0.0008333333333333334  * 3 # this is roughly 90 meters as .0008333333333333334 DD is about 30 meters

    '''
    Initializes a point using provided lng and lat
    '''
    def __init__(self, lng, lat):

        if lng < -180.0: # if initial lng value too far west
            lng = -180.0 # set point to west limit
        elif lng > 180.0: # if initial lng value too far east
            lng = 180.0 # set point to east limit

        if lat < -90.0: # if initial lng value too far south
            lat = -90.0 # set point to south limit
        elif lat > 90.0: # if initial lng value too far north
            lat = 90.0 # set point to north limit

        self.lng = lng
        self.lat = lat

        # We will now add the 'neighboring' points as tuples
        self.neighbors = set()

        left = lng - self.POINT_STEP
        right = lng + self.POINT_STEP

        # take care of out of bounds issues
        if left < -180.0: # if you've gone too far west
            left = 180.0 - (POINT_STEP - (lng - (-180.0))) # then your coordinate can be specified as the residual of your step less than 180 DD East
        if right > 180.0: # if you've gone too far east 
            right = -180.0 + (POINT_STEP - (180.0 - lng)) # then your then your coordinate can be specified as the residual of your step added to -180.0 DD West

        top = lat + self.POINT_STEP
        bottom = lat - self.POINT_STEP

        # take care of out of bounds issues
        if bottom < -90.0: # if you've gone too far south, well this is quite and issue; basically, we just want to set everything to -90.0 from here on out
            bottom = -90.0
        if top > 90.0: # if you've gone too far north
            top = 90.0 # same as above

        # Neighbors will be orthogonal and diagonal neighbors 
        # add everything but the Point itself to the neighbors list

        # top row of neighbors
        self.neighbors.add((left, top))
        self.neighbors.add((self.lng, top))
        self.neighbors.add((right, top))

        # middle row
        self.neighbors.add((left, lat))
        self.neighbors.add((right, lat))

        # bottom row
        self.neighbors.add((left, bottom))
        self.neighbors.add((self.lng, bottom))
        self.neighbors.add((right, bottom))

    '''
    Print function for Point
    '''
    def __str__(self):
        return "( {} , {} )".format(self.lng, self.lat)

    '''
    Reproduce function for Point
    '''
    def __repr__(self):
         return "( {} , {} )".format(self.lng, self.lat)

    '''
    Hash function for Point
    '''
    def __hash__(self):
        return hash(repr(self))
    '''
    Equivalence function for Point
    '''
    def __eq__(self, other):
        return (abs(self.lng - other.lng) < 0.0008333333333333334) and (abs(self.lat - other.lat) < 0.0008333333333333334)  # I've defined this as such in order to deal with round off errors
    '''
    Less than function for point
    '''
    def __lt__(self, other):
        return math.sqrt((self.lng**2) + (self.lat**2)) < math.sqrt((other.lng**2) + (other.lat**2))
    
    '''
    Sets the elevation value of the point

    @ parm elevation is the elevation value you would like to provide
    '''
    def setElevation(self, elevation):
        self.elevation = elevation

    '''
    Casts the neighbor tuples as Points themselves
    
    @ returns the neighbors as points
    '''
    def neighborsAsPoints(self):
        neighbors_as_points = set()
        for point in self.neighbors: # for each tuple in neighbors
            neighbors_as_points.add(Point(point[0], point[1]))
        return neighbors_as_points

'''
The rateView() algorithm core. This function is called recursively to compute the view score for a given point

Note: One of the big things I worked on was providing consistent scores, which was problematic at first on account of round-off errors
and set element access randomness. At present, I seem to have a solution that is working, though perhaps not altogether very rapid.

@ param point is the lng, lat point for which the algorithm is computing the score
@ param step is the step in recursion of the algorithm
@ param max_steps is the the maximum number of recursive steps to be taken
Note: I could have implemented this using a while-loop, and perhaps should have done; that said, the algorithm runs fast enough for present purposes, and I believe
is a little more intuitive and a little less cumbersome when implemented recursively
@ param k is an optional parameter which helps to reduce the score

@ returns the view score for your point
'''
def _rateView(point, step, max_steps, point_set, *, k=5000):

    if(step > max_steps): # if we've gone far enough
        return 0 # we've reached then edge of our search, so don't add to the sum

    neighbors = point.neighborsAsPoints() - point_set # get all neighbors save those already visited
    point_set.update(neighbors) # we need to add these neighbors to the point_set before looping through each one; otherwise, we may visit points which should already have been marked visited in future steps of recursion

    neighbors = list(neighbors) # I've done this to try and reduce the variation in the score 

    neighbors.sort() # this is also done in aims of providing consistency

    elev_dif_sum = 0

    for neighbor in neighbors: # for each neighbor

        # make requests for elevation values
        r = requests.get('https://apps.gis.ucla.edu/elevation/api/v1/lookup?locations=%s,%s'%(neighbor.lat, neighbor.lng))
        result = r.json()
        elevation = result['results'][0]['elevation'] # get elevation
        elevation*=3.28084 # convert to feet
        neighbor.setElevation(elevation)
        elev_dif = point.elevation - neighbor.elevation

        # this is the tricky part: we would like to normalize the score; however, if we do so based on elevation alone, we will give a great advantage to
        # points at lower elevation
        # for the moment, I've chosen to not really normalize this score, so much as reduce it
        # I've selected a k of 5000 feet as an arbitrarily large elevation difference which the user can overwrite if they please
        neighbor_elev_dif_sum = elev_dif/(k*(step + 1)) # the elev_dif_sum to be added is normalized by the elevation of the point and weighted by the step in recursion

        elev_dif_sum += neighbor_elev_dif_sum + _rateView(neighbor, step + 1, max_steps, point_set, k=k) # get all of the values from the Points on this branch of recursion

    return elev_dif_sum

'''
The rateView() wrapper function.

@ param lng is the longitude value for the point being rated
@ param lat is the latitiude value for the point being rated
@ param max_steps is the maximum number of recursive steps the user would like to have taken in computing the score
Note: the max_steps also essentially accounts for how far the user would like to consider their view; e.g., 6 max steps is a maximum of 6 90 meter steps
away from the start point or 540 meters

@ returns the view score for the point
'''
def rateView(lng, lat, max_steps):

    start_point = Point(lng, lat) # the start point has a value of zero
    r = requests.get('https://apps.gis.ucla.edu/elevation/api/v1/lookup?locations=%s,%s'%(start_point.lat, start_point.lng))
    result = r.json()
    elevation = result['results'][0]['elevation'] # get elevation
    elevation*=3.28084 # convert to feet
    start_point.setElevation(elevation)

    point_set = set()

    point_set.add(start_point)

    return _rateView(start_point, 0, max_steps, point_set) # run_recursion
