###############################################################################################################
## Runfinder is built for use within a web environment
## Calls to runfinder are intended only to be made via its core algortihm findRun()
###############################################################################################################
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
from osgeo import ogr
import datetime
import math
import numpy as np
import ee
import requests
import zipfile
import tempfile
import config
import os
import shutil

###############################################################################################################
# CONSTANTS
RUN_FINDER_EPSG = 3857


###############################################################################################################
#CLASS: DEMMeta
# a lighweight class containing dem raster metadata
class DEMMeta:
    '''
    Initializes a DEMMeta object using a filepath

    @param filepath is the filepath of the dem
    '''
    def __init__(self, output_dir):
        self.output_dir = output_dir
        raster = gdal.Open(filepath)
        self.x_max = raster.RasterXSize - 1
        self.y_max = raster.RasterYSize - 1

###############################################################################################################
# CLASS: GoogleDEM
# fetches a DEM using Earth Engine
class GoogleDEM:
    '''
    initializes GoogleDEM using set of bounds
    @param bounds are the bounds of the DEM request
    '''
    def __init__(self, bounds, output_dir):

        credentials = ee.ServiceAccountCredentials(config.getServiceAccount(), config.getKeyPath())

        ee.Initialize(credentials) # initialize the GEE session

        area = ee.Geometry.Rectangle(bounds)

        img = ee.Image('USGS/SRTMGL1_003').clip(area)

        url = img.getDownloadURL({'name': 'dem', 'crs': 'EPSG:3857' , 'scale': 30}) # get the url of the GEE request with EPSG 3857 as projection

        self.output_dir = output_dir

        filename = output_dir + '/DEM.zip'

        # Download the zip of requested using the GEE url
        r = requests.get(url, stream=True)
        with open(filename, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=1024):
                fd.write(chunk)

        # Extract the GeoTIFF for the zipped download
        z = zipfile.ZipFile(filename)
        z.extractall(output_dir)

###############################################################################################################
# CLASS: DEM
# a multiple utility dem raster class
class DEM:
    '''
    initializes DEM with a meta object, a gdalraster and a gdal_array

    @param filepath is the filepath for the DEM
    '''
    def __init__(self, output_dir):
        if not filepath:
            return
        self.dem_meta = DEMMeta(output_dir)
        self.gdalraster = gdal.Open(output_dir + '/dem.elevation.tif')
        self.gdalraster.GetRasterBand(1).SetNoDataValue(-32767)
        self.gdalarray = gdal_array.LoadFile(output_dir + '/dem.elevation.tif')
        self.x_offset, self.pixel_width, self.rotation_1, self.y_offset, self.rotation_2, self.pixel_height = self.gdalraster.GetGeoTransform() # this is one solution to the projection reprojection problem

    '''
    Gets pixel x and y given a coordinate pair

    @param lng is the longitude value
    @param lat is the latitude value

    @returns the pixel x and y
    '''
    def pixelFromLngLat(self, lng, lat, *, crs=4326):

        # first we need to convert from 4326
        source = osr.SpatialReference()
        source.ImportFromEPSG(crs)

        target = osr.SpatialReference()
        target.ImportFromEPSG(RUN_FINDER_EPSG)

        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(lng, lat)

        transformation = osr.CoordinateTransformation(source, target)

        # transform point
        point.Transform(transformation)

        lng = point.GetX()
        lat = point.GetY()

        # compute approximate pixel 
        x = int((lng - self.x_offset)/self.pixel_width)
        y = int((self.y_offset - lat)/(self.pixel_height * -1))

        return x, y

    '''
    Gives the lng, lat pair for a given pixel

    @pixel  is the pixel being looked at
    
    @returns the lng, lat pair (DD)
    '''
    def lngLatFromPixel(self, pixel):


        lng = pixel.x * self.pixel_width + self.x_offset + (self.pixel_width/2)
        lat = pixel.y * self.pixel_height + self.y_offset + (self.pixel_height/2)

        return lng, lat

    '''
    Gets value from DEM given pixel

    @param pixel is the pixel for which a value is requested
    @returns the value at said pixel
    '''
    def getValue(self, pixel):
        return self.gdalarray[pixel.y, pixel.x]

    '''
    reprojects DEM and slope rasters (if slope rasters are present)

    @returns whether slope rasters were reprojected
    '''
    def quickReproject(self):
        # Transform to WGS 84 / Pseudo-Mercator
        output_raster = self.dem_meta.output_dir + "/projected_dem.elevation.tif"
        gdal.Warp(output_raster, self.dem_meta.output_dir + '/dem.elevation.tif', dstSRS='EPSG:3857') 
        self.gdalraster = gdal.Open(output_raster) # set DEM raster to the projected raster 
        self.gdalarray = gdal_array.LoadFile(output_raster) # set DEM array to the projected array
        
        # if slope has already been added, go ahead and update it
        try:
            if self.slope_filepath != None: # if the slope has already been set
                self.addSlope()
                return True
        except:
            return False

    '''
    adds slope raster to DEM object
    '''
    def addSlope(self):
        self.slope_filepath = self.dem_meta.output_dir + '/slope.elevation.tif'
        gdal.DEMProcessing(self.slope_filepath, self.dem_meta.output_dir + '/dem.elevation.tif', 'slope')
        self.gdalraster_slope = gdal.Open(self.slope_filepath)
        self.gdalarray_slope = gdal_array.LoadFile(self.slope_filepath)
        return self.gdalraster_slope, self.gdalarray_slope

    '''
    the core recursive algorithm for Run Finder; finds a swath of pixels which represents a good run

    @param pixels is the set of good pixels; initially this will only be the start pixel; this parameter changes over the recursion
    @run_set is the set of good pixels, representing the swath; this parameter changes over the recursion
    @param k is the the constant by which the absolute maximum elevation will be determined; this parameter changes over the recursion
    @param R is the constant by which the relative maximum elevation will be determined; this parameter remains constant over the recursion
    @param max_elev is the absolute maximum elevation which will be applied at the relevant step of recursion; this parameter changes over the recursion
    @param rel_max_elev is the relative maximum elevation which will be applied at the relevant step of recursion; this parameter changes over the recursion
    @param MAX_ELEV_DROP_REDUCER is a constant which can be used to limit the drop of the; this parameter remains constant over the recursion
    @param MAX_SLOPE is the maximum slope allowed for any pixel; this parameter remains constant over the recursion
    @param SLOPE_LIKE is the maximum difference between consecutive slope values allowed; this parameter remains constant over the recursion

    @returns run_set
    '''
    def findRunSwath(self, pixels, run_set, k, R, max_elev, rel_max_elev, MAX_ELEV_DROP_REDUCER, MAX_SLOPE, SLOPE_LIKE):

        try: # check to see that the raster has the needed properties
            self.dem_meta
            self.gdalraster
            self.gdalarray
            self.gdalraster_slope
            self.gdalarray_slope
        except: # ends the function and returns false
            return False

        if len(pixels) < 1: # if there aren't pixels in the run_set break recursion
            return run_set

        new_value = list(pixels)[0].value + 1 # determine the value that will be used to set the raster

        # initialize the sets that will be added to
        good_neighbors = set()
        good_neighbors_elev = set() # this is for determining rel_max_elev
        
        for pixel in pixels: # for each pixel in the set of pixels
            neighbors = pixel.neighborsAsPixels(new_value) # get all the pixel's neighbors
            value_p = self.gdalarray[pixel.y, pixel.x] # get the elevation value of the pixel
            value_p_s = self.gdalarray_slope[pixel.y, pixel.x] # get the slope value of the pixel
            for neighbor in neighbors: # for each neighbor of each pixel
                value_n = self.gdalarray[neighbor.y, neighbor.x] # get the elevation value of the neighbor
                value_n_s = self.gdalarray_slope[neighbor.y, neighbor.x] # get the slope value of the neighbor
                if (value_n < value_p) & (value_n < max_elev) & (abs(value_p_s - value_n_s) < SLOPE_LIKE) & (value_n_s < MAX_SLOPE) & (neighbor not in run_set): # if the pixel matches all criteria
                    good_neighbors.add(neighbor) # add it to the set of good neighbors
                    run_set.add(neighbor) # add it to the run_set
                    good_neighbors_elev.add(value_n) # add it  to the elevation set

        max_elev = (k**(len(good_neighbors)/MAX_ELEV_DROP_REDUCER)) * max_elev # drop the max elevation

        if len(good_neighbors_elev) > 0:
            rel_max_elev = min(good_neighbors_elev) * R

        # run recursion
        self.findRunSwath(good_neighbors, run_set, k, R, max_elev, rel_max_elev, MAX_ELEV_DROP_REDUCER, MAX_SLOPE, SLOPE_LIKE)


    '''
    Finds the farthest pixel from the a given pixel in the run_set

    @param pixel is the comparison pixel
    @param run_set is the run_set
    '''
    def farthestSwathPixel(self, pixel, run_set):
        if len(run_set) < 1: # if the runsent is empty
            return False # don't bother
        # find farthest pixel
        farthest_pixel = pixel
        for swath_pixel in run_set:
            if swath_pixel.getEuclidDist(pixel) > farthest_pixel.getEuclidDist(pixel):
                farthest_pixel = swath_pixel
        return farthest_pixel

    '''
    Finds the 'best run' in a run_set; was originally implemented recursively; now is implemented using a while loop

    @param start_pixel is the starting point of the run; this parameter will remain constant
    @param current_pixel is the current pixel considered at the relevant step of in the loop; initially this should be set to the farthest pixel from the start
    @param best_run_set is the set of pixels in the best_run
    @param FLOWY is a boolean: if FLOWY = True, a run with similar slope values will be computed; if FLOWY = False, a run with highly varied slope values will be computed; this parameter will remain constant
    @param WINDY is a boolean: if WINDY = True, a run which moves slowly back towards the start pixel will be computed; if WINDY = Fales, a run which heads directly back towards the start pixel will be computed; this parameter will remains constant
    @param TYPE_WEIGHT allows the user to modulate the emphasis on FLOWINESS or WINDINESS; should be a value between 0 and 1
    @param STEEPNESS_WEIGHT allows the user to modulate the steepness of the best_run

    @returns the best run as a list and produces shapefile output
    '''
    def findBestRun(self, START_PIXEL, RUN_SET, FLOWY, WINDY, TYPE_WEIGHT, STEEPNESS_WEIGHT):

        # this is a back-link method START_PIXEL actually refers to the endpoint of this algorithm
        # START_PIXEL, is, however, the start point of the run

        try: # check to see that the raster has the needed properties
            self.dem_meta
            self.gdalraster
            self.gdalarray
            self.gdalraster_slope
            self.gdalarray_slope
        except: # ends the function and returns false
            return False

        if len(RUN_SET) < 1:
            return False

        best_run_set = set()

        farthest_pixel = self.farthestSwathPixel(START_PIXEL, RUN_SET)
        current_pixel = farthest_pixel

        while current_pixel != START_PIXEL: # so long as we have not reached the start pixel

            # add the current pixel to the runset
            best_run_set.add(current_pixel)

            neighbors = (current_pixel.neighborsAsPixels(1) & RUN_SET) - best_run_set # we only want numbers that are in the run_set, but that have not yet been added to the best_run_set
            best_neighbor = None

            if len(neighbors) < 1: # you may have reached a deadend, try to find a way out
                get_out_neighbors = current_pixel.neighborsAsPixels(1) & RUN_SET
                if len(get_out_neighbors) < 2: # if you've only got one way out, use it
                    best_neighbor = list(get_out_neighbors)[0]
                else: # if you've got options, pick the one closest to the start pixel, we could do this using set min or max, but computational efficiency gain would be minimal on account of the size limit of the set of neighbors at 6
                    smallest_distance = None
                    closest_pixel = None
                    for swath_pixel in get_out_neighbors:
                        distance = swath_pixel.getEuclidDist(START_PIXEL)
                        if (smallest_distance is None) or (distance < smallest_distance):
                            smallest_distance = distance
                            closest_pixel = swath_pixel
                    best_neighbor = closest_pixel # your best neighbor is the one closest to the start point
            elif len(neighbors) < 2: # if you you've only got one neighbor, that's your best neighbor
                best_neighbor = list(neighbors)[0] # here?
            else: # find best neighbor
                value_p = self.gdalarray[current_pixel.y, current_pixel.x] # get the elevation value of the current pixel
                value_p_s = self.gdalarray_slope[current_pixel.y, current_pixel.x] # get the slope value of the current pixel
                value_p_d = current_pixel.getEuclidDist(START_PIXEL) # get the distance to the start pixel for the current pixel
                highest_neighbor_score = None # will detail the score below
                for neighbor in neighbors:
                    value_n = self.gdalarray[neighbor.y, neighbor.x] # get the elevation value of the neighbor
                    value_n_s = self.gdalarray_slope[neighbor.y, neighbor.x] # get the slope value of the neighbor
                    value_n_d = neighbor.getEuclidDist(START_PIXEL) # get the distance to the initial pixel of the destination pixel

                    # Calculate the score
                    # We'll use a metric wherein the user can select between flowly or gnarly (slope-related) and windy or fast (distance-related)
                    slope_score = 0
                    distance_score = 0
                    if FLOWY:# If the user wants a flowy run, assign a better score to the pixel the lower the difference in slope
                        slope_score = -1*abs(value_n_s - value_p_s)/(value_p_s +.000001) # .00001 for division by zero errors
                    else: # If the user wants a gnarly run, assign a better score to the pixel the higher the difference in slope 
                        slope_score = abs(value_n_s - value_p_s)/(value_p_s + .000001) 

                    if WINDY: # if the user wants a windy run, give better scores to pixels at a greater distance from the start
                        distance_score = value_n_d/(value_p_d + .000001)
                    else: # if the user wants  a fast run, giver better scores to pixels at a lesser distance distance from the start
                        distance_score = value_n_d/(value_p_d + .000001)

                    # compute final score; the last bit is to force the algorithm to run uphill to degree specified by user input
                    neighbor_score = (slope_score * TYPE_WEIGHT) + (distance_score * (1 - TYPE_WEIGHT)) * ((value_n/value_p) * STEEPNESS_WEIGHT)

                    if (highest_neighbor_score is None) or (neighbor_score > highest_neighbor_score):
                        highest_neighbor_score = neighbor_score
                        best_neighbor = neighbor

               
            current_pixel = best_neighbor
        else:
            best_run_set.add(START_PIXEL) # add the start pixel once you've made it back
        


        # we need to refine the best_run_set down to a one-pixel-thick ordered set
        # I'm going to do this solely by comparing elevation which should refine things
        # The thing we really want is for this to run pretty much just downhill as much as possible, i.e., uphill in the context of the algorithm
        refined_run_set = set()
        refined_run_list = list() # we'll need to use a list so we can preserve the order for a return value

        run = ogr.Geometry(ogr.wkbLineString) # initialize the run wkt geometry

        current_pixel = farthest_pixel
        while current_pixel != START_PIXEL: # so long as we're not back to the start

            refined_run_set.add(current_pixel) # add the current pixel to the refined_run_set...
            refined_run_list.append(current_pixel) # and the list (to preserve order)...
            lng, lat = self.lngLatFromPixel(current_pixel) # and the ogr Geometry
            run.AddPoint(lng, lat)

            neighbors = (current_pixel.neighborsAsPixels(1) & best_run_set) - refined_run_set # all neighbors in best_run_set and not in the refined_run_set

            best_neighbor = list(neighbors)[0]
            best_elevation = self.getValue(best_neighbor)

            for neighbor in neighbors:
                value_n = self.gdalarray[neighbor.y, neighbor.x] # get the elevation value of the neighbor
                elev_dif = value_n

                if (    (len((best_neighbor.neighborsAsPixels(1) & best_run_set) - refined_run_set) < 1 and len((neighbor.neighborsAsPixels(1) & best_run_set) - refined_run_set) > 0) # if the current best_neighbor doesn't have neighbors and the considered neighbor has good neighbors of its own to move to...
                        or (value_n > best_elevation and len((neighbor.neighborsAsPixels(1) & best_run_set) - refined_run_set) > 0) # or if the considered neighbor is higher in elevation and has good neighbors to move to...
                        or (current_pixel == START_PIXEL) # or if we'ved made it back to the start
                    ):
                    best_neighbor  = neighbor
                    best_elevation = value_n

            current_pixel = best_neighbor
        else: # add the START_PIXEL on loop exit
            refined_run_set.add(START_PIXEL)
            refined_run_list.append(START_PIXEL)
            lng, lat = self.lngLatFromPixel(START_PIXEL)
            run.AddPoint(lng, lat)

        # write to run to wkt
        with open(self.dem_meta.output_dir + "/run.wkt", 'w') as file:
            file.write(run.ExportToWkt())

        # output wkt to shape
        sref = osr.SpatialReference()
        sref.ImportFromEPSG(RUN_FINDER_EPSG)

        driver = ogr.GetDriverByName("ESRI Shapefile")
        output = driver.CreateDataSource(self.dem_meta.output_dir + "/run.shp")

        output_layer = output.CreateLayer("run", sref, geom_type=ogr.wkbLineString)

        # Add attribute fields
        f_defn = ogr.FieldDefn("ID", ogr.OFTInteger)
        f_defn.SetWidth(10)
        output_layer.CreateField(f_defn)

        with open(self.dem_meta.output_dir + '/run.wkt') as file:
            for id, row in enumerate(file.readlines()): # we'll do all lines of the file just in case anything weird happens
                line = ogr.CreateGeometryFromWkt(row)
                feature = ogr.Feature(output_layer.GetLayerDefn())
                feature.SetGeometry(line)
                feature.SetField("ID", id) # A field with an unique id.
                output_layer.CreateFeature(feature)
                feature.Destroy()
        output.Destroy()

        with open(self.dem_meta.output_dir + '/run.wkt') as wkt:
            # prepare display output
            source = osr.SpatialReference()
            source.ImportFromEPSG(3857)

            target = osr.SpatialReference()
            target.ImportFromEPSG(4326)

            transform = osr.CoordinateTransformation(source, target)
            display_wkt = ogr.CreateGeometryFromWkt(wkt)
            display_wkt.Transform(transform)
            with open(self.dem_meta.output_dir + '/display_run.wkt', 'w') as display:
                display.write(display_wkt.ExportToWkt())

        return refined_run_list # for legacy purposes


    '''
    Produce final raster

    @run_set is the set of run pixels which will be converted into a raster
    '''
    def setToRaster(self, run_set, name):
        # create an array with dimensions of image
        arr = np.zeros([self.gdalraster.RasterYSize, self.gdalraster.RasterXSize], np.float64)
        
        for pixel in run_set:
            arr[pixel.y, pixel.x] = pixel.value

        # prepare target crs
        target = osr.SpatialReference()
        target.ImportFromEPSG(RUN_FINDER_EPSG)

        # set driver
        driver = gdal.GetDriverByName('GTiff')
        filepath = output_raster = self.dem_meta.output_dir + "/" + name + "_dem.elevation.tif"
        output = driver.Create(filepath, self.gdalraster.RasterXSize, self.gdalraster.RasterYSize, 1, gdal.GDT_Float64)

        # write to tiff
        output.SetGeoTransform(tuple(self.gdalraster.GetGeoTransform()))
        output.SetProjection(target.ExportToWkt())
        output.GetRasterBand(1).WriteArray(arr)
        output.GetRasterBand(1).SetNoDataValue(-9999)
        output = None

        if name == "swath": # if the raster computed represents the swath
            self.gdalraster_run_swath = gdal.Open(filepath)
            self.gdalarray_run_swath = gdal_array.LoadFile(filepath)
            return self.gdalraster_run_swath, self.gdalarray_run_swath
        elif name == "best": # if the raster computed represents the best run; (built for legacy)
            self.gdalraster_best_run = gdal.Open(filepath)
            self.gdalarray_best_run = gdal_array.LoadFile(filepath)
            return self.gdalraster_best_run, self.gdalarray_best_run
        else:
            return gdal.Open(self.run_filepath), gdal_array.LoadFile(self.run_filepath) 

 


###############################################################################################################
# CLASS: PixelLite
# a lighweight version of a pixel for quick computation 
class PixelLite:
    '''
    Initializes the pixel using a pixel x-value, a pixel y-value, a value for the pixel and the link back to a raster meta
    '''
    def __init__(self, x, y, value, raster_meta):
        # set x, y and value
        self.x = x
        self.y = y
        self.value = value

        # set max_x and max_y
        self.raster_meta = raster_meta
        self.max_x = raster_meta.x_max
        self.max_y = raster_meta.y_max

        # set error flag
        self.flag = False

        # the set of all neighbors
        self.neighbors = set()

        if self.x > self.max_x | self.x < 0 | self.y > self.max_y | self.y < 0: # something has gone terribly wrong
            self.flag = True
            return

        # neighbor columns
        if self.x > 0:
            left = self.x - 1
        else:
            left = self.x

        center = self.x

        if self.x < self.max_x:
            right = self.x + 1
        else:
            right = self.x

        # neighbor columns
        if self.y > 0:
            bottom = self.y - 1
        else:
            bottom = self.y

        middle = self.y

        if self.y < self.max_y:
            top = self.y + 1
        else:
            top = self.y

        # top row of neighbors
        self.neighbors.add((left, top))
        self.neighbors.add((center, top))
        self.neighbors.add((right, top))

        # middle row of neigbors
        self.neighbors.add((left, middle))
        self.neighbors.add((right, middle))

        # bottom row of neighbors
        self.neighbors.add((left, bottom))
        self.neighbors.add((center, bottom))
        self.neighbors.add((right, bottom))

        # the method I'm using will sometimes add the center, middle element as I've used that as a default to set values to if there is an out of bounds error. So let's just pop that off...
        self.neighbors.discard((center, middle))

    '''
    Gets all neighboring pixels as objects of type PixelLite
    Note: PixelLite objects cannot be initialized with pixels as neighbors or else all neighbors for a raster will be generated at once
    This method, further, provides additional control over pixel value

    @param args[0] will be the value with which all pixels are initialized, if desired

    @returns the set of neighboring pixels
    '''
    def neighborsAsPixels(self, *args):
        if len(args) < 1: # if a value has not been provided to initialize the neighbors with
            neighbors_as_pixels = set()
            for pixel in self.neighbors:
                neighbors_as_pixels.add(PixelLite(pixel[0], pixel[1], self.value, self.raster_meta)) # initialize the neighbors with this pixel's value
            return neighbors_as_pixels
        else: # if a value has been provided to initialize the neighbors with
            neighbors_as_pixels = set()
            for pixel in self.neighbors:
                neighbors_as_pixels.add(PixelLite(pixel[0], pixel[1], args[0], self.raster_meta)) # initialize the neighbors with said value
            return neighbors_as_pixels

    '''
    Gets the Euclidean distance to another pixel

    @param pixel is the comparison pixel

    @returns the Euclidean distance to the comparison pixel
    '''
    def getEuclidDist(self, pixel):
        return math.sqrt((self.x - pixel.x)**2 + (self.y - pixel.y)**2)

    '''
    Print function for PixelLite
    '''
    def __str__(self):
        return "( {} , {} )".format(self.x, self.y)

    '''
    Reproduce function for PixelLite
    '''
    def __repr__(self):
        return "( {} , {}, {} )".format(self.x, self.y, self.raster_meta)

    '''
    Hash function for PixelLite
    '''
    def __hash__(self):
        return hash(repr(self))
    '''
    Comparison function for PixelLite
    '''
    def __eq__(self, other):
        return repr(self) == repr(other)
    '''
    Associates a raster with the PixelLite object

    @param raster is the raster to associate
    '''
    def assocRaster(self, raster):
        self.raster = raster
    
    '''
    Associates a gdalarray with the PixelLite object
    '''
    def assocGdalArray(self, g_array):
        self.g_array = g_array

###############################################################################################################
# CLASS: UserSettings
# a class for storing user settings
class UserSettings:

    '''
    Initializes a user settings object

    @params k, R, MAX_ELEV_DROP_REDUCER, MAX_SLOPE and SLOPE_LIKE are relevant to the @func DEM.findRunSwath()
    @params FLOWY, WINDY, TYPE_WEIGHT and STEEPNESS_WEIGHT are all relevant to the @func DEM.findBestRun()
    '''
    def __init__(self, *, k=.99999, R=1.2, MAX_ELEV_DROP_REDUCER=10, MAX_SLOPE=45, SLOPE_LIKE=7, FLOWY=False, WINDY=False, TYPE_WEIGHT=.5, STEEPNESS_WEIGHT = 1.0):
        
        self.k = k
        self.R = R
        self.MAX_ELEV_DROP_REDUCER = MAX_ELEV_DROP_REDUCER
        self.MAX_SLOPE = MAX_SLOPE
        self.SLOPE_LIKE = SLOPE_LIKE

        self.FLOWY = FLOWY
        self.WINDY = WINDY
        self.TYPE_WEIGHT = TYPE_WEIGHT
        self.STEEPNESS_WEIGHT = STEEPNESS_WEIGHT



###################################################################################################################
###################################################################################################################
###################################################################################################################
# CORE: findRun
# the core algorithm to be called by module user

'''
Finds a run given a start point and a set of bounds

@param start_point is the start point lng, lat pair
@param bounds are the bounds of the area to be looked at for a run
@param user_settings are the settings provided by the user to modulate the algorithm

@returns the output directory
'''
def findRun(start_point, bounds, *, user_settings=UserSettings()):

    
    OUTPUT_DIR = tempfile.mkdtemp(dir='static/output/')

    #ensure that the output dir has write access
    os.chmod(OUTPUT_DIR, 0o777)

    dem = DEM(GoogleDEM(bounds).output_dir)

    dem.addSlope()

    px, py = dem.pixelFromLngLat(start_point[0], start_point[1])

    pixel = PixelLite(px, py, 1, dem.dem_meta)
    lng, lat = dem.lngLatFromPixel(pixel)
    pixels = set([pixel, ])
    run_swath = set([pixel, ])
    best_run = set()

    dem.findRunSwath(pixels, run_swath, 
                     user_settings.k, user_settings.R, dem.getValue(pixel), dem.getValue(pixel), 
                     user_settings.MAX_ELEV_DROP_REDUCER, user_settings.MAX_SLOPE, user_settings.SLOPE_LIKE)

    dem.findBestRun(pixel,
                    run_swath, 
                    user_settings.FLOWY, user_settings.WINDY, 
                    user_settings.TYPE_WEIGHT, user_settings.STEEPNESS_WEIGHT)

    dem.setToRaster(run_swath, "swath")

    shutil.make_archive(OUTPUT_DIR + '-output', 'zip', OUTPUT_DIR)

    return OUTPUT_DIR

