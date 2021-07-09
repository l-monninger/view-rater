##############
#Title: Finding City Elevations Using Requests
#Author: Liam Monninger SID 004748349
#Purpose: Practicing working with API's and producing usable data
##############
from osgeo import ogr
import csv
import requests

filepath = "tl_2017_06_place/tl_2017_06_place.shp" # Identify the location of the shapefile to be used

driver = ogr.GetDriverByName("ESRI Shapefile") # Create an OGR Driver instance, with the driver type "ESRI Shapefile"
shapefile = driver.Open(filepath, 0) # Open the shapefile using the driver (read-only)
layer = shapefile.GetLayer() # Get the data layer from the shapefile (there should be only one with shapefiles, so this is simple)

with open("city_elevations.csv", "w", newline='') as output: # Open a CSV file for writing with the DictWriter class
    fieldnames = ['name', 'lng', 'lat', 'elev_ft'] # set header names
    writer = csv.DictWriter(output, fieldnames=fieldnames) # initialize writer
    writer.writeheader() # write in header names
    for feature in layer: # Begin iterating over features in the layer
        name = feature.GetField('NAME')

        geom = feature.GetGeometryRef() # Get the geometry of the feature
        centroid = geom.Centroid() # Find the point centroid of the feature
        lng = centroid.GetX() # get longitude
        lat = centroid.GetY() # get latitude

        # make requests for elevation values
        r = requests.get('https://apps.gis.ucla.edu/elevation/api/v1/lookup?locations=%s,%s'%(lat, lng))
        result = r.json()
        elevation = result['results'][0]['elevation'] # get elevation
        elevation*=3.28084 # convert to feet

        print("{} ({},{}) -- {:.0f} ft".format(name, lat, lng, elevation)) # print out info
        writer.writerow({'name': name, 'lng': lng, 'lat': lat, 'elev_ft': elevation}) # write data to appropriate columns

