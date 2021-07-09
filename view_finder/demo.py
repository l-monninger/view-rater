##############
#Title: tifview Demo
#Author: Liam Monninger SID 004748349
#Purpose: Demos tifview module (final project).
##############
try: # check that tifview is an accessible module
    import tifview as tv
except:
    print("It doesn't seem like you have tifview in the same folder as this script.")

print("Hello and welcome to TIFVIEW, my final project!\n"
      "I'll get into the demos shortly. But, first a brief note on the development of this project:\n\n"
      "Originally, I had planned on developing a module that would use the requests API and the service\n"
      "put together by Dr. Burkhart to build a DEM. I thought this might be slow, but I didn't realize\n"
      "just how slow...\n")

see_dem = input("Press 'Y' if you'd like to see a small example of this process: ")

if see_dem == 'Y' or see_dem == 'y':
    bounds = tv.Bounds(40.02, 40.02, 40.03, 40.03) # initialize bounds

    try:
        tv.DEM(bounds, 'test')
        print("Sample DEM output to current working directory as test.tif.elevation.\n"
              "Check out how small the DEM is. And, it took that long!")
    except:
        print("Sample DEM failed. Ensure your working directory has write permissions.")

print("\nI briefly considered using grequests to do this as quickly as possible. But, in the end, I decided\n"
      "to do something DEM-related which (I hope) is a bit more interesting.\n\n"
      "Thus, without further ado, allow me to introduce: How's my view?--the quick and easy view rater. Simply\n"
      "provide a longitude and latitude coordinate pair, and I'll tell give your spot a rating.\n"
      "(Note: As you'll see in the comments in my module, the scores here are rather tricky to normalize.\n"
      "Generally speaking, they seem to fall in the range -3 to 3.)")

print("\nHere's a comparison test between a point on the top of Mount Shasta and a point on the floor of\n" 
    "California's North Valley:")
try:
    print("Mt. Shasta Score: ", tv.rateView(-122.195105, 41.409213, 6))
    print("Valley Floor Score: ", tv.rateView(-122.328904, 39.446331, 6))
    print("\nOk, so gimmicks aside, the cool thing about this function is, that you cold potentially produce a raster\n"
      "which has pixels with values representing the quality of view, i.e., this could be a pretty nifty little site\n"
      "suitability algorithm.")
except:
    print("An error occurred. Please check your internet connection.")
