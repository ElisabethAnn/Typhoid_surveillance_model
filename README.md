# Typhoid_surveillance_model

Files needed:
Download shapefiles of drainage lines, drainage points, and, if available, blue lines from the community. es.world provides these shapefiles for select cities. If no blue lines or drainage lines are available one option is to use street maps from osmdata but a new spatial script will need to be created. The build_drainage_map_final.Rmd relies on the drainage points provided by es.world shapefiles to determine the direction of flow. 

See https://es.world/help/content/methodology/10_main.html for methodology on how drainage lines, drainage points, and flow direction are determined for Novel-t's maps. 

Put shapefiles, .Rmd, and both .py scripts into one directory. 

Script 1: build_drainage_map_final.Rmd
If drainage lines, drainage points, and blue lines are available, run build_drainage_map_final.Rmd to create a system map and save a .csv of the input table required for the model. 

Script 2: detectProb.py
This script contains the function for calculating detection probability. It should be stored in the working directory when running run_model.py.

Script 3: dyn_model.py
This script contains all the modules required for the main model. It should be storedin the working directory when running run_model.py. 
