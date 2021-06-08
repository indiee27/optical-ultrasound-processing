# imports
from utils import all #import the functions from the utils script

# name the directory
directory = ("C:/Users/indie/Documents/GitHub/optical-ultrasound-processing/tungsten wire and finger imaging - debugging")

# list the file prefixes
file_list = file_listing_func(directory)

# main loop
for item_prefix in file_list:

    # hydrophone parameters (dictionary) using the function from utils
    hydrophone_parameters = params_import(item_prefix)

    # import list of image files
    images = data_import(directory, item_prefix)

    # second loop
    for file in images:
        data = read_file(file)
