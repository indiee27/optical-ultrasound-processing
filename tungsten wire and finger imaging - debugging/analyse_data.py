# script for analysing image data combining cross-talk removal and image reconstruction

# %%

# imports
import os

# import data and set parameters
directory = 'C:/Users/indie/Documents/GitHub/optical-ultrasound-processing/tungsten wire and finger imaging - debugging/'

file_list = list()

for fname in os.listdir(directory):
    if fname.endswith('.txt'):
        file_list.append(fname)
        continue
    else:
        continue

print(file_list)

# %%

def import_data(file):
    fileObj = open(file, "r")
    data = fileObj.read().splitlines()
    return data

# %%

# set constants
Nx = 200
Nt = data.shape[0] # this needs to be inside the loop so that data is quantified