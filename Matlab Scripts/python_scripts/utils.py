# script for analysing image data combining cross-talk removal and image reconstruction
'''
script for analysing image data
'''
#%%
%matplotlib qt 

# %%

# imports
import os
import numpy as np 

#%%

# create prefix list from directory

def file_listing_func(directory):
    '''
    Import the list of file prefixes from the directory
    '''
    file_list = []
    
    for fname in os.listdir(directory):
        if fname.endswith('_hdr.txt'):
            file_list.append(fname)
        else:
            continue

    prefix_list = []
    for file in file_list:
        prefix_list.append(file.replace('_hdr.txt', ''))   

    return prefix_list

directory = ("C:/Users/indie/Documents/GitHub/optical-ultrasound-processing/tungsten wire and finger imaging - debugging")

prefix = (file_listing_func(directory))

#%%

# import data and set parameters
fileprefix = prefix[0]

def params_import(fileprefix):
    '''
    Import the parameters from the text files as a dictionary
    PARAMS: fileprefix: string location of file
    OUTPUT: param_dict: dictionary of key parameters
    '''
    filehdr = fileprefix + '_hdr.txt'
    
    f = open(filehdr, 'r')

    line_names = []
    for line in f:
        line = line.rstrip()
        line = line.split()
        line_names.append(line)

    key_list = []
    value_list = []
    for item in line_names:
        length = len(item)
        key = item[0]
        key_list.append(key)
        value = item[length-1]
        value_list.append(int(value))

    param_dict = dict(zip(key_list, value_list))

    return param_dict

# %%

# define function to read data
def read_file(filepath):
    '''
    function to read source file containing data
    '''
    data = np.genfromtxt(filepath)
    return data

# %%

# import list of images from prefix
def data_import(directory,prefix):
    '''
    Function to import all the txt files per prefix
    '''
    array_list = []
    for fname in os.listdir(directory):
        if fname.startswith(prefix + '_image') and fname.endswith('.txt'):
            array_list.append(fname)

    return array_list


# %%
test_prefix = prefix[0]
test1 = data_import(directory,test_prefix)
print(test1)
# %%

# defining constants
def constants():