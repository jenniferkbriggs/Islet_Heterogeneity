#LoadData.py

# This is a function file which loads data. Here you can preset the path to give quick access to data or you can call a gui
# Jennifer Briggs - 2022. 

# Optional argument: path to timeseries file
# %%import packages
from Anne_code_spiketime_detection import OpenPrefined
import easygui
import pickle
import pandas as pd
import numpy as np
import h5py

def LoadData(path = 0, USE_CONFIGURED_ISLETS = 'False', TEST_FILE = ''):
    # if path is not passed in
    if path == 0:
        path = easygui.fileopenbox('Select Time signal file')

    #load file based on file type
    if path[-3:int(len(path))] == 'csv':
        ca = pd.read_csv(path, index_col=False)

    elif path[-3:int(len(path))] == 'mcr': #Anne's electrode data
        sd, td = OpenPrefined(path)
        #Electrode data is stored by x,y values. We will index them first along the x, then y. 
        # For example, [1,1] = 0, [1, 2] = 1, [2, 1] = len(y)
        numcell = (td.max_x - td.min_x+1)*(td.max_y - td.min_y+1)
        time = td.getElectrode(td.min_x, td.min_y)
        time = len(time.values)
        dat = np.empty([time, numcell])

        print('Reshaping Values')
        loc = np.empty([2,numcell])
        dat = dict()
        i = 0

        # if USE_CONFIGURED_ISLETS:
        #     islet_configuration = read_islet_configuration(TEST_FILE)
        #     for islet in islet_configuration:
        #         for e in islet:
        #             electrode = td.getElectrode(e[0],e[1])
        #             #dat[:,i] = list(electrode.values)
        #             dat.update({str(x) + ',' + str(y): list(electrode.values)})
        #             loc[:,i] = [x,y]
        #             i = i+1        
        # else:
        for y in range(td.min_y, td.max_y+1):
            for x in range(td.min_x, td.max_x+1):
                electrode = td.getElectrode(x,y)
                    #dat[:,i] = list(electrode.values)
                dat.update({str(x) + ',' + str(y): list(electrode.values)})
                loc[:,i] = [x,y]
                i = i+1
        
    elif path[-3:int(len(path))] == '.h5':
        with h5py.File(path, "r") as f: #this code is very specific for Anne Gresch's electrode data. Will need updated for another file
            # Print all root level object names (aka keys) 
            # these can be group or dataset names 
            print("Keys: %s" % f.keys())
            # get first object name/key; may or may NOT be a group
            a_group_key = list(f.keys())[0]

            # If a_group_key is a group name, 
            # this gets the object names in the group and returns as a list
            nextname = list(f[a_group_key])
            # preferred methods to get dataset values:
            nextname2 = list(f[a_group_key][nextname[0]])
            data = f[a_group_key][nextname[0]][nextname2[0]][()]
            metadata = f[a_group_key][nextname[0]][nextname2[1]][()]


        min_y = metadata['Region.Top'][0]
        min_x = metadata['Region.Left'][0]
        max_y = metadata['Region.Bottom'][0]
        max_x = metadata['Region.Right'][0]

        numcell = (max_x - min_x+1)*(max_y - min_y+1)

        loc = np.empty([2,numcell])
        dat = dict()
        i = 0

        for x in range(data.shape[1]):
            for y in range(data.shape[2]):

                #dat[:,i] = list(electrode.values)
                dat.update({str(min_x+x) + ',' + str(min_y+y): list(data[:,x,y])})
                loc[:,i] = [x+min_x,y+min_y]
                i = i+1


    else:
        raise Exception('File type not yet valid for this analysis')

        #fs = td.tickrate 

        #timeall = np.arange(0,int(time/fs),1/fs)
    timeall = np.arange(0, int(np.shape(data)[0]),1)
    ca = pd.DataFrame(dat)
    ca['Time'] = timeall

    print('Loaded Data Correctly')
    return ca 
# %%


# %%
