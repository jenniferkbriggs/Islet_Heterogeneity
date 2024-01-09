# %% 
# This code is made to extract the functional network from Calcium timecourses. 
# Jennifer Briggs 2022

## ---- Options for you to change -------
# set to true if you'd like to see and save figures. set to false if you don't need figures
fig_on = True

# if you want to predefine a savepath. If not, comment out this line by putting at # in front!
global savepath
#savepath = '/Users/briggjen/Documents/GitHub/Islet_Analysis/Examples/SignalProcessing/Slow_'
path = '/Users/briggjen/Documents/GitHub/Islete_Heterogeneity/Examples/Calcium.csv'

# savepath = 'E:/Promotion_Postdoc/CMOS/CMOS_Daten/210705_Aktivitaet/HDF5_3985'
# path = 'E:/Promotion_Postdoc/CMOS/CMOS_Daten/210705_Aktivitaet/HDF5_3985/210705_3985_G8_I1.h5'

# How do you want to define the threshold? (Either number_of_connections, scalefreeish, setthreshold)
threshold_opts = 'number_of_connections'

#for threshold_opts = 'setthreshold'
threshold_set = 0.7 #change if you choose to set the threshold manually

#for threshold_opts = 'number_of_connections'
k = 6

#threshold_opts = 'scalefreeish'
min_connect = 3 #minimum average connections for the scale free threshold
max_connect = 12


USE_CONFIGURED_ISLETS = 'FALSE'
TEST_FILE = ''



# %% Import packages
import csv
from collections import OrderedDict
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.optimize as op
#import easygui #for selecting files using gui
import pandas as np
import math
import tkinter as tk
#from tkinter import ttk
#from tkinter.messagebox import askyesno
from Extract_Functional_Net import *
from LoadData import *

# %%  Load timeseries file
if 'path' in locals():
    ca = LoadData(path, USE_CONFIGURED_ISLETS, TEST_FILE)
else:
    ca = LoadData()

try: #if time is in the first axis, we save it and remove
    time = ca.Time
    ca = ca.drop('Time', axis=1)
except:
    print('No time avaliable')
    timeopt = "No"


# %% Compute the correlation matrix
cor_mat = ca.corr() #computes correlation matrix

if fig_on:
    f = plt.figure(figsize=(19, 15))
    plt.matshow(cor_mat, fignum=f.number)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    plt.title('Correlation Matrix', fontsize=16)
    if 'savepath' not in locals():
        print('Select folder to save data in')
        savepath = easygui.diropenbox('Select Folder to Save Data In')
        savepath = savepath + savepath[0] #adds slash
    try: 
        plt.savefig(savepath + 'Corrmat.png')
    except:
        print('Save path is not working')
        savepath = easygui.diropenbox('Select Folder to Save Data In')
        savepath = savepath + savepath[0] #adds slash
        plt.savefig(savepath + 'Corrmat.png')


    plt.clf

cor_active = [np.mean(cor)>0.1 for cor in cor_mat.values]
cor_mat = cor_mat.iloc[cor_active, cor_active]

# set diagonals equal to zero:
cor_mat = cor_mat.where(cor_mat.values != np.diag(cor_mat),0,cor_mat.where(cor_mat.values != np.flipud(cor_mat).diagonal(0),0,inplace=True))

# %% Computing the network -- need to code in how to find the threshold (8 or power law)
# NOT WORKING
#If how to set threshold is not predefined, choose how to set through gui 
# if 'threshold_opts' not in locals():
#     root = tk.Tk()
#     # click event handler
#     def b_degree():
#         threshold_opts = 'number_of_connections'
#         min_connect = 5
#         max_connect = 20
#         print('done')
#         root.destroy()
#         return threshold_opts
#     def b_scalefree():
#         threshold_opts = 'scalefreeish'
#         print('done')
#         root.destroy()
#         return threshold_opts
#     top = ttk.Frame(root)
#     bottom = ttk.Frame(root)
#     top.pack(side=tk.TOP)
#     bottom.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
#     # create the widgets for the top part of the GUI,
#     # and lay them out
#     b = ttk.Button(root, text="Predefined average degre", command=b_degree)
#     c = ttk.Button(root, text="Scale Free (ish)",  command=b_scalefree)
#     b.pack(in_=top, side=tk.LEFT)
#     c.pack(in_=top, side=tk.LEFT)

#     # start the app
#     root.mainloop()


# %%
if threshold_opts == 'number_of_connections':
    # Specify average number of connections:
    thr = thr_based_on_degree(cor_mat, k)
elif threshold_opts == 'scalefreeish':
    maxbnds = float(thr_based_on_degree(cor_mat, min_connect)) #because the minimum connection gives the largest threshold
    minbnds = float(thr_based_on_degree(cor_mat, max_connect))
    bnds = (minbnds, maxbnds)
    x0 = np.mean(bnds)
    # Run optimization
    final_symp = op.minimize(makegraph_err,x0, args=(cor_mat), method = 'Nelder-Mead', bounds = ((minbnds, maxbnds),))
    #find treshold
    thr = final_symp.x[0]
elif threshold_opts == 'setthreshold':
    thr = threshold_set

G = makegraph(thr, cor_mat)


if fig_on: 
    from pyvis.network import Network
    net = Network(notebook = True)
    net.from_nx(G)
    net.show(savepath + "network.html")

    if type(list(G.nodes)[0]) is str: #node names are the positions
        positions = {i:list(map(int, i.split(","))) for i in G.nodes}

    deg_seq = lookatnetwork(G, savepath, positions)


#save adjacency list
nx.write_adjlist(G, savepath + "Network.adjlist")
# %% Extract network stats

# Find hubs (which have degree greater than 60%)
degrees = [deg for (node, deg) in G.degree()]
sixtypercentdegree = round(max(degrees)*0.6) 
deg = dict(G.degree)
hubs = {k:v for (k,v) in deg.items() if v >= sixtypercentdegree}

# Centrality hubs
deg = (nx.degree_centrality(G))
sixtypercentdegree = (max(list(deg.values()))*0.6) 
centralhubs = {k:v for (k,v) in deg.items() if v >= sixtypercentdegree}




net_stats = {'Average_Correlation': np.mean(list(cor_mat.values)),
'Degree': 2*G.number_of_edges()/G.number_of_nodes(), 
'Clustering': nx.average_clustering(G), 
'Global_Efficiency': nx.global_efficiency(G),
'Local_efficiency': nx.local_efficiency(G),
'Hubs': list(hubs.keys()),
'Most_Central': list(centralhubs.keys()),
'Threshold': thr
}

with open(savepath + 'Net_stats.csv', 'w') as f:  
    writer = csv.writer(f)
    for k, v in net_stats.items():
       writer.writerow([k, v])
# %%
