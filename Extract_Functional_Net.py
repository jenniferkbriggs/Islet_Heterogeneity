#Extract_Functional_Net.py>
# %%
# This code contains functions to run Run_Network.py
# Jennifer Briggs 2022
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import scipy.optimize as op
import math

# %% Functions and Classes
def lookatnetwork(G, savepath, positions):
    degree_sequence = sorted((d for n, d in G.degree()), reverse=True)
    dmax = max(degree_sequence)

    fig = plt.figure("Degree of a random graph", figsize=(8, 8))
    # Create a gridspec for adding subplots of different sizes
    axgrid = fig.add_gridspec(5, 4)

    ax0 = fig.add_subplot(axgrid[0:3, :])
    Gcc = G.subgraph(nx.draw(G, pos=positions, node_size = 15))

    ax0.set_title("Visualization of network")
    ax0.set_axis_off()

    ax1 = fig.add_subplot(axgrid[3:, :2])
    ax1.plot(degree_sequence, "b-", marker="o")
    ax1.set_title("Degree Rank Plot")
    ax1.set_ylabel("Degree")
    ax1.set_xlabel("Rank")

    ax2 = fig.add_subplot(axgrid[3:, 2:])
    ax2.bar(*np.unique(degree_sequence, return_counts=True))
    ax2.set_title("Degree histogram")
    ax2.set_xlabel("Degree")
    ax2.set_ylabel("# of Nodes")

    fig.tight_layout()
    plt.savefig(savepath + 'NetworkStats.png', facecolor=fig.get_facecolor(), edgecolor='none')

    return degree_sequence #export degree sequence



def makegraph(thr, cor_mat):
    #Finding optimal threshold based on log log fit
    G = nx.Graph() #make an empty graph
    G.add_nodes_from(list(cor_mat)) #add node for each cell
    #loop over every combination of nodes and add an edge if the correlation is greater than threhsold
    for x in G.nodes:
        for y in G.nodes:
            if cor_mat.loc[x,y] > thr and x != y and G.has_edge(x, y) == False:
                G.add_edge(x,y)

    return G


## To compute thresholds
def loglogcost(x, y):
    # remove any zeros - THIS STEPS ASSUMES THE ONLY THE CONNECTED NODES FIT A POWER LAW
    slope, intercept, r, p, std_err = stats.linregress(np.log10(x), np.log10(y))
    err = (r + 1) #r should be -1 so error is r - (-1)
    return err

def makegraph_err(thr, cor_mat):
    #Finding optimal threshold based on log log fit
    G = nx.Graph() #make an empty graph
    G.add_nodes_from(list(cor_mat)) #add node for each cell
    #loop over every combination of nodes and add an edge if the correlation is greater than threhsold
    for x in G.nodes:
        for y in G.nodes:
            if cor_mat.loc[x,y] > thr and x != y and G.has_edge(x, y) == False:
                G.add_edge(x,y)

    deg_dict = dict(G.degree()) #export degree of each node into a dictonary
    deg_hist = np.histogram(list(deg_dict.values()))
    deg_hist_bins = list(deg_hist[1])
    deg_hist_vals = list(deg_hist[0])

    x_omit = [i for i,x in enumerate(deg_hist_vals) if x == 0 ]

    for i in np.flip(x_omit):
        deg_hist_bins.pop(i)
        deg_hist_vals.pop(i)

    err = loglogcost(deg_hist_bins[1:],deg_hist_vals)

    return (err)


def thr_based_on_degree(cor_mat, k):
    #This function finds the threshold for a corrlation matrix (cor_mat) that gives the network a desired average degree k.
    #k = desired average degree
    #m = k*n/2 : where m is the total number of edges and n is the total number of nodes
    #p = 2m/(n(n-1)): the percent of edges we want compared to the number of possible edges
    #Once p is found, we sort all non-self correlations and find the threshold which returns p percent.

    p = 2*k/(len(cor_mat)-1)
    all_cors = np.reshape(list(cor_mat.values), (1,-1)) #1d list of all correlations
    all_cors_sort = sorted(all_cors[0], reverse=True)
    m = (k*len(cor_mat)) #deisred number of edges
    last_val = p*len(cor_mat) #value to pick for threshold
    thr = all_cors_sort[math.ceil(m)]

    return thr



# %%
