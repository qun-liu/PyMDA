#!/bin/env python
# free to use or modify
#  Authors: 
#  Qun Liu: qun.liu@gmail.com
#  Lina Takamaru:  llt45@cornell.edu
#  References:  J. Appl. Cryst. (2020). 53 https://doi.org/10.1107/S160057671901673X
from scipy.cluster import vq
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree, fcluster
from matplotlib.pyplot import *
import sys
from numpy import genfromtxt
from scipy import io
import os

def llf(id):
    return str(id+1)

def get_cluster_elements(clust):
    # return ids for elements in a cluster sub-tree
    #if clust.id>0:
    if clust.is_leaf():
        # positive id means that this is a leaf
        return [clust.get_id()]
    else:
        # check the right and left branches
        cl = []
        cr = []
        if clust.get_left()!=None:
            cl = get_cluster_elements(clust.get_left())
        if clust.get_right()!=None:
            cr = get_cluster_elements(clust.get_right())
        return cl+cr

def ucr(data, nClusters, dir_output):
	data = genfromtxt(data,comments='#', delimiter =",")
    	n_clusters = int(nClusters)
    	nsets =len(data) - 1
    	li = linkage(data,method='ward', metric='euclidean')
    	tree, mapper = to_tree(li,rd=True)
    	clusters=fcluster(li,n_clusters,criterion='maxclust')
    	
	if os.path.exists(dir_output + '/cluster.txt'):
		os.remove(dir_output + '/cluster.txt')
	
	clust = open(dir_output+"/cluster.txt", "a")
	for i_cluster in range(1,n_clusters+1):
        	for leaves in range(nsets+1):
          	  if clusters[leaves] == i_cluster:
                	print "cluster" + str(i_cluster) + "q",  leaves+1
			clust.write("cluster" + str(i_cluster) + "q " + str(leaves+1) + "\n")
	clust.close()
	
    	test=dendrogram(li,color_threshold=200,orientation='right', leaf_label_func=llf)
    	#xlim((5, 0))
    	xmin,xmax = xlim()
    	matplotlib.rcParams.update({'lines.color': '#112233'})
    	ax = gca()
    	for line in ax.xaxis.get_ticklines():
       	 	line.set_markeredgewidth(2)
    	savefig(dir_output+'/unit_cell_single.png')#,transparent=True)
	os.rename(dir_output+'/unit_cell_single.png', dir_output + '/c_unit_cell_single.png')	
	
    	#show()
