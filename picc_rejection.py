#!/bin/env python
# free to use or modify
#  Authors: 
#  Qun Liu: qun.liu@gmail.com
#  Lina Takamaru:  llt45@cornell.edu
#  References:  J. Appl. Cryst. (2020). 53 https://doi.org/10.1107/S160057671901673X

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from matplotlib.pyplot import *
import sys
from numpy import genfromtxt,fromstring
from scipy import io
import os

def llf(id):
    return str(id+1)

def pr(datafile, clusts, dir_output):
#if __name__ == "__main__":
	data1 = 1.0-fromstring(''.join(open(datafile, 'r').read().splitlines()),sep=" ")
	n_clusters = int(clusts)
	li1 = linkage(data1,method='average')
	clusters=fcluster(li1,n_clusters,criterion='maxclust')
	nsets = len(clusters)-1

	if os.path.exists(dir_output + '/cluster.txt'):
		os.remove(dir_output + '/cluster.txt')

	clusts = open(dir_output + "/cluster.txt", "a")
	for i_cluster in range(1,n_clusters+1):
		for leaves in range(nsets+1):
			if clusters[leaves] == i_cluster:
				print "clusterb" + str(i_cluster),  leaves+1
				clusts.write("clusterb" + str(i_cluster) + " " + str(leaves+1) + "\n")	
	clusts.close()
	#print("li1:", li1)
	fig = figure()	
	test1=dendrogram(li1,color_threshold=200,orientation='right', leaf_label_func=llf)
	ax1 = gca()
	for line in ax1.xaxis.get_ticklines():
		line.set_markeredgewidth(1)
	#axhline(y=0, xmin=0,xmax=1)
	#xlim((0.3, 0))
	xmin,xmax = xlim()
	locs,labels=xticks()
	#print "test", xmin, xmax, locs
	matplotlib.rcParams.update({'lines.color': 'black'})
	savefig('1-picc.pdf')#,transparent=True)
	fig.show()

if __name__ == "__main__":
	pr()
