import os
import numpy as np
from scipy.stats.stats import pearsonr

os.chdir("/Users/Tuan/workspace/GenomeMDS/output/1mb/genome")

dist = np.loadtxt("distanceHeatmap.txt")

dist = dist[:24,:24] #remove the last col, row of 

#print dist.shape

contact = np.loadtxt("interContactMap.txt", delimiter = ",")
#print contact.shape

print pearsonr(dist.ravel(), contact.ravel())

#(-0.38603258124520456, 6.5788843823539546e-22)
