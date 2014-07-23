import sys, os, time
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib


output_path = os.getcwd().replace("Scripts","Output/")

f_brca = open(output_path+"brca/normal_rnd_betas_1.txt","r").readlines()
f_coad = open(output_path+"coad/normal_rnd_betas_1.txt","r").readlines()

barcodes = []
row_num = len(f_brca)+len(f_coad)-2
betas = np.zeros([row_num,27578])
i = 0
for line in f_brca:
	line = line.strip().split('\t')
	if line[0] == "Samples":
		continue
	else:
		barcodes.append(line[0])
		betas[i,:] = np.array(line[1:])
		i += 1
for line in f_coad:
	line = line.strip().split('\t')
	if line[0] == "Samples":
		continue
	else:
		barcodes.append(line[0])
		betas[i,:] = np.array(line[1:])
		i += 1
print betas

labels_d = []
for i in range(0,len(f_brca)):
	labels_d.append("B"+str(i))
for i in range(0,len(f_coad)):
	labels_d.append("C"+str(i))


betas_dist = pdist(betas,"euclidean")
print betas_dist
betas_linkage = linkage(betas_dist,method="average")
dendrogram(betas_linkage,labels=labels_d)
matplotlib.pylab.show()
#matplotlib.pylab.savefig("test.png")



