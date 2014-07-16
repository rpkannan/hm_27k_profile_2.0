import sys, os, time
import numpy as np
import matplotlib.pyplot as plt


output_path = os.getcwd().replace("Scripts","Output/")

# Maps genes to probe ids: based on HM27K ADF
gs_to_probe_file = open(os.getcwd().replace("hm_27k_profile_2.0/Scripts","Annotation/Misc/gene_symbol_to_probe.txt"),"r").readlines()
gsp_dict = {}
for line in gs_to_probe_file:
	line = line.strip().split('\t')
	if line[0] == "Samples":
		continue
	else:
		gsp_dict[line[0]] = line[1:]
gene_symbols = sorted(gsp_dict.keys())

cancer = "gbm"
methyl_brc = open(output_path+cancer+'/methyl_barcodes.txt',"r").readlines()
f_1 = open(output_path+cancer+'/all_raw_betas.txt',"r").readlines()
f_2 = open(output_path+cancer+'/all_gene_ratios.txt',"r").readlines()

# For sample1, plot gene symbol ratio vs avg probe value
# X Value: gene_symbols
# Y Value: Avg. beta value

probes_27 = f_1[0].strip().split('\t')[1:]
gene_s = f_2[0].strip().split('\t')[1:]

betas = f_1[1].strip().split('\t')[1:]
ratios = f_2[1].strip().split('\t')[1:]


x_val = []
y_val = []
for g in gene_symbols:
	if g in gene_s:
		temp = []
		for p in gsp_dict[g]:
			b = betas[probes_27.index(p)]
			if b == "nan":
				continue
			else:
				temp.append(float(b))
		if len(temp) == 0:
			continue
		else:
			x_val.append(float(ratios[gene_s.index(g)]))
			y_val.append(sum(temp)/len(temp))
	else:
		x_val.append(0)
		y_val.append(0)

plt.scatter(x_val,y_val)
plt.savefig("test.png")

		
		
		
		
		
				
	