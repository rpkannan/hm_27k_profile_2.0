import sys, os, time
import numpy as np
#from numpy.linalg import norm 
import pandas as pd 

'''
f_1 = open(os.getcwd().replace("Scripts","Output/brca/")+"normal_raw_betas.txt","r").readlines()
f_2 = open(os.getcwd().replace("Scripts","Output/brca/")+"normal_dis_betas.txt","r").readlines()
f_3 = open(os.getcwd().replace("Scripts","Output/brca/")+"normal_rnd_betas.txt","r").readlines()

print len(f_1[0].strip().split('\t'))
print len(f_2[0].strip().split('\t'))
print len(f_3[0].strip().split('\t'))

f_4 = open(os.getcwd().replace("Scripts","Output/coad/")+"normal_raw_betas.txt","r").readlines()
f_5 = open(os.getcwd().replace("Scripts","Output/coad/")+"normal_dis_betas.txt","r").readlines()
f_6 = open(os.getcwd().replace("Scripts","Output/coad/")+"normal_rnd_betas.txt","r").readlines()

print len(f_4[0].strip().split('\t'))
print len(f_5[0].strip().split('\t'))
print len(f_6[0].strip().split('\t'))
'''



f_3 = open(os.getcwd().replace("Scripts","Output/brca/")+"normal_rnd_betas.txt","r").readlines()

probe_dict = {}
f_line = []
for line in f_3:
	line = line.strip().split('\t')
	if line[0] == "Samples":
		f_line = line[1:]
	else:
		probe_dict[line[0]] = np.array(line[1:],float)

df = pd.DataFrame.from_dict(probe_dict,orient="index")
df.columns = f_line

snp_probes = []
total = 0
for c in df.columns:
	if np.sum(np.isnan(df.loc[:,c]) == True) == 27:
		total += 1
		snp_probes.append(c)
print total
f_out = open(os.getcwd().replace("Scripts","Annotation/jhu_HM27_adf/snp_probes.txt"),"w")
f_out.write('\t'.join(snp_probes))
f_out.close()

'''
total1 = 0
for c in df.columns:
	if np.sum(np.isnan(df.loc[:,c]) == True) > 0:
		total1 += 1
print total1 - total

f_6 = open(os.getcwd().replace("Scripts","Output/coad/")+"normal_rnd_betas.txt","r").readlines()

probe_dict = {}
f_line = []
for line in f_6:
	line = line.strip().split('\t')
	if line[0] == "Samples":
		f_line = line[1:]
	else:
		probe_dict[line[0]] = np.array(line[1:],float)

df = pd.DataFrame.from_dict(probe_dict,orient="index")
df.columns = f_line

total = 0
for c in df.columns:
	if np.sum(np.isnan(df.loc[:,c]) == True) == 37:
		total += 1
print total

total1 = 0
for c in df.columns:
	if np.sum(np.isnan(df.loc[:,c]) == True) > 0:
		total1 += 1
print total1 - total
'''

'''
f_3 = open(os.getcwd().replace("Scripts","Output/brca/")+"normal_rnd_betas.txt","r").readlines()
f_4 = open(os.getcwd().replace("Scripts","Output/brca/")+"normal_rnd_betas_1.txt","r").readlines()
probe_dict_3 = {}
f_line = []
for line in f_3:
	line = line.strip().split('\t')
	if line[0] == "Samples":
		f_line = line[1:]
	else:
		probe_dict_3[line[0]] = np.array(line[1:],float)

total = 0
total2 = 0
for key in probe_dict_3.keys():
	for i in probe_dict_3[key]:
		if np.isnan(i):
			total += 1
		elif i == 0.5:
			total2 += 1
print total
print total+total2

probe_dict_4 = {}
f_line = []
for line in f_4:
	line = line.strip().split('\t')
	if line[0] == "Samples":
		f_line = line[1:]
	else:
		probe_dict_4[line[0]] = np.array(line[1:],float)

total = 0
for key in probe_dict_4.keys():
	for i in probe_dict_4[key]:
		if i == 0.5:
			total += 1
print total
'''
