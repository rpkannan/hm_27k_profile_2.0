import sys, os, time
import numpy as np

###################### MUST BE SET BEFORE USING OTHER METHODS ############################
input_path = "None"
output_path = "None"
annot_path = "None"

# User must set file path in which manifest is found
def set_input_path(path):
	global input_path
	input_path = path
	return

# User must set file path in which output files should be placed
def set_output_path(path):
	global output_path
	output_path = path
	return

def set_annotation_path(path):
	global annot_path
	annot_path = path
	return

#################################### DIAGNOSTICS #########################################
def get_input_path():
	global input_path
	return input_path

def get_output_path():
	global output_path
	return output_path
	
def get_annotation_path():
	global annot_path
	return annot_path	

def get_barcodes(geneexp_dict,file_name):
	global output_path
	barcodes = []
	f_out = open(output_path+file_name,"w")
	f_out.write(",".join(geneexp_dict.keys()))
	return geneexp_dict.keys()

################################# GLOBAL OBJECTS #########################################
gene_symbols = open(str(os.getcwd()).replace("hm_27k_profile_2.0/Scripts","Annotation/unc_Agilent/gene_symbols.txt"),"r").readline().strip().split('\t')

############################### DATA EXTRACTION METHODS ##################################

def load_manifest_1():
	global input_path
	manifest = {}
	man_file = open(input_path+"file_manifest.txt","r").readlines()
	for line in man_file:
		spl = line.strip().split('\t')
		plat = spl[2]
		level = spl[3]
		if (plat == 'AgilentG4502A_07_1' and level == '3'):
			brc = spl[5]
			manifest[brc] = spl[6]
	return manifest

def load_manifest_2():
	global input_path
	manifest = {}
	man_file = open(input_path+"file_manifest.txt","r").readlines()
	for line in man_file:
		spl = line.strip().split('\t')
		plat = spl[2]
		level = spl[3]
		if (plat == 'AgilentG4502A_07_2' and level == '3'):
			brc = spl[5]
			manifest[brc] = spl[6]
	return manifest
	
def extract_ratios_1(fileName):
	global input_path
	ge_dict = {}
	f_in = open(input_path+"Expression-Genes/UNC__AgilentG4502A_07_1/Level_3/"+fileName,"r").readlines()
	for line in f_in:
		line = line.strip().split('\t')
		if line[0] == "Hybridization REF" or line[0] == "Composite Element REF":
			continue
		else:
			ge_dict[line[0]] = line[1]
	return ge_dict

def extract_ratios_2(fileName):
	global input_path
	ge_dict = {}
	f_in = open(input_path+"Expression-Genes/UNC__AgilentG4502A_07_2/Level_3/"+fileName,"r").readlines()
	for line in f_in:
		line = line.strip().split('\t')
		if line[0] == "Hybridization REF" or line[0] == "Composite Element REF":
			continue
		else:
			ge_dict[line[0]] = line[1]
	return ge_dict
		
def build_geneexp_dict(manifest1,manifest2,file_name):
	global output_path
	geneexp_dict = {}
	for brc in sorted(manifest1.keys()):
		ratios_d = extract_ratios_1(manifest1[brc])
		r_ls = []
		for g in sorted(ratios_d.keys()):
			 r_ls.append(ratios_d[g])
		geneexp_dict[brc] = r_ls
	for brc in sorted(manifest2.keys()):
		ratios_d = extract_ratios_2(manifest2[brc])
		r_ls = []
		for g in sorted(ratios_d.keys()):
			 r_ls.append(ratios_d[g])
		geneexp_dict[brc] = r_ls
	
	if file_name == "NA":
		return geneexp_dict
	else:
		f_out = open(output_path+file_name,"w")
		header = ["Samples"] + gene_symbols
		f_out.write('\t'.join(header)+'\n')
		for brc in sorted(geneexp_dict.keys()):
			f_out.write(brc+'\t'+'\t'.join(geneexp_dict[brc])+'\n')
		f_out.close()
		return geneexp_dict
'''
	
file_in1 = open("/Users/rathikannan/Documents/tcga_gene_exp/gbm/Expression-Genes/UNC__AgilentG4502A_07_1/Level_3/US14702406_251584710166_S01_GE2-v5_91_0806.txt_lmean.out.logratio.gene.tcga_level3.data.txt","r").readlines()
for i in range(0,5):
	print file_in1[i].strip().split('\t')

genes1 = []
for line in file_in1:
	line = line.strip().split('\t')
	if line[0] == "Hybridization REF" or line[0] == "Composite Element REF":
		continue
	else:
		genes1.append(line[0])
genes1 = np.array(genes1)
		
file_in2 = open("/Users/rathikannan/Documents/hm_27k_profile_2.0/Input/gbm/Expression-Genes/UNC__AgilentG4502A_07_2/Level_3/US45102955_251780410344_S01_GE2-v5_95_Feb07.txt_lmean.out.logratio.gene.tcga_level3.data.txt","r").readlines()
genes2 = []
for line in file_in2:
	line = line.strip().split('\t')
	if line[0] == "Hybridization REF" or line[0] == "Composite Element REF":
		continue
	else:
		genes2.append(line[0])
genes2 = np.array(genes2)

gene_symbol_to_probe = open("/Users/rathikannan/Documents/Annotation/JHU_HM27_adf/gene_symbol_to_probe.txt","r").readlines()
gsp_dict = {}
for line in gene_symbol_to_probe:
	line = line.strip().split('\t')
	gsp_dict[line[0]] = line[1:]

gene_symbols = np.array(gsp_dict.keys())

#print len(np.intersect1d(gene_symbols,genes1))
#print len(np.intersect1d(gene_symbols,genes2))
#print len(np.intersect1d(genes1,genes2))

file_out = open("/Users/rathikannan/Documents/Annotation/unc_Agilent/gene_symbols.txt","w")
file_out.write('\t'.join(sorted(genes1)))
file_out.close()
'''