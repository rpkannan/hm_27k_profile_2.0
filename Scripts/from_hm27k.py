import sys, os, time
import numpy as np

###################### MUST BE SET BEFORE USING OTHER METHODS ############################
version = "None"
up_path = "None"
out_path = "None"

# User must define file format version: "A" or "B"
def set_version(v):
	global version
	if v == "A" or v == "B":
		version = v
	else:
		print "Error: please enter valid file format version: A or B"
	return

# User must set file path in which manifest is found
def set_up_path(path):
	global up_path
	up_path = path
	return

# User must set file path in which output files should be placed
def set_out_path(path):
	global out_path
	out_path = path
	return

#################################### DIAGNOSTICS #########################################
def get_up_path():
	global up_path
	return up_path

def get_out_path():
	global out_path
	return out_path	

def get_barcodes(manifest,file_name):
	global out_path
	f_out = open(out_path+file_name,"w")
	f_out.write(",".join(manifest.keys()))
	return manifest.keys()


################################# GLOBAL OBJECTS #########################################
probes_27 = open(str(os.getcwd()).replace("hm_27k_profile_2.0/Scripts","Annotation/jhu_HM27_adf/probes_27.txt"),"r").readline().strip().split('\t')

############################### DATA EXTRACTION METHODS ##################################
'''
Reads file manifest and returns a dictionary
Obeys file structure generated from TCGA data matrix
Output:
Dictionary with Keys: barcode Values: file names
'''
def load_manifest():
	global up_path
	manifest = {}
	man_file = open(up_path+"file_manifest.txt","r").readlines()
	for line in man_file:
		spl = line.strip().split('\t')
		plat = spl[2]
		level = spl[3]
		if (plat == 'HumanMethylation27' and level == '3'):
			brc = spl[5]
			manifest[brc] = spl[6]
	return manifest

'''
Given a float raw beta value, returns a discretized float value.
'''
def discretize(beta):
	if (np.isnan(beta)):
		return np.nan
	elif (beta <= 0.2):
		return 0.0
	elif (beta > 0.6):
		return 2.0
	else:
		return 1.0

'''
For use with version A file type. 
Arguments:
fileName: file name without file path
discretize: True or False
Returns: dictionary of probe: float beta value entries
'''
def extract_betas_A(fileName,discrete):
	global up_path
	if discrete:
		file = open(up_path+"/DNA_Methylation/JHU_USC__HumanMethylation27/Level_3/"+fileName,"r").readlines()
		betas = {}
		for line in file:
			line = line.strip().split('\t')
			if line[2] == "beta value":
				continue
			else:
				if line[2] == "NA":
					betas[line[1]] = np.nan
				else:
					betas[line[1]] = discretize(float(line[2]))
		return betas
	elif not discrete:
		file = open(up_path+"/DNA_Methylation/JHU_USC__HumanMethylation27/Level_3/"+fileName,"r").readlines()
		betas = {}
		for line in file:
			line = line.strip().split('\t')
			if line[2] == "beta value":
				continue
			else:
				if line[2] == "NA":
					betas[line[1]] = np.nan
				else:
					betas[line[1]] = float(line[2])
		return betas
	else:
		sys.exit("Error: please enter valid value for discretize: True or False")

'''	
For use with version B file type.	
Arguments:
fileName: file name without file path
discretize: True or False
Returns: dictionary of probe: beta value entries
'''
def extract_betas_B(fileName,discrete):
	global platform, up_path
	if discrete:
		file = open(up_path+"/DNA_Methylation/JHU_USC__HumanMethylation27/Level_3/"+fileName,"r").readlines()
		betas = {}
		for line in file:
			line = line.strip().split('\t')
			if line[0] == "Hybridization REF" or line[0] == "Composite Element REF":
				continue
			else:
				if line[1] == "NA":
					betas[line[0]] = np.nan
				else:
					betas[line[0]] = discretize(float(line[1]))
		return betas
	elif not discrete:
		file = open(up_path+"/DNA_Methylation/JHU_USC__HumanMethylation27/Level_3/"+fileName,"r").readlines()
		betas = {}
		for line in file:
			line = line.strip().split('\t')
			if line[0] == "Hybridization REF" or line[0] == "Composite Element REF":
				continue
			else:
				if line[1] == "NA":
					betas[line[0]] = np.nan
				else:
					betas[line[0]] = float(line[1])
		return betas
	else:
		sys.exit("Error: please enter valid value for discretize: True or False")

'''
Arguments:
manifest: dictionary
file_name: file name without file path
discrete: True or False; discretize or not
Returns: A dictionary with keys: barcode, values: list of float beta values indexed by probes_27.txt

'''
def build_probe_dict(manifest, file_name, discrete):
	
	global out_path, probes_27, version
	samples = sorted(manifest.keys())
	probe_dict = {}
	
	for s in samples:
		if version == "A":
			betas = extract_betas_A(manifest[s],discrete)
			sorted_betas = []
			for p in probes_27:
				if p in betas:
					sorted_betas.append(betas[p])
				else:
					sorted_betas.append(np.nan)
			probe_dict[s] = np.array(sorted_betas)
		else:
			betas = extract_betas_B(manifest[s],discrete)
			sorted_betas = []
			for p in probes_27:
				if p in betas:
					sorted_betas.append(betas[p])
				else:
					sorted_betas.append(np.nan)
			probe_dict[s] = np.array(sorted_betas)
	if file_name == "NA":
		return probe_dict
	else:
		f_out = open(out_path+file_name,"w")
		f_line = ["Samples"]+probes_27
		f_out.write('\t'.join(f_line)+'\n')
		barcodes = sorted(probe_dict.keys())
		for b in barcodes:
			f_out.write(b+'\t')
			f_out.write('\t'.join(str(a) for a in probe_dict[b]))
			f_out.write('\n')
		f_out.close()
		return probe_dict

