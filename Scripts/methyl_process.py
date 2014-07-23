import sys, os, time
import numpy as np

output_path = "None"
snp_probes = open(os.getcwd().replace("Scripts","Annotation/jhu_HM27_adf/snp_probes.txt"),"r").readline().strip().split('\t')

def set_path(path):
	global output_path
	output_path = path
	
def load_barcodes(filename):
	global output_path
	barcodes = open(output_path+filename,"r").readline().strip().split(",")
	return barcodes

def check_brc(brc):
	brc_s = brc.split('-')
	if brc_s[3][0:2] in ['01','02','03','04','05','06','07','08','09']:
		return "tumor"
	elif brc_s[3][0:2] in ['10','11','12','13','14','15','16','17','18','19']:
		return "normal"
	else:
		return "other"

'''
Given a float raw beta value, returns a discretized float value.
'''
def discretize_list(beta_ls):
	new_ls = []
	for beta in beta_ls:
		if (beta == "nan"):
			new_ls.append(np.nan)
		elif (float(beta) <= 0.2):
			new_ls.append(0.0)
		elif (float(beta) > 0.6):
			new_ls.append(2.0)
		else:
			new_ls.append(1.0)
	return np.array(new_ls)

def replace_NA(beta_ls,discrete):
	new_ls = []
	if discrete:
		for beta in beta_ls:
			beta = float(beta)
			if np.isnan(beta):
				new_ls.append(1.0)
			else:
				new_ls.append(beta)
	else:
		for beta in beta_ls:
			beta = float(beta)
			if np.isnan(beta):
				new_ls.append(0.5)
			else:
				new_ls.append(beta)
	return new_ls
	

def iround(num):
	if num == "nan":
		return np.nan
	else:
		return round(float(num),2)

def modify_ls(beta_ls,discretize,round,replace_na):
	if discretize and replace_na:
		new_ls = 

'''		
Input:
filename: "all_raw_betas.txt"
option: "tumor", "normal", "all"
discretize: True/False
round: True/False
replace_na: 0- True/False
file_out: if "NA", output not written to file
Output:
probe_dict: keys: barcodes values: discretized numpy list
'''
def load_discrete_betas(file_in,option,discretize,round,replace_na,file_out):
	global output_path
	probe_dict = {}
	f_in = open(output_path+file_in,"r").readlines()
	if option == "tumor":
		for line in f_in:
			line = line.strip().split('\t')
			if line[0] == "Samples":
				continue
			else:
				if check_brc(line[0]) == "tumor":
					if discretize:
						if replace_na:
							probe_dict[line[0]] = replace_NA(discretize_list(line[1:]),True)
						else:
							probe_dict[line[0]] = discretize_list(line[1:])
					else:
						if round:
							if replace_na:
								temp = np.array(line[1:])
								probe_dict[line[0]] = replace_NA([iround(i) for i in temp],False)
							else:
								temp = np.array(line[1:])
								probe_dict[line[0]] = [iround(i) for i in temp]
						else:
							if replace_na:
								probe_dict[line[0]] = replace_NA(np.array(line[1:]),False)
							else:
								probe_dict[line[0]] = np.array(line[1:])
	elif option == "normal":
		for line in f_in:
			line = line.strip().split('\t')
			if line[0] == "Samples":
				continue
			else:
				if check_brc(line[0]) == "normal":
					if discretize:
						if replace_na:
							probe_dict[line[0]] = replace_NA(discretize_list(line[1:]),True)
						else:
							probe_dict[line[0]] = discretize_list(line[1:])
					else:
						if round:
							if replace_na:
								temp = np.array(line[1:])
								probe_dict[line[0]] = replace_NA([iround(i) for i in temp],False)
							else:
								temp = np.array(line[1:])
								probe_dict[line[0]] = [iround(i) for i in temp]
						else:
							if replace_na:
								probe_dict[line[0]] = replace_NA(np.array(line[1:]),False)
							else:
								probe_dict[line[0]] = np.array(line[1:])				
	elif option == "all":
		for line in f_in:
			line = line.strip().split('\t')
			if line[0] == "Samples":
				continue
			else:
				if discretize:
					if replace_na:
						probe_dict[line[0]] = replace_NA(discretize_list(line[1:]),True)
					else:
						probe_dict[line[0]] = discretize_list(line[1:])
				else:
					if round:
						if replace_na:
							temp = np.array(line[1:])
							probe_dict[line[0]] = replace_NA([iround(i) for i in temp],False)
						else:
							temp = np.array(line[1:])
							probe_dict[line[0]] = [iround(i) for i in temp]
					else:
						if replace_na:
							probe_dict[line[0]] = replace_NA(np.array(line[1:]),False)
						else:
							probe_dict[line[0]] = np.array(line[1:])
	else:
		print "Error: please enter valid option"
		return probe_dict
	if file_out == "NA":
		return probe_dict
	else:
		f_out = open(output_path+file_out,"w")
		f_out.write(f_in[0])
		for pr in sorted(probe_dict.keys()):
			f_out.write(pr+'\t'+'\t'.join(str(a) for a in probe_dict[pr])+'\n')
		f_out.close()
		return probe_dict

def remove_snp(beta_ls):
	global snp_probes
	
	


###################################### MAIN METHOD #######################################	
main_path = os.getcwd().replace("Scripts","Output/")
cancers = ["brca","coad","gbm","kirc","luad","lusc","ov","ucec"]
cancers = ["brca","coad"]
for cancer in cancers:
	set_path(main_path+cancer+"/")
	#load_discrete_betas("all_raw_betas.txt","normal",True,True,False,"normal_dis_betas.txt")
	#load_discrete_betas("all_raw_betas.txt","normal",False,False,False,"normal_raw_betas.txt")
	#load_discrete_betas("all_raw_betas.txt","normal",False,True,False,"normal_rnd_betas.txt")
	load_discrete_betas("all_raw_betas.txt","normal",True,True,True,"normal_dis_betas_1.txt")
	load_discrete_betas("all_raw_betas.txt","normal",False,False,True,"normal_raw_betas_1.txt")
	load_discrete_betas("all_raw_betas.txt","normal",False,True,True,"normal_rnd_betas_1.txt")
	print cancer
print "C'est fini!"



	
