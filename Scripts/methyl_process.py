import sys, os, time
import numpy as np

output_path = "None"
snp_probes = open(os.getcwd().replace("Scripts","Annotation/jhu_HM27_adf/snp_probes.txt"),"r").readline().strip().split('\t')
probes_27 = open(os.getcwd().replace("Scripts","Annotation/jhu_HM27_adf/probes_27.txt"),"r").readline().strip().split('\t')
snp_mask = np.in1d(probes_27,snp_probes) #True: probe is snp False: probe is not snp


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

def remove_snp(beta_ls):
	global snp_mask
	new_ls = []
	for i,j in zip(beta_ls,snp_mask):
		if not j: new_ls.append(i)
	return new_ls
	
def modify_ls(beta_ls,beta_option,replace_na):
	new_ls = []
	if beta_option == "discrete":
		for beta in beta_ls:
			if beta == "nan":
				if replace_na: new_ls.append(1.0)
				else: new_ls.append(np.nan)
			elif float(beta) <= 0.2:
				new_ls.append(0.0)
			elif float(beta) > 0.6:
				new_ls.append(2.0)
			else:
				new_ls.append(1.0)
	elif beta_option == "round":
		for beta in beta_ls:
			if beta == "nan":
				if replace_na: new_ls.append(0.5)
				else: new_ls.append(np.nan)
			else:
				new_ls.append(round(float(beta),2))
	elif beta_option == "original":
		for beta in beta_ls:
			if beta == "nan":
				if replace_na: new_ls.append(0.5)
				else: new_ls.append(np.nan)
			else:
				new_ls.append(float(beta))
	else:
		print "Error: please enter valid beta_option: discrete, round, original"
		return
	return new_ls
			

'''
Input
file_in: 
tissue_option: "tumor","normal","all"
beta_option: "discrete","round","original"
replace_na: True/False
rm_snp: True/False
file_out:
'''
def load_betas(file_in,tissue_option,beta_option,replace_na,rm_snp,file_out):
	# check inputs
	if tissue_option not in ["tumor","normal","all"]:
		print "Error: please enter valid tissue_option: tumor, normal, all"
		return
	if beta_option not in ["discrete","round","original"]:
		print "Error: please enter valid beta_option: discrete, round, original"
		return
	if not isinstance(replace_na,bool):
		print "Error: please enter valid replace_na: True, False"
		return
	if not isinstance(rm_snp,bool):
		print "Error: please enter valid rm_snp: True, False"
		return
	
	# build probe_dict
	global output_path,snp_mask
	probe_dict = {}
	f_in = open(output_path+file_in,"r").readlines()
	if rm_snp:
		for line in f_in[1:]:
			line = line.strip().split('\t')
			if check_brc(line[0]) == tissue_option: probe_dict[line[0]] = modify_ls(remove_snp(line[1:]),beta_option,replace_na)
	else:
		for line in f_in[1:]:
			line = line.strip().split('\t')
			if check_brc(line[0]) == tissue_option: probe_dict[line[0]] = modify_ls(line[1:],beta_option,replace_na)
	
	# write to file
	if file_out == "NA":
		return probe_dict
	else:
		f_out = open(output_path+file_out,"w")
		if rm_snp:
			header = []
			first_ln = f_in[0].strip().split('\t')
			for i,j in zip(first_ln[1:],snp_mask):
				if not j: header.append(i)
			f_out.write("Samples"+'\t'+'\t'.join(header))
		else:
			f_out.write(f_in[0])
		for pr in sorted(probe_dict.keys()):
			f_out.write(pr+'\t'+'\t'.join(str(a) for a in probe_dict[pr])+'\n')
		f_out.close()
		return probe_dict

###################################### MAIN METHOD #######################################	
main_path = os.getcwd().replace("Scripts","Output/")
cancers = ["brca","coad","gbm","kirc","luad","lusc","ov","ucec"]
cancers = ["brca","coad"]
for cancer in cancers:
	set_path(main_path+cancer+"/")
	load_betas("all_raw_betas.txt","normal","round",True,True,"normal_rnd_betas.txt")
	load_betas("all_raw_betas.txt","normal","discrete",True,True,"normal_dis_betas_tes.txt")
	load_betas("all_raw_betas.txt","tumor","round",True,True,"tumor_rnd_betas.txt")
	load_betas("all_raw_betas.txt","tumor","discrete",True,True,"tumor_dis_betas_tes.txt")
	load_betas("all_raw_betas.txt","all","round",True,True,"all_rnd_betas.txt")
	load_betas("all_raw_betas.txt","all","discrete",True,True,"all_dis_betas_tes.txt")


	print cancer
print "C'est fini!"



	
