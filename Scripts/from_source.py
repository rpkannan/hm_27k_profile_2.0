import sys, os, time
import numpy as np
import from_hm27k
import from_geneexp


main_up_path = os.getcwd().replace("hm_27k_profile_2.0/Scripts","tcga_methyl_27/")
main_out_path = os.getcwd().replace("Scripts","Output/")

cancers = ["brca","coad","gbm","kirc","luad","lusc","ov","ucec"]



##################### Read and Store All Raw Values from HM27K Files #####################
for cancer in cancers:
	
	print "Started: "+str(cancer)
	
	from_hm27k.set_version("B")
	from_hm27k.set_up_path(main_up_path+str(cancer)+"/")
	from_hm27k.set_out_path(main_out_path+str(cancer)+"/")
	
	manifest = from_hm27k.load_manifest()
	probe_dict = from_hm27k.build_probe_dict(manifest,"all_raw_betas.txt",False)
	barcodes = from_hm27k.get_barcodes(probe_dict,"methyl_barcodes.txt")

	print "Finished: "+str(cancer)

'''
################# Read and Store All Raw Values from Gene-Exp Files ######################
main_input_path = os.getcwd().replace("hm_27k_profile_2.0/Scripts","tcga_gene_exp/")
main_output_path = os.getcwd().replace("Scripts","Output/")

cancers = ["gbm"]

for cancer in cancers:
	print "Started: "+str(cancer)
	
	from_geneexp.set_input_path(main_input_path+cancer+"/")
	from_geneexp.set_output_path(main_output_path+cancer+"/")
	
	manifest1 = from_geneexp.load_manifest_1()
	manifest2 = from_geneexp.load_manifest_2()
	geneexp_dict = from_geneexp.build_geneexp_dict(manifest1,manifest2,"all_raw_ratios.txt")
	print "Finished: "+str(cancer)

main_up_path = os.getcwd().replace("hm_27k_profile_2.0/Scripts","tcga_methyl_27/")
main_out_path = os.getcwd().replace("Scripts","Output/")
cancers = ["gbm"]
for cancer in cancers:
	print "Started: "+str(cancer)
	
	from_hm27k.set_version("B")
	from_hm27k.set_up_path(main_up_path+str(cancer)+"/")
	from_hm27k.set_out_path(main_out_path+str(cancer)+"/")
	
	manifest = from_hm27k.load_manifest()
	probe_dict = from_hm27k.build_probe_dict(manifest,"all_raw_betas.txt",False)
	barcodes = from_hm27k.get_barcodes(probe_dict,"methyl_barcodes.txt")
	
	main_up_path = os.getcwd().replace("hm_27k_profile_2.0/Scripts","tcga_gene_exp/")
	from_geneexp.set_input_path(main_up_path+str(cancer)+"/")
	from_geneexp.set_output_path(main_out_path+str(cancer)+"/")
	
	manifest1 = from_geneexp.load_manifest_1()
	manifest2 = from_geneexp.load_manifest_2()
	geneexp_dict = from_geneexp.build_geneexp_dict(manifest1,manifest2,"all_gene_ratios.txt")

	print "Finished: "+str(cancer)
'''