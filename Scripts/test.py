import sys, os, time

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
