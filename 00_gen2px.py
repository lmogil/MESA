#!/usr/bin/python
'''This python script takes the MESA IMPUTE2 prob (gen) files as input, 
filters by info score (>0.8) and maf (>0.01),
keeps individuals with expression data,
and makes output files for each autosome for PrediXcan:
imp2_<pop>_chr<chr>_ref1kg_dosage.txt.gz
samples_gen3_<pop>.txt
dose allele is Allele2(B), see https://github.com/hakyimlab/PrediXcan/
'''

#each chr takes 1-12 hrs depending on cluster use

import gzip
import sys

pop = sys.argv[1]
chr = sys.argv[2]
date = sys.argv[3]

zfillchr = chr.zfill(2) #mesa files use chr01, etc., but PrediXcan doesn't

chrfile = "/data/hwheeler1/MesaSHARe/Imputation/Generation_3/" + pop + " imputation v.3." + \
pop.lower() + ".1 (1000G)/MESA_SHARe_imputation_" + pop + "_1kg_" + date + "/imputation_data/imp2_" + \
pop + "_chr" + zfillchr + "_ref1kg.gz"

#make info score dict, only include SNPs with MAF>0.01
#keys -> rsid, values -> info
infodict = {}
infofile = "/data/hwheeler1/MesaSHARe/Imputation/Generation_3/" + pop + " imputation v.3." + \
pop.lower() + ".1 (1000G)/MESA_SHARe_imputation_" + pop + "_1kg_" + date + "/SNP_info/imp2_" +  pop + "_chr" + zfillchr + "_ref1kg_info"

for line in open(infofile):
    arr = line.strip().split()
    (blank, rsid, pos, exp_freq_a1, info) = arr[0:5]
    if blank == "snp_id":
        next
    elif float(exp_freq_a1) > 0.01 and (1 - float(exp_freq_a1)) > 0.01 and float(info) > 0.8: 
        infodict[rsid] = info

# get a list of ids that have exp data
expidfile = "/data/hwheeler1/mesaQC/MESA_samples_with_exp.txt"
expidlist = open(expidfile).read().split()

# get list of impute2 sample ids, need FID IID format
samplefile = "/data/hwheeler1/MesaSHARe/Imputation/Generation_3/" + pop + " imputation v.3." + \
pop.lower() + ".1 (1000G)/MESA_SHARe_imputation_" + pop + "_1kg_" + date + "/imputation_sample_file_gen3_" + pop + ".txt"
samplelist = open(samplefile).read().split()
outsamples = open("/data/hwheeler1/mesaQC/predixcan_dosages/samples_gen3_" + pop + ".txt","w")
newlist = [sam + ' ' + sam for sam in samplelist if sam in expidlist]
outsamples.write('\n'.join(newlist) + '\n')
outsamples.close()

# get dosage file data
#SNP ID, RS #, base-pair position, allele coded A, and allele coded B
outdosage = gzip.open("/data/hwheeler1/mesaQC/predixcan_dosages/imp2_" + pop + \
"_chr" + chr + "_ref1kg_dosage.txt.gz","wb")
for line in gzip.open(chrfile):
    arr = line.strip().split()
    (blank, rsid, pos, alleleA, alleleB) = arr[0:5]
    if rsid in infodict:
        info = infodict[rsid]
        gtProbs = arr[5:]
        assert len(gtProbs) % 3 == 0, "gtProbs not divisible by 3"
        gtProbs = map(float,gtProbs)
        #make alleleB dosage list
        dosagerow = []
        for i in range(0,len(gtProbs)-2,3):
            dos = gtProbs[i+1] + gtProbs[i+2]*2
	    dosagerow.append(round(dos,4))
        freqB = round(sum(dosagerow)/(len(dosagerow)*2),4) #calc  alleleB freq (I found that ALT is not always the minor allele)
        #keep dosages with exp data, dosagerow and samplelist are same length/order
        dosagerow = [dosagerow[i] for i in range(len(dosagerow)) if samplelist[i] in expidlist]
        dosages = ' '.join(map(str,dosagerow))
        output = 'chr' + chr + ' ' + rsid + ' ' + pos + ' ' + alleleA + ' ' + alleleB + ' ' + str(freqB) + ' ' + dosages + '\n'
        outdosage.write(output)

outdosage.close()

