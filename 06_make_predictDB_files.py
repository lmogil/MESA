import gzip
import sys
from collections import defaultdict

pop=sys.argv[1]
chr=sys.argv[2]

#pop='CAU'
#chr='22'


dosage = "/home/lauren/mesa_predixcan_dosages/imp2_"+pop+"_chr"+chr+"_ref1kg_dosage.txt.gz"
meqtl="/home/lauren/meqtl_genotypes/imp2_"+pop+"_chr"+chr+"_ref1kg_dosage.txt"

dosedict = defaultdict(set)
for line in open(meqtl):
        if line.startswith('rsid') == False:
                arr = line.strip().split()
                rsid = arr[0]
                d=arr[1:]
                d = str(d[0:])
                d = d.replace("[", "")
                d = d.replace("]", "")
                dosedict[rsid].add(d)

anot_file=open("/home/lauren/files_for_PredDB/"+pop+"_"+chr+"_annot.txt","wb")
snp_file=open("/home/lauren/files_for_PredDB/"+pop+"_"+chr+"_snp.txt","wb")

anot_file.write('Chr'+'\t'+'pos'+'\t'+'var_id'+'\t'+'ref_a1'+'\t'+'alt_a2'+'\t'+'snp_id'+'\t'+'RSID_dbSNP137'+'\t'+'num_alt'+'\n')

        
for line in gzip.open(dosage):
        arr=line.strip().split()
        (c, rs, pos, a1, a2, maf)=arr[0:6]
        varid= chr+'_'+pos+'_'+a1+'_'+a2+'_'+'b37'
        snpid='snp'+'_'+chr+'_'+pos
        if rs in dosedict:
                dose=dosedict[rs]
                dose=map(str, dose)
                dose=str(dose)
                dose = dose.replace("[", "")
                dose = dose.replace("]", "")
                dose = dose.replace("'", "")
                dose = dose.replace('"', "")
                dose = dose.replace(",", "\t")
                snp= varid+'\t'+dose+'\n'
                annot= c+'\t'+pos+'\t'+varid+'\t'+a1+'\t'+a2+'\t'+snpid+'\t'+rs+'\t'+'1'+'\n'
                anot_file.write(annot)
                snp_file.write(snp)
        
        
        
anot_file.close()
snp_file.close()
