##done on stritch
##GWAS QC for MESA african population
##same QC done for hispanic and european populations

#1. sex check and filter out missing genotypes 
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp --check-sex --missing --out AFA_gtool_conversion_to_ped/afa_all_genexp

#24144 MB RAM detected; reserving 12072 MB for main workspace.
#39295080 variants loaded from .bim file.
#234 people (102 males, 132 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 234 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.964319.
#--missing: Sample missing data report written to
#AFA_gtool_conversion_to_ped/afa_all_genexp.imiss, and variant-based missing
#data report written to AFA_gtool_conversion_to_ped/afa_all_genexp.lmiss.
#39295080 variants and 234 people pass filters and QC.
#Note: No phenotypes present.
#--check-sex: 657490 Xchr and 0 Ychr variant(s) scanned, 2 problems detected.
#Report written to AFA_gtool_conversion_to_ped/afa_all_genexp.sexcheck 

#2. Recalculate individual call rates after removing SNPs with call rates <99%
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp --geno 0.01 --make-bed --out AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01

#24144 MB RAM detected; reserving 12072 MB for main workspace.
#39295080 variants loaded from .bim file.
#234 people (102 males, 132 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 234 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.964319.
#17211163 variants removed due to missing genotype data (--geno).
#22083917 variants and 234 people pass filters and QC.
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01 --missing  --out AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01
#24144 MB RAM detected; reserving 12072 MB for main workspace.
#22083917 variants loaded from .bim file.
#234 people (102 males, 132 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 234 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.997789.
#--missing: Sample missing data report written to
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.imiss, and variant-based
#missing data report written to
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.lmiss.

#3. Calculate HWE statistics to flag SNPs...use later for filtering 
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01 --hardy  --out AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01

#4. LD prune data (rm 1 SNP if r2>0.3 in 50 SNP window) for relationship check and heterozygosity calculation
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01 --indep-pairwise 50 5 0.3 --out AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01

#22083917 variants loaded from .bim file.
#234 people (102 males, 132 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 234 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.997789.
#22083917 variants and 234 people pass filters and QC.
#Note: No phenotypes present.
#Pruned 1519373 variants from chromosome 1, leaving 189843.
#Pruned 2836343 variants from chromosome 2, leaving 461052.
#Pruned 1397155 variants from chromosome 3, leaving 185471.
#Pruned 1343854 variants from chromosome 4, leaving 167797.
#Pruned 1276204 variants from chromosome 5, leaving 175132.
#Pruned 1083543 variants from chromosome 7, leaving 136385.
#Pruned 1101325 variants from chromosome 8, leaving 143005.
#Pruned 823069 variants from chromosome 9, leaving 101197.
#Pruned 944709 variants from chromosome 10, leaving 126335.
#Pruned 953393 variants from chromosome 11, leaving 122248.
#Pruned 915341 variants from chromosome 12, leaving 119717.
#Pruned 687289 variants from chromosome 13, leaving 92150.
#Pruned 621497 variants from chromosome 14, leaving 76789.
#Pruned 549302 variants from chromosome 15, leaving 64487.
#Pruned 583908 variants from chromosome 16, leaving 60811.
#Pruned 492551 variants from chromosome 17, leaving 45690.
#Pruned 532480 variants from chromosome 18, leaving 62454.
#Pruned 353907 variants from chromosome 19, leaving 26816.
#Pruned 419070 variants from chromosome 20, leaving 48628.
#Pruned 241842 variants from chromosome 21, leaving 25015.
#Pruned 224520 variants from chromosome 22, leaving 23933.
#Pruned 632942 variants from chromosome 23, leaving 95345.
#Pruning complete.  19533617 of 22083917 variants removed.
#Marker lists written to
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.prune.in and
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.prune.out

#5. Relationship check with a min PI_HAT of 0.05 yields IBD calculation 
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01 --extract AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.prune.in --genome --min 0.05 --out AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3

#22083917 variants loaded from .bim file.
#234 people (102 males, 132 females) loaded from .fam.
#--extract: 2550300 variants remaining.
#Using up to 8 threads (change this with --threads).
#Before main variant filters, 234 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.995916.
#2550300 variants and 234 people pass filters and QC.
#Note: No phenotypes present.
#Excluding 95345 variants on non-autosomes from IBD calculation.
#IBD calculations complete.  
#Finished writing
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.genome

#6. computes observed and expected autosomal homozygous genotype counts for each sample, and reports method-of-moments F coefficient estimates 
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01 --het --out AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01

#22083917 variants loaded from .bim file.
#234 people (102 males, 132 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 234 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.997789.
#22083917 variants and 234 people pass filters and QC.
#Note: No phenotypes present.
#--het: 5080387 variants scanned, report written to
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.het


#7. Extract the pruned in SNPs and remove related individuals 
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01 --extract AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.prune.in --remove AFA_gtool_conversion_to_ped/afa_rels.txt --make-bed --out AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels
#FID 25700 removed. Related
#22083917 variants loaded from .bim file.
#234 people (102 males, 132 females) loaded from .fam.
#--extract: 2550300 variants remaining.
#--remove: 233 people remaining.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 233 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate in remaining samples is 0.995917.
#2550300 variants and 233 people pass filters and QC.
#Note: No phenotypes present.
#--make-bed to
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.bed +
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.bim +
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.fam

#8.  Check heterozygosity (across all autosomal SNPs) -- look at that distribution across individuals to check for and rm outliers (F: mean +/-3 sd)
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels --het --out AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.h.het

#2550300 variants loaded from .bim file.
#233 people (102 males, 131 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 233 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.995917.
#2550300 variants and 233 people pass filters and QC.
#Note: No phenotypes present.
#--het: 2452252 variants scanned, report written to
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.h.het

#9. remove outliers from het check
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels --remove AFA_gtool_conversion_to_ped/afa.het_outlier.txt --make-bed --out AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout

#outlier hets removed. 3. 24906, 25490, 25038

#2550300 variants loaded from .bim file.
#233 people (102 males, 131 females) loaded from .fam.
#--remove: 230 people remaining.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 230 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate in remaining samples is 0.995948.
#2550300 variants and 230 people pass filters and QC.
#Note: No phenotypes present.
#--make-bed to
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout.bed +
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout.bim +
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout.fam

#10. Merge population with hapmap data **can be tricky...remove problem SNPs
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout --bmerge pedmap/hapmap_and_merge/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig --make-bed --out AFA_gtool_conversion_to_ped/afa_hap_merge

#2550300 markers loaded from
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout.bim.
#1499880 markers to be merged from
#pedmap/hapmap_and_merge/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.bim.
#Of these, 1045437 are new, while 454443 are present in the base dataset.
#176 more multiple-position warnings: see log file.
#Error: 114396 variants with 3+ alleles present.
#* If you believe this is due to strand inconsistency, try --flip with
#  AFA_gtool_conversion_to_ped/afa_hap_merge-merge.missnp.
#  (Warning: if this seems to work, strand errors involving SNPs with A/T or C/G
#  alleles probably remain in your data.  If LD between nearby SNPs is high,
#  --flip-scan should detect them.)
#* If you are dealing with genuine multiallelic variants, we recommend exporting
#  that subset of the data to VCF (via e.g. '--recode vcf'), merging with
#  another tool/script, and then importing the result; PLINK is not yet suited
#  to handling them.

#11. **best way to merge with hapmap 
./plink --bfile pedmap/hapmap_and_merge/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig --exclude AFA_gtool_conversion_to_ped/afa_hap_merge-merge.missnp --make-bed --out AFA_gtool_conversion_to_ped/hap1

./plink --bfile AFA_gtool_conversion_to_ped/hap1 --exclude AFA_gtool_conversion_to_ped/afa_hap_merge.missnip --make-bed --out AFA_gtool_conversion_to_ped/hap_for_merge

#12. make bed/bim/fam of population merged with hapmap pops  
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout --bmerge AFA_gtool_conversion_to_ped/hap_for_merge --make-bed --out AFA_gtool_conversion_to_ped/afa_hap_merged

#230 people loaded from
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout.fam.
#391 people to be merged from AFA_gtool_conversion_to_ped/hap_for_merge.fam.
#Of these, 391 are new, while 0 are present in the base dataset.
#2550300 markers loaded from
#AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout.bim.
#1376092 markers to be merged from
#AFA_gtool_conversion_to_ped/hap_for_merge.bim.
#Of these, 1045437 are new, while 330655 are present in the base dataset.
#Warning: Variants 'rs7555212' and 'AFFX-SNP_6564577__rs7555212' have the same position.
#Warning: Variants 'rs10915322' and 'AFFX-SNP_6869948__rs10915322' have the same position.
#Warning: Variants 'rs9439505' and 'AFFX-SNP_4600__rs9439505' have the same position.
#2865 more same-position warnings: see log file.
#Performing single-pass merge (621 people, 3595737 variants).
#Merged fileset written to AFA_gtool_conversion_to_ped/afa_hap_merged-merge.bed
#+ AFA_gtool_conversion_to_ped/afa_hap_merged-merge.bim +
#AFA_gtool_conversion_to_ped/afa_hap_merged-merge.fam .
#3595737 variants loaded from .bim file.
#621 people (298 males, 323 females) loaded from .fam.
#391 phenotype values loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 621 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Warning: 1676 het. haploid genotypes present (see
#AFA_gtool_conversion_to_ped/afa_hap_merged.hh ); many commands treat these as missing.
#Total genotyping rate is 0.466668.
#3595737 variants and 621 people pass filters and QC.
#Among remaining phenotypes, 0 are cases and 391 are controls.  (230 phenotypes are missing.)
#--make-bed to AFA_gtool_conversion_to_ped/afa_hap_merged.bed +
#AFA_gtool_conversion_to_ped/afa_hap_merged.bim +
#AFA_gtool_conversion_to_ped/afa_hap_merged.fam

#13. remove sex and other chromosomes that arent autosomal and LD prune and filter by MAF > 0.05
./plink --bfile AFA_gtool_conversion_to_ped/afa_hap_merged --geno 0.01 --maf 0.05 --autosome --indep-pairwise 50 5 0.3 --out AFA_gtool_conversion_to_ped/afa_hap_mergeLDpn

#3462931 variants loaded from .bim file.
#621 people (298 males, 323 females) loaded from .fam.
#391 phenotype values loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 621 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.467119.
#3259294 variants removed due to missing genotype data (--geno).
#3872 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#199765 variants and 621 people pass filters and QC.
#Among remaining phenotypes, 0 are cases and 391 are controls.  (230 phenotypes are missing.)
#Pruned 10058 variants from chromosome 1, leaving 11799.
#Pruned 2431 variants from chromosome 2, leaving 4810.
#Pruned 1766 variants from chromosome 3, leaving 3774.
#Pruned 1478 variants from chromosome 4, leaving 3283.
#Pruned 1525 variants from chromosome 5, leaving 3502.
#Pruned 6880 variants from chromosome 7, leaving 8061.
#Pruned 7179 variants from chromosome 8, leaving 7790.
#Pruned 5699 variants from chromosome 9, leaving 6998.
#Pruned 6636 variants from chromosome 10, leaving 7954.
#Pruned 6589 variants from chromosome 11, leaving 7269.
#Pruned 5956 variants from chromosome 12, leaving 7389.
#Pruned 4683 variants from chromosome 13, leaving 5726.
#Pruned 4050 variants from chromosome 14, leaving 4927.
#Pruned 3640 variants from chromosome 15, leaving 4956.
#Pruned 3447 variants from chromosome 16, leaving 5043.
#Pruned 2688 variants from chromosome 17, leaving 4319.
#Pruned 3478 variants from chromosome 18, leaving 4914.
#Pruned 1476 variants from chromosome 19, leaving 2776.
#Pruned 3008 variants from chromosome 20, leaving 4056.
#Pruned 1604 variants from chromosome 21, leaving 2407.
#Pruned 1406 variants from chromosome 22, leaving 2335.
#Pruning complete.  85677 of 199765 variants removed.
#Marker lists written to AFA_gtool_conversion_to_ped/afa_hap_mergeLDpn.prune.in
#and AFA_gtool_conversion_to_ped/afa_hap_mergeLDpn.prune.out

#14. extract pruned in SNPs from previous step Ready for smartpca

./plink --bfile AFA_gtool_conversion_to_ped/afa_hap_merged --geno 0.01 --maf 0.05 --autosome --extract AFA_gtool_conversion_to_ped/afa_hap_mergeLDpn.prune.in --recode --out AFA_gtool_conversion_to_ped/afa_hap_merged_pruned

#3462931 variants loaded from .bim file.
#621 people (298 males, 323 females) loaded from .fam.
#391 phenotype values loaded from .fam.
#--extract: 114088 variants remaining.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 621 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.998386.
#0 variants removed due to missing genotype data (--geno).
#0 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#114088 variants and 621 people pass filters and QC.
#Among remaining phenotypes, 0 are cases and 391 are controls.  (230 phenotypes are missing.)
#--recode ped to AFA_gtool_conversion_to_ped/afa_hap_merged_pruned.ped +
#AFA_gtool_conversion_to_ped/afa_hap_merged_pruned.map ... done

#15. make fam file in correct format for smartpca 

awk '{print $1,$2,$3,$4,$5,1}' AFA_gtool_conversion_to_ped/afa_hap_merged.fam > AFA_gtool_conversion_to_ped/afa_hap_merged_pruned.fam

#16. command to run smartpca by EIGENSTRAT

./smartpca -p afamerge.par

#parameter file: ../afamerge.par
### THE INPUT PARAMETERS
##PARAMETER NAME: VALUE
#genotypename: /home/lmogil/data/convert_to_plink/AFA_gtool_conversion_to_ped/afa_hap_merged_pruned.ped
#snpname: /home/lmogil/data/convert_to_plink/AFA_gtool_conversion_to_ped/afa_hap_merged_pruned.map
#indivname: /home/lmogil/data/convert_to_plink/AFA_gtool_conversion_to_ped/afa_hap_merged_pruned.fam
#evecoutname: /home/lmogil/data/convert_to_plink/AFA_gtool_conversion_to_ped/afamerged5.pca.evec
#evaloutname: /home/lmogil/data/convert_to_plink/AFA_gtool_conversion_to_ped/afamerged5.eval
#numoutevec: 10
#numoutlieriter: 0
#numoutlierevec: 2
#outliersigmathresh: 6
## smartpca version: 13050
#norm used

#genetic distance set from physical distance
#number of samples used: 621 number of snps used: 114088
#Using 8 threads, and partial sum lookup algorithm.
#total number of snps killed in pass: 0  used: 114088

#17. recode for afa alone (no merge with hapmap) 
./plink --bfile AFA_gtool_conversion_to_ped/afa_all_genexp.geno0.01.LD0.3.rmrels.rmout --recode --out AFA_gtool_conversion_to_ped/afa_pca

#2550300 variants loaded from .bim file.
#230 people (101 males, 129 females) loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 230 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.995948.
#2550300 variants and 230 people pass filters and QC.
#Note: No phenotypes present.
#--recode ped to AFA_gtool_conversion_to_ped/afa_pca.ped +
#AFA_gtool_conversion_to_ped/afa_pca.map


