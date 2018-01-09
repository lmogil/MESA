from __future__ import division
import gzip
import re
import sys
from collections import defaultdict

#c = "22"
c = sys.argv[1]

afa_dosage ="/home/lauren/mesa_dosages/imp2_AFA_chr"+c+"_ref1kg_dosage.txt.gz"
his_dosage="/home/lauren/mesa_dosages/imp2_HIS_chr"+c+"_ref1kg_dosage.txt.gz"



refdict = defaultdict(set)
for i in gzip.open(afa_dosage):
	afad = i.strip().split() #split by white space
	(chra, rsida, posa, a1a, a2a,mafa) = afad[0:6] #assign first 5 indices
	adose = afad[6:]
	adose = str(afad[6:])
	adose = adose.replace("[", "")
	adose = adose.replace("]", "")
	refdict[rsida].add(adose)
	
ofile=open("/home/lauren/mesa_dosages/imp2_AFHI_chr"+c+"_ref1kg_dosage.txt","wb")
for j in gzip.open(his_dosage):
	caud = j.strip().split() #split by white space
	(chrc, rsidc, posc, a1c, a2c,mafc) = caud[0:6] #assign first 5 indices
	cdose =caud[6:]
	cdose =str(caud[6:])
	cdose = cdose.replace("[", "")
	cdose = cdose.replace("]", "")
	cdose = cdose.replace("'", "")
	if caud[1] in refdict:
		ad = refdict[caud[1]]
		ad=map(str, ad)
		ad=str(ad)
		ad=ad.replace("[", "")
		ad=ad.replace("]", "")
		ad=ad.replace('"', "")
		ad=ad.replace("'", "")
        ad=ad.replace(",", " ")
		out= chrc+ ' '+rsidc+ ' '+posc+ ' '+a1c+ ' '+a2c+ ' '+mafc+' '+ad+ ' '+cdose+ '\n'
        o.append(out)

ofile.write(out)

	
ofile.close()	
