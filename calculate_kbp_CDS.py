 #!/usr/bin/python3

"""

Code to calculate kbp CDS per species from a folder of multiple protein alignments.

"""

import os
from Bio import SeqIO
import csv

#create empty dictionary to input taxids as keys and total kbp CDS as values
CDS = {}
#iterate through all alignment files in the working directory
for filename in os.listdir("os.curdir"):
    if filename.endswith(".faa"):
        for record in SeqIO.parse(filename, "fasta"):
            x = record.seq
            y = str(x)
            z = y.replace("-","")
            kbp = len(z)
            a = record.id
            b = str(a)
            taxid = b.split('|')[0]
            if taxid in CDS.keys():
                c = CDS.get(taxid)
                total = c + kbp
            else:
                total = kbp
            CDS.update({taxid : total})
    else:
        continue
print ("Dict taxid-kbp CDS are : ") 
for i in CDS : 
            print(i, CDS[i]) 

with open('kbp_CDS.csv', 'w') as f:  
    w = csv.DictWriter(f, CDS.keys())
    w.writeheader()
    w.writerow(CDS)
