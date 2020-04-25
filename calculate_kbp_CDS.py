 #!/usr/bin/python3

"""

Code to calculate bp CDS per species from a folder of multiple protein alignments.

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
            bp = len(z)
            a = record.id
            b = str(a)
            taxid = b.split('|')[0]
            if taxid in CDS.keys():
                c = CDS.get(taxid)
                total = c + bp
            else:
                total = bp
            CDS.update({taxid : total})
    else:
        continue
print ("Dict taxid-bp CDS are : ") 
for i in CDS : 
            print(i, CDS[i]) 

with open('bp_CDS.csv', 'w') as f:  
    w = csv.DictWriter(f, CDS.keys())
    w.writeheader()
    w.writerow(CDS)
