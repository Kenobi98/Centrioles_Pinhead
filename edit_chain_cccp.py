###########################################################
####$$Edit the pdb files obtained from CCCP$$####
#@@Author: Kartik Majila

import sys

file = sys.argv[1]       #Input the file name 

#Read the pdb file
with open (file, "r") as f:
    pdb = f.readlines()

#This script will overwrite the existing PDB file
w = open(file, "w")
for i in range(len(pdb)):
    if pdb[i].startswith("ATOM"):

        #In CCCP files chain name is present at 72 position in each line
        chain = pdb[i][72] 

        #Add the chain name at position 22 as in PDB files usually
        a = pdb[i][:21] + chain + pdb[i][22:]
        w.writelines(a)
    else:
        w.writelines(pdb[i])
w.close()




