import pandas as pd
import os
from chimera import runCommand as rc
from chimera import replyobj

data = pd.read_excel('/Users/oliveradams/Desktop/BSc Project/relevant_PDBs.xls')

for i in range(len(data)-1):
    replyobj.status("processing" + data['PDB_ID'][i])
    rc("open " + data['PDB_ID'][i])
    rc("addh")
    rc("sel :" + data['RES. No.'][i]+"."+data['RESCHAIN'][i]+"|:"+(data['RES. No.'][i]-1)+"."+data['RESCHAIN'][i]+"@C|:"+(data['RES. No.'][i]+1)+"."+data['RESCHAIN'][i]+"@N")
    rc("write selected format pdb 0 arg_" + data[PDB_ID])