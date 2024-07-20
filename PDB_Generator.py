def endcap_cat(inputfile, residue_code):
    "Selects the relevant section of the pdb for the cation residue and applies the end caps. The output is a string of lines to be put into the modified pdb."
    with open(inputfile, 'r+') as file:
        file = file.readlines()
        num_of_lines = 0
        
        if residue_code == 'ARG':
            for line in file:
                num_of_lines = num_of_lines + 1
                #select the ARG section of the PDB + end caps to be editted
                if line[17:20]==residue_code:
                    #account for other potential ARG/LYS/HIS residues
                    if int(file[num_of_lines-2][24:27])==(int(file[num_of_lines-1][24:27])-1):
                        #Number of atoms in ARG residue fixed therefore lines selected using this info
                        residue = (file[num_of_lines-4:num_of_lines+24])
                        residue_number = int(file[num_of_lines][24:27])
        
        if residue_code == 'LYS':
            for line in file:
                num_of_lines = num_of_lines + 1
                #select the LYS section of the PDB + end caps to be editted
                if line[17:20]==residue_code:
                    #account for other potential ARG/LYS/HIS residues
                    if int(file[num_of_lines-2][24:27])==(int(file[num_of_lines-1][24:27])-1):
                        #Number of atoms in LYS residue fixed therefore lines selected using this info
                        residue = (file[num_of_lines-4:num_of_lines+22])
                        residue_number = int(file[num_of_lines][24:27])
        
        if residue_code == 'HIS':
            for line in file:
                num_of_lines = num_of_lines + 1
                #select the HIS section of the PDB + end caps to be editted
                if line[17:20]==residue_code:
                    #account for other potential ARG/LYS/HIS residues
                    if int(file[num_of_lines-2][24:27])==(int(file[num_of_lines-1][24:27])-1):
                        #Number of atoms in HIS residue fixed therefore lines selected using this info
                        residue = (file[num_of_lines-4:num_of_lines+18])
                        residue_number = int(file[num_of_lines][24:27])
        
        capped_residue = '' 
        num_of_lines_2 = 0
        for line in residue:
            num_of_lines_2 = num_of_lines_2 + 1
            #replace C on residue number - 1 with H plus remove other atoms in this residue
            if int(line[24:27]) == residue_number - 1 and line[77]=='C' and (residue[num_of_lines_2][13:15]=='H1' or int(residue[num_of_lines_2][24:27])==residue_number):
            #num_of_lines_2 term needed to account for the neighbouring residue potentially being captured by Chimera's zr function
                line=list(line)
                line[77]='H'
                line[13]='H'
                line=''.join(line)
                capped_residue = capped_residue + line
            if int(line[24:27]) == residue_number:
                capped_residue = capped_residue + line
            #replace N on residue number + 1 with H plus remove other atoms in this residue
            if int(line[24:27]) == residue_number + 1 and line[77]=='N' and int(residue[num_of_lines_2-2][24:27])==residue_number:
                line=list(line)
                line[77]='H'
                line[13]='H'
                line=''.join(line)
                capped_residue = capped_residue + line
   
    return capped_residue

def endcap_pi(inputfile):
    "Selects the relevant section of the pdb for the pi residue and applies the end cap. The output is a list of lines to be put into the modified pdb."
    with open(inputfile, 'r+') as file:
        file = file.readlines()
        residue = ''
        for line in file:
            #Need to specify atom at LINK lines have residue code at the same index
            if line[:4]=='ATOM':
                #Select the TRP/TYR/PHE residue
                if line[17:20]=='TRP' or line[17:20]=='TYR' or line[17:20]=='PHE':
                    residue = residue + line
   
       #Next section addressess issue of a few atoms from an unwanted TRP/TYR/PHE residue being present in pdb
        #NB this is a result of a TRP/TYR/PHE residue neighbouring the relevant cation residue
        residue = residue.splitlines()
        num_TRP = 0; res_TRP = ''
        num_TYR = 0; res_TYR = ''
        num_PHE = 0; res_PHE = ''
        for line in residue:
            if line[17:20]=='TRP':
                num_TRP = num_TRP + 1
                res_TRP = res_TRP + line + '\n'
            if line[17:20]=='TYR':
                num_TYR = num_TYR + 1
                res_TYR = res_TYR + line + '\n'
            if line[17:20]=='PHE':
                num_PHE = num_PHE + 1
                res_PHE = res_PHE + line + '\n'
        #By comparing the number of times each type of residue occurs the residue that is present only a few times can be removed
        if num_TRP>num_TYR and num_TRP>num_PHE:
            residue = res_TRP.splitlines()
        if num_TYR>num_TRP and num_TYR>num_PHE:
            residue = res_TYR.splitlines()
        if num_PHE>num_TYR and num_PHE>num_TRP:
            residue = res_PHE.splitlines()
   
        capped_residue=''
        for line in residue:
            #Make the alpha carbon a hydrogen atom
            if line[13:15]=='CA':
                line=list(line)
                line[77]='H'
                line[13:15]='H '
                line=''.join(line)
                capped_residue = capped_residue + line + '\n'  
                alpha_number = int(line[9:11])
        for line in residue:
            #List of conditionals to remove unwanted N, C and H
            if int(line[9:11])!=alpha_number and int(line[9:11])!=alpha_number-1 and int(line[9:11])!=alpha_number+1 and line[13:15]!='HN' and line[13:16]!='H  ' and line[13:16]!='HA ':
                capped_residue = capped_residue + line + '\n'
        capped_residue = capped_residue.splitlines()
    return capped_residue

