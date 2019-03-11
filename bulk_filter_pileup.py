#Import prerequisite libraries to perform installation
import argparse
import numpy as np
parser = argparse.ArgumentParser()

parser.add_argument("mut_pileup_csv", help="Bulked mutant pileup csv file for BSA", type=str)
parser.add_argument("wt_pileup_csv", help="Bulked wildtype pileup csv file for BSA", type=str)

args = parser.parse_args()

mut_pileup = args.mut_pileup_csv
wt_pileup = args.wt_pileup_csv

wt_nmu_pileup=[]
mut_nmu_pileup=''
mut_filt_pileup=''

with open(wt_pileup) as pileup:                       #Read through wild-type bulk SNP pileup calls...
    for line in pileup:
        ref=line.split(',')[2]
        var=line.split(',')[3:7]
        cov=line.split(',')[7]
        if ref=='A':                                      #If reference allele is 'A'.... test is variant is 'G'
            C = line.split(',')[4]
            G = line.split(',')[5]
            T = line.split(',')[6]
            frq=int(G)/int(cov)
            var=[C,G,T]
            max_var=var.index(max(var))                   #Store the position of the most frequent variant (G=1)
            if (max_var==1) and (int(frq)>=0.5):          #and if max_var equals 1 and frequency meets cut off
                wt_nmu_pileup.append(line)                #then store the current pileup up line...
        elif ref=='T':                                #If reference allele is 'T'.... test is variant is 'C'
            A = line.split(',')[3]
            C = line.split(',')[4]
            G = line.split(',')[5]
            frq=int(C)/int(cov)
            var=[A,C,G]
            max_var=var.index(max(var))                   #Store the position of the most frequent variant (C=1)
            if (max_var==1) and (int(frq)>=0.5):          #and if max_var equals 1 and frequency meets cut off
                wt_nmu_pileup.append(line)                #then store the current pileup up line...

with open(mut_pileup) as pileup:                       #Read through wild-type bulk SNP pileup calls...
    for line in pileup:
        chr=line.split(',')[0]
        pos=line.split(',')[1]
        ref=line.split(',')[2]
        var=line.split(',')[3:7]
        cov=line.split(',')[7]
        if ref=='A':                                      #If reference allele is 'A'.... test is variant is 'G'
            C = line.split(',')[4]
            G = line.split(',')[5]
            T = line.split(',')[6]
            frq=int(G)/int(cov)
            var=[C,G,T]
            max_var=var.index(max(var))                   #Store the position of the most frequent variant (G=1)
            if (max_var==1):                              #and if max_var equals 1
                shared=0
                mut_nmu_pileup=mut_nmu_pileup+line        #Storing all NMU SNPs for plotting regardless of filtering here...
                for wt_var in wt_nmu_pileup:
                    wt_chr=wt_var.split(',')[0]
                    wt_pos=wt_var.split(',')[1]
                    if (chr==wt_chr) and (pos==wt_pos):
                        shared=shared+1
                if (shared==0 and frq>=0.99):             #and this allele is fixed and not shared with the wildtype list...
                    print(chr+' '+pos+' G allele is dominant variant: '+str(frq))
                    mut_filt_pileup=mut_filt_pileup+line  #then store the current pileup up line...
        elif ref=='T':                                    #If reference allele is 'T'.... test is variant is 'C'
            A = line.split(',')[3]
            C = line.split(',')[4]
            G = line.split(',')[5]
            frq=int(C)/int(cov)
            var=[A,C,G]
            max_var=var.index(max(var))                   #Store the position of the most frequent variant (G=1)
            if (max_var==1):                              #and if max_var equals 1
                shared=0
                mut_nmu_pileup=mut_nmu_pileup+line        #Storing all NMU SNPs for plotting regardless of filtering here...
                for wt_var in wt_nmu_pileup:
                    wt_chr=wt_var.split(',')[0]
                    wt_pos=wt_var.split(',')[1]
                    if (chr==wt_chr) and (pos==wt_pos):
                        shared=shared+1
                if (shared==0 and frq>=0.99):             #and this allele is fixed and not shared with the wildtype list...
                    print(chr+' '+pos+' C allele is dominant variant: '+str(frq))
                    mut_filt_pileup=mut_filt_pileup+line  #then store the current pileup up line...

file1 = open('mutant_specific_filt_variants.complete.pileup.csv' , "w")
file1.write(mut_filt_pileup)                         #Write wildtype specific variants to new csv file...
file1.close()

file1 = open('mutant_NMU_variants.complete.pileup.csv' , "w")
file1.write(mut_nmu_pileup)                         #Write all putative NMU variants identified to new csv file...
file1.close()