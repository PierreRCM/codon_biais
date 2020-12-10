from CAI import CAI
from Bio import SeqIO
import glob
import os
import pandas as pd

'''Apply CAI function to all genes files, turn CAI array into a CSV file
   Which contains CAI for each gene'''

reference = [seq.seq for seq in SeqIO.parse("E_C.fasta", "fasta")] # reference array


all_files = glob.glob(os.getcwd() + "/gene_EC_split/*")  # Select all genes files from my path

all_prot_names = []
all_CAI = []

for file in all_files:

    sequence = SeqIO.read(file, "fasta")  # read fasta file which contain nucleotide sequence
    prot_name = file.split("/")[-1]  # recover prot name from file name
    print(prot_name)
    all_prot_names.append(prot_name)
    a_cai = CAI(sequence, reference=reference)  # compute CAI
    print(a_cai)
    all_CAI.append(a_cai)


df = pd.DataFrame(index=all_prot_names, data=all_CAI, columns=["CAI"])  # Turn list into a Dataframe
df.to_csv(os.getcwd()+"/CAI_EC")  # Turn dataframe into a CSV file

