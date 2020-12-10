from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import numpy as np
from Bio import SeqIO
import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

ref = ['ttt', 'tct', 'tat', 'tgt', 'ttc', 'tcc', 'tac', 'tgc', 'tta', 'tca', 'taa', 'tga', 'ttg', 'tcg', 'tag', 'tgg',
       'ctt', 'cct',
       'cat', 'cgt', 'ctc', 'ccc', 'cac', 'cgc', 'cta', 'cca', 'caa', 'cga', 'ctg', 'ccg', 'cag', 'cgg', 'att', 'act',
       'aat', 'agt',
       'atc', 'acc', 'aac', 'agc', 'ata', 'aca', 'aaa', 'aga', 'atg', 'acg', 'aag', 'agg', 'gtt', 'gct', 'gat', 'ggt',
       'gtc', 'gcc',
       'gac', 'ggc', 'gta', 'gca', 'gaa', 'gga', 'gtg', 'gcg', 'gag', 'ggg']


def upper_ref(ref):
    """Turn ref codons into capital letter"""
    ref_upped = []
    for codon in ref:
        codon = codon.upper()
        ref_upped.append(codon)

    return ref_upped


def count_codons(sequence, ref):
    """Input: sequence FASTA, just count each codons in the sequence
       Output: 64-dim array"""
    codons_array = np.zeros([64])
    sequence = ["".join(list(sequence[i:i + 3])) for i in range(0, len(sequence), 3)]
    i = 0
    for codon in sequence:
        try:
            indice_codon = ref.index(codon)
            codons_array[indice_codon] += 1
        except ValueError:
            #H Value error if nucleotide is not ATGC, sometimes N,S ....
            # in this we just skip the codon, it's not serious as it rare event
            i += 1
            print(i)
    return codons_array


ref = upper_ref(ref)

reference = [seq.seq for seq in SeqIO.parse("E_C_full.fasta", "fasta")]

all_files = glob.glob(os.getcwd()+ "/gene_EC_split/*") # Select all gene files
N = len(all_files)
x = np.zeros([N, 64])
prot_names = []

for i, file in enumerate(all_files):

    sequence = SeqIO.read(file, "fasta")
    row = count_codons(sequence, ref)
    prot_names.append(str(file.split("/").pop(-1)))
    x[i] = row/np.sum(row)

x = StandardScaler().fit_transform(x) # Apply transformation to scale data
my_pca = PCA(n_components=3)  # Create a PCA object, choosing 3 PC
pc = my_pca.fit_transform(x) # apply the transformation
# Turn pc score into a dataframe (looks like a excel file)
principalDf = pd.DataFrame(data=pc, columns=['principal component 1', 'principal component 2', 'principal component 3'])

df_weights = pd.DataFrame(my_pca.components_, columns=ref) # turn codons weight into df



color = {0:"r", 1:"g", 2:"b", 3: "y"}
kmean = KMeans(n_clusters=4, random_state=0).fit(pc)  # apply clustering method

prediction = kmean.predict(pc)
label_color = [color[int(i)] for i in prediction]


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(principalDf.iloc[:, 0], principalDf.iloc[:, 1], principalDf.iloc[:, 2], s=1.5, color=label_color)
plt.xlabel("pc1")
plt.ylabel("pc2")
ax.set_zlabel("pc3")
plt.show()



