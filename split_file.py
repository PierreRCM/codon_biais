import os
import glob

"""split every gene into file, from a full genome"""


def check_in_file(all_file, name):

    n_all_file = []

    for file in all_file:

        my_file = file.split("/")[-1]
        n_all_file.append(my_file)

    if name in n_all_file:
        return True

    return False


i = 1
gene = ""

with open(os.getcwd()+"/E_C.fasta") as file:
    k=0
    for j, letter in enumerate(file.read()):
        gene += letter
        if j != 0:
            if letter == ">":  # > this letter characterize the beginning of a gene

                all_files = glob.glob(os.getcwd() + "/gene_EC_split_highly_expressed/*")
                # This part select the name of the protein within the the gene file
                name_split = gene.split("[protein=")
                name_part = name_split[1]
                name = name_part.split("]")
                prot = name[0]

                if "/" in prot:
                    prot_split = prot.split("/")
                    prot = "".join(prot_split)

                while check_in_file(all_files, prot):
                    # Most of gene don't have name and are called "hypothetical protein" we need to avoid
                    # overwriting this kind of gene by giving another name
                    # It's a super slow algorythm which can be easely improved
                    k += 1
                    if "_" in prot:
                        a = prot.split("_")
                        a[1] = str(k)
                        prot = "_".join(a)
                    else:
                        prot += "_1"

                with open(os.getcwd()+ "/gene_EC_split_highly_expressed/" + str(k), "w+") as n_f:

                    gene = gene[:-1]
                    n_f.write(gene)

                k += 1
                gene = ">"