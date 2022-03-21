import pandas as pd
import numpy as np
import h5py as h5
import re
from sample_list import samp_dict
import pprint
import pickle


species = "human"
version = "11"

single_cell_prob_thresh = 0.5

hf = h5.File(species+"_matrix_v"+version+".h5", "r")
samples = np.array([
   hf["meta"]["samples"]["geo_accession"],
   hf["meta"]["samples"]["singlecellprobability"]
]).T

genes = np.array([
   hf["meta"]["genes"]["genes"]
]).T

samples[:, 0] = samples[:, 0].astype(str) # Convert Geo Accession into strings
genes = genes[:].astype(str) # Convert genes into strings

for organ, gsm_list in samp_dict.items():
    gene_expression = list()
    column_index = 0
    filtered_sample_list = list()
    for sample, sc_prob in samples[:, 0:2]:
        if sample in gsm_list and sc_prob < 0.5:
            filtered_sample_list.append(sample)
            gene_expression.append(hf["data"]["expression"][0:len(genes), column_index])
        column_index += 1

    gene_expression_pd = pd.DataFrame(gene_expression).transpose()
    gene_expression_pd.columns = filtered_sample_list
    gene_expression_pd.index = genes[:,0]
    gene_expression_pd.index.name = "gene"
    #gene_expression_pd.to_pickle("{}_gene_exp.pkl".format(organ)) 
    gene_expression_pd.to_csv("{}_gene_exp.tsv".format(organ), sep="\t")
    print("Done: {}".format(organ))

hf.close()



