import pandas as pd
import numpy as np
import h5py as h5
import re
import pprint

species = "human"
version = "11"

single_cell_prob_thresh = 0.5

f = h5.File(species+"_matrix_v"+version+".h5", "r")
gse_scprob = np.array([
    f["meta"]["samples"]["series_id"], 
    f["meta"]["samples"]["geo_accession"],
    f["meta"]["samples"]["singlecellprobability"],
    f["meta"]["samples"]["title"]
]).T
f.close()

gse_scprob[:, 0:2] = gse_scprob[:, 0:2].astype(str) # Convert Series Id and Geo Accession into strings
gse_scprob[:, 3] = pd.Series(gse_scprob[:, 3]).str.decode("utf-8") # Decode values

pprint.pprint(gse_scprob[:, 0:2])

# Identify samples w/single cell prob > thresh, Return matching indexes in Numpy Array
single_cell_samp = np.argwhere(gse_scprob[:, 2] > single_cell_prob_thresh)       

# Identify studies corresponding to sc samples, return list of GSE
#print(gse_scprob[single_cell_samp][:,:,2])
single_cell_study = np.unique(gse_scprob[single_cell_samp][:,:,0])               
#print(single_cell_study)

# Boolean mask s.t. {T = bulk, F = sc}
bulk_study_bool = np.isin(gse_scprob[:, 0], single_cell_study, invert = True)    
#print(bulk_study_bool)

# Index corresponding to bulk RNA-seq (T)
bulk_study_idx = np.arange(0, len(bulk_study_bool))[bulk_study_bool]     # Index won't be recorded in array if their is a F in isBulk array 
#pprint.pprint(bulk_study_idx)
#print(bulk_study_idx.shape)

bulk_study_idx2 = np.arange(0, len(bulk_study_bool)) 
pprint.pprint(bulk_study_idx2)
#print(bulk_study_idx2.shape)

# Filtering out scRNA-seq via boolean indexing and appending corresponding h5 index
bulk_study_meta = np.append(gse_scprob[bulk_study_bool],                       
                            bulk_study_idx[:, np.newaxis], axis = 1)             


'''
pd.DataFrame(bulk_study_meta).to_csv(
    species+"_bulk_study_meta.csv", 
    header=["series_id", "geo_accession", "singlecellprobability", "sample_title", "h5_idx"])
'''