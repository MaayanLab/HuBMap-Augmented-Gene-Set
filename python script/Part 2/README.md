## Python Script

### What does it do?
- expression.py

Extract tissue specific expression data from h5 file and filter by RNA-seq samples. 

- cor.R

Normalized using counts adjusted with TMM factors (CTF) and log transforms using asinh. Finally generate pairwise correlation matrix. 

## Usage
- expression.py  
```bash
python3 expression.py
```
- cor.R
```bash
Rscript cor.R <Expression Data Directory>
```
- sample.py

Contains dictionary variable where the key is the tissue and the value is the GSM sample 

### Python Dependencies
1. pandas 
2. numpy 
3. h5py 
4. re
5. pickle
6. Optionally: pprint

### R Dependencies
1. tidyverse
2. edgeR
3. qpgraph
4. DESeq2

## Data
H5 file can be found at [ARCHS4 Download](https://maayanlab.cloud/archs4/download.html)
