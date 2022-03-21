## Python Script

### What does it do?
- expression.py
Extract tissue specific expression data from h5 file and filter by RNA-seq samples. 

- cor.R
Normalized using counts adjusted with TMM factors (CTF) and log transforms using asinh. Finally generate pairwise correlation matrix. 

## Data
H5 file can be found at [ARCHS4 Download](https://maayanlab.cloud/archs4/download.html)

## Usage
- expression.py  
```bash
python3 expression.py
```
- cor.R
```bash
Rscript cor.R <Expression Data Directory>
```

### Dependencies
1. pandas 
2. numpy 
3. h5py 
4. re
5. pickle
6. Optionally: pprint
