## Python Script

### What does it do?
Given a gene correlation matrix in the form of an feather object, find the top 100 coexpressed gene per gene and store in a dictionary.
In the dictionary the key is the gene and value is a list of the top 100 gene with correlation value closest to 1.

Then given a gene set text file where the first column is the cell type and the following columns are the gene set list associated with the cell type,
find the top 100 average genes that coexpresses with the gene set list. Higher weight is given to genes that correlated to multiple genes in the gene set list

This script generates a tsv file where the first column is the cell type and the following columns are the original gene set list concatenated with the top 100 coexpressed genes

## Usage
- Option 1: If coexpression pickle object is used as input, will decrease run time.
```bash
python3 add_genes.py -g <gene set text file> -p <correlation matrix in the form of a pickle object>
```
- Option 2: If no coexpression pickle object is supplied
```bash
python3 add_genes.py -g <GMT file> -f <pair-wise correlation matrix in the form of a feather object> 
```

### Dependencies
1. pyarrow.feather as feather
2. os 
3. argparse 
4. pickle
5. Optionally: pprint

## Running Multiple Gene Sets Text Files
This repository contains a bash script call run.sh where given an absolute path directory that contains all the gene set text files will run add_genes.py on all the text files in the given directory. The directory can be changed by reassigning a new path to the file_dir variable in the bash script. The bash will only work if the coexpression pickle object is supplied. The coexpression pickle object can be generated using Option 2 in Usage.
