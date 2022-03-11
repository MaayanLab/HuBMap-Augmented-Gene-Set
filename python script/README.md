## Python Script

#### What does it do?
Given a gene correlation matrix in the form of an feather object, find the top 100 coexpressed gene per gene and store in a dictionary.
In the dictionary the key is the gene and value is a list of the top 100 gene with correlation value closest to 1.

Then given a gene set text file where the first column is the cell type and the following columns are the gene set list associated with the cell type,
find the top 100 average genes that coexpresses with the gene set list. Higher weight is given to genes that correlated to multiple genes in the gene set list

This script generates a tsv file where the first column is the cell type and the following columns are the original gene set list concatenated with the top 100 coexpressed genes

#### Dependencies
1. pyarrow.feather as feather
2. os 
3. argparse 
4. pickle
