'''
Usage: 
Option 1: If coexpression pickle object is used as input
python3 adding_genes.py -g <gene set text file> -p <correlation matrix in the form of a pickle object>

Option 2: If no coexpression pickle object is supplied
python3 adding_genes.py -g <GMT file> -f <pair-wise correlation matrix in the form of a feather object> 

Given a gene correlation matrix in the form of an feather object, find the top 100 coexpressed gene per gene and store in a dictionary.
In the dictionary the key is the gene and value is a list of the top 100 gene with correlation value closest to 1.

Then given a gene set text file where the first column is the cell type and the following columns are the gene set list associated with the cell type,
find the top 100 average genes that coexpresses with the gene set list. Higher weight is given to genes that correlated to multiple genes in the gene set list

This script generates a tsv file where the first column is the cell type and the following columns are the original gene set list concatenated with the top 100 coexpressed genes
'''
import pyarrow.feather as feather
import os 
import argparse 
import pickle
#import pprint
#import sys

class CommandLine():
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Input correlation feather object and gene set library'
        )

        self.parser.add_argument(
            '-g', type=str, help='Specify gene set library', required=True
        )

        self.parser.add_argument(
            '-f', type=str, help='Specify correlation feather object', required=False
        )

        self.parser.add_argument(
            '-p', type=str, help='Specify correlation pickle object', required=False
        )

        self.args = self.parser.parse_args()

def buildCoexpressionDict(read_df):
    gene_list = read_df.columns.to_list()
    coexpression_dict = dict()

    for row_number in range(1, len(read_df)):

        # Iterate through each row of the data frame
        subset_df = read_df.iloc[row_number-1:row_number, 0:len(read_df.columns)]

        # Sum the columns in the subsetted data frame
        sum_df = subset_df.sum()

        # Sort the data frame in descesding order
        sorted_df = subset_df[sum_df.sort_values(ascending=False).index[:]]

        # Find the top 100 correlation values (values closest to 1)
        top_df = sorted_df.iloc[0, 1:101]

        # Find the associated genes with the top 100 correlation values and output to list
        top_genes_list = top_df.index.to_list()
        top_genes_value_list = top_df.to_list()

        gene_correlation_value = list()
        for index in range(0, len(top_genes_list)):
            gene_correlation_value.append([top_genes_list[index], top_genes_value_list[index]])
        
        coexpression_dict[gene_list[row_number - 1]] = gene_correlation_value
    
    return(coexpression_dict)

# Assign gene with a correlation weight based on the number of correlation and the average correlation values
def averageList(correlation_value_list):
        return((sum(correlation_value_list) / len(correlation_value_list)) + len(correlation_value_list) - 1)

def returnTopGenes(gene_marker_list, coexpression_dict):
    # Dictionary where the key is gene and the value is a list of correlation value
    #print(gene_marker_list)
    gene_correlation_dict = dict()
    count_gene_marker_not_coexpression_dict = 0

    for gene_marker in gene_marker_list:
        if gene_marker in coexpression_dict:
            for gene, correlation_value in coexpression_dict[gene_marker]:
                if gene not in gene_correlation_dict:
                    gene_correlation_dict[gene] = list()
                    gene_correlation_dict[gene].append(correlation_value)
                else:
                    gene_correlation_dict[gene].append(correlation_value)
        else:
            count_gene_marker_not_coexpression_dict += 1

    # Dictionary where the key is gene and the value is the average correlation value
    gene_correlation_average_dict = dict()

    for gene, correlation_value_list in gene_correlation_dict.items():
        gene_correlation_average_dict[gene] = averageList(correlation_value_list)

    #pprint.pprint(gene_correlation_dict)
    #pprint.pprint(gene_correlation_average_dict)

    top_100_gene_list = sorted(gene_correlation_average_dict.items(), key=lambda x:-x[1])[:100]
    return top_100_gene_list, count_gene_marker_not_coexpression_dict

def parseTSV(file_name):
    fh = open(file_name, "r")
    gene_marker_dict = dict()
    for line in fh:
        tab_list = line.strip("\n").split("\t")
        try:
            cell_type = tab_list[0]
            gene_list = tab_list[1]
            gene_marker_dict[cell_type] = gene_list
        except:
            pass

    return(gene_marker_dict)

def main():
    myCommandLine = CommandLine()
    gene_marker_dict = parseTSV(myCommandLine.args.g)
    coexpression_dict = dict()
    
    if not myCommandLine.args.p:
        read_df = feather.read_feather(myCommandLine.args.f)
        coexpression_dict = buildCoexpressionDict(read_df)
    else:
        pickle_dict = open(myCommandLine.args.p, "rb")
        coexpression_dict = pickle.load(pickle_dict)

    final_dict = dict()
    counts = 0
    for cell_type, gene_set_str in gene_marker_dict.items():
        gene_set_list = gene_set_str.split(", ")
        top_coexpressed_genes_list, counts = returnTopGenes(gene_set_list, coexpression_dict)
        final_gene_set_list = gene_set_list + top_coexpressed_genes_list
        final_dict[cell_type] = final_gene_set_list 
        counts += counts

    #print(counts) # Display the number of genes in gene set not found in correlation matrix
    #pprint.pprint(final_dict)

    # Change the current working directory to output directory
    os.chdir('{}/outdir'.format(os.getcwd()))

    pre_process_organ_output_name = myCommandLine.args.g[:-4] 
    pre_process_organ_output_name = pre_process_organ_output_name.split("/")
    organ_output_name = pre_process_organ_output_name[len(pre_process_organ_output_name) - 1]
    output_name = "augment_{}.tsv".format(organ_output_name)
    
    fh = open(output_name, "w")
    # Write out
    for cell_type, gene_list in final_dict.items():
        fh.write("{}".format(cell_type))
        for gene in gene_list:
            if type(gene) == str:
                fh.write("\t{}".format(gene))
            else:
                fh.write("\t{}".format(gene[0]))
        fh.write("\n")

if __name__ == "__main__":
	main()

'''
# Unit test
#gene_marker_list = ['A1BG', 'A1CF', 'A2M', 'A2ML1'] #======================================
#gene_marker_dict = {'cell_type': ['A1BG', 'A1CF', 'A2M', 'A2ML1']} #======================================

#pprint.pprint(coexpression_dict)
'''

'''
# For creating coexpression dictionary in the form of a pickle object
# Used to save computation time
# create a binary pickle file 
f = open("coexpression_dict.pkl", "wb")

# write the python object (dict) to pickle file
pickle.dump(coexpression_dict, f)

# close file
f.close()

sys.exit()
'''
