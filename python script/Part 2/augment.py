import pandas as pd
import os 
import argparse 
import pickle
import sys
import re
#import pprint

class CommandLine():
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description='Input correlation feather object and gene set library'
        )

        self.parser.add_argument(
            '-g', type=str, help='Specify directory pathway for gmt data', required=True
        )

        self.parser.add_argument(
            '-c', help='Specify directory pathway for correlation data', required=False
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

# Create Gene Marker Dictionary. Key = Cell Type, Value = Gene list
def parse_gmt_tsv(file_name, organ_name):
    fh = open(file_name, "r")
    gene_marker_dict = dict()
    for line in fh:
        tab_list = line.strip("\n").split("\t")
        try:
            cell_type = tab_list[0]
            gene_list = tab_list[1]
            gene_marker_dict[cell_type] = gene_list
            #gene_marker_dict[cell_type+":"+organ_name.capitalize()] = gene_list
        except:
            pass

    return(gene_marker_dict)

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def main():
    myCommandLine = CommandLine()
    cur_dir = os.getcwd()

    # Essential Dictionary
    gene_marker_dict = dict()
    coexpression_dict = dict()
    organ_list = list()
    
    # Build Coexpression Dictionary and store as a pickle object. Key = Organ, Value = Dictionary(Key = Gene, Value = Top 100 Coexpressed Genes)
    if not myCommandLine.args.p:
        print("Building Coexpression Dictionary...")
        corr_dir = myCommandLine.args.c
        corr_dir_list = os.listdir(corr_dir)
        os.chdir(corr_dir)

        for file in corr_dir_list:
            organ_name = '_'.join(file.split("_")[:-1]).lower() # Standardize organ names 
            print("For... {}: {}".format(file, organ_name))
            coexpression_dict[organ_name] = buildCoexpressionDict(pd.read_csv(file, sep='\t'))
            organ_list = coexpression_dict.keys()
            print("{} - {}".format(organ_list, len(coexpression_dict[organ_name])))

        # Create pickle file
        print("Writing to pickle")
        os.chdir(cur_dir)
        f = open("coexpression_dict.pkl", "wb")
        # Write coexpression_dict to pickle file
        pickle.dump(coexpression_dict, f)
        # Close pickle file
        f.close()
    else:
        pickle_dict = open(myCommandLine.args.p, "rb")
        coexpression_dict = pickle.load(pickle_dict)

        coexpression_dict["bone_marrow"] = coexpression_dict["bone"]
        del coexpression_dict["bone"]
        organ_list = coexpression_dict.keys()
        print("Organs in coexpression dictionary: {}".format(organ_list))
    
    # Build Gene Marker Dictionary. Key = Organ, Value = Dictionary(Key = Cell Type, Value = Gene Marker)
    print("Building Gene Marker Dictionary...")
    gmt_dir = myCommandLine.args.g
    gmt_dir_list = os.listdir(gmt_dir)

    gmt_list = ['_'.join(file.split("_")[:-2]).lower() for file in gmt_dir_list]
    mapping_list = intersection(organ_list, gmt_list)

    organs_not_in_coexpression_dict = list()

    for file in gmt_dir_list:
        organ_name = '_'.join(file.split("_")[:-2]).lower() # Standardize organ names
        if organ_name in mapping_list:
            gene_marker_dict[organ_name] = parse_gmt_tsv(cur_dir+"/"+gmt_dir+"/"+file, organ_name)
        else:
            organs_not_in_coexpression_dict.append(organ_name)
    print("GMT files organs not in coexpression dictionary: {}".format(organs_not_in_coexpression_dict))
    
    # Add top genes to gene set associated with a cell type
    print("Finding top genes for...")
    final_dict = dict()
    counts = 0
    for organ, gene_dict in gene_marker_dict.items():
        print(organ)
        for cell_type, gene_set_str in gene_dict.items():
            gene_set_list = gene_set_str.split(";")    # Separate Genes in Gene Set
            #pprint.pprint(coexpression_dict[organ])
            top_coexpressed_genes_list, counts = returnTopGenes(gene_set_list, coexpression_dict[organ])
            final_gene_set_list = gene_set_list + top_coexpressed_genes_list
            if organ not in final_dict.keys():
                final_dict[organ] = dict()
            final_dict[organ][cell_type] = set(final_gene_set_list)     # Remove duplicated genes in gene set 
            counts += counts

    #print(counts) # Display the number of genes in gene set not found in coexpression dictionary
    #pprint.pprint(final_dict)
    
    # Check if output directory exist
    output_path = '{}/outdir'.format(os.getcwd())
    if os.path.isdir(output_path):
        # Change the current working directory to output directory
        os.chdir(output_path)
    else:
        # Create output directory if it does not exist
        os.mkdir(output_path)
        os.chdir(output_path)
    
    #pre_process_organ_output_name = myCommandLine.args.g[:-4] 
    #pre_process_organ_output_name = pre_process_organ_output_name.split("/")
    #organ_output_name = pre_process_organ_output_name[len(pre_process_organ_output_name) - 1]
    output_name = "CellMarker_Augmented_2022_table.tsv"
    
    # Write out
    print("Writing Output...")

    fh = open(output_name, "w")
    for organ, cell_gene_dict in final_dict.items():
        for cell_type, gene_list in cell_gene_dict.items():
            fh.write("{}\t".format(cell_type))
            for count, gene in enumerate(gene_list):
                # Case for first gene
                if count == 0:
                    if type(gene) == str:
                        fh.write("{}".format(gene))
                    else:
                        fh.write("{}".format(gene[0]))
                
                # Case for the rest of the genes
                if count > 0 and type(gene) == str:
                    fh.write(";{}".format(gene))
                elif count > 0:
                    fh.write(";{}".format(gene[0]))
            fh.write("\n")

if __name__ == "__main__":
	main()

'''
# Unit test
#gene_marker_list = ['A1BG', 'A1CF', 'A2M', 'A2ML1'] #======================================
#gene_marker_dict = {'cell_type': ['A1BG', 'A1CF', 'A2M', 'A2ML1']} #======================================
#pprint.pprint(coexpression_dict)
'''
