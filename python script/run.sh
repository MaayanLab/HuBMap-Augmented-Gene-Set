#!/bin/bash

# Input directory that contains gene set text files
file_dir="/home/bnguy2/correlation_matrix/indir"

for file in "$file_dir"/*.tsv
	do
        command="-g $file -p coexpression_dict.pkl"
        eval "python3 open_feather.py $command"
done