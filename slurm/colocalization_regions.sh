#!/bin/bash

# script arguments
working_dir=$1
input_fpath=$2
output_fpath=$3
job_id=$4

# if output file exists, remove it
if [ -f "$output_fpath" ]; then
	rm $output_fpath
fi


# summary stats file endings
append="preprocessed.sorted.tsv"

# move to where files are
cd $working_dir

# Loop over an analysis
while IFS= read -r line; do
	# Extract columns
	gwas=$(echo "$line" | awk '{print $1}')
	gwasb=$(echo "$line" | awk '{print $2}')
	chr=$(echo "$line" | awk '{print $3}')
	start=$(echo "$line" | awk '{print $4}')
	end=$(echo "$line" | awk '{print $5}')

	gwas_fname=${gwas}_${append}
	gwasb_fname=${gwasb}_${append}

	# Pull summary stats and ensure they are sorted
	awk -v awk_chr="$chr" -v awk_start="$start" -v awk_end="$end" '
                $1 == awk_chr && $2 >= awk_start && $2 <= awk_end {
                        key=$1 "," $2;
                        if(!seen[key]++) {
                                print;
                        }
                }
        ' ${gwas_fname} > tmp_files/tmp_${job_id}.txt
	sort -k2,2 tmp_files/tmp_${job_id}.txt | cut -f2,5,6 > tmp_files/tmp_${job_id}_sorted.txt
	echo `wc -l tmp_files/tmp_${job_id}_sorted.txt`

	awk -v awk_chr="$chr" -v awk_start="$start" -v awk_end="$end" '
		$1 == awk_chr && $2 >= awk_start && $2 <= awk_end {
                        key=$1 "," $2;
                        if(!seen[key]++) {
                                print;
                        }
                }
        ' ${gwasb_fname} > tmp_files/tmp_${job_id}.txt
        sort -k2,2 tmp_files/tmp_${job_id}.txt | cut -f2,5,6 > tmp_files/tmp_${job_id}_sortedb.txt
        echo `wc -l tmp_files/tmp_${job_id}_sortedb.txt`

	# merge stats across studies
	join -1 1 -2 1 tmp_files/tmp_${job_id}_sorted.txt tmp_files/tmp_${job_id}_sortedb.txt > ${output_fpath}2 && mv ${output_fpath}2 ${output_fpath}

done < "$input_fpath"

# sort output file numerically
sort -k1n ${output_fpath} > ${output_fpath}2 && mv ${output_fpath}2 ${output_fpath}

# Clean up
rm tmp_files/tmp_${job_id}.txt
rm tmp_files/tmp_${job_id}_sorted.txt
rm tmp_files/tmp_${job_id}_sortedb.txt
