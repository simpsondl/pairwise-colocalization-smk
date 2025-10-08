#!/bin/bash
#SBATCH --mem=8GB
#SBATCH --time=05:00:00
#SBATCH --output=coloc_logs/log_%A_%a.out

# Submit with
### sbatch --export=start=xx,end=yy submit_colocalization.s

module load R

# array parameters
start=$(( (SLURM_ARRAY_TASK_ID - 1) * 400 + 1 ))
end=$(( (SLURM_ARRAY_TASK_ID) * 400 ))

###### HARDCODED WARNING ######
if [ $end -gt 69256 ]; then
	end=69256
fi

# Loop through analysis
for i in $(seq $start $end)
do
	echo $i
	# get pair for analysis
	awk -v awk_i=$i '$6 == awk_i' pairwise_coloc.txt > analysis_config_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

	# retrieve_summary_stats
	bash colocalization_regions.sh \
	preprocessed_v2 \
	analysis_config_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt \
	coloc_${i}_sumstats_output.txt \
	${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
done

# Clean up
rm analysis_config_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt

for i in $(seq $start $end)
do
	echo "Starting colocalization"
	Rscript run_coloc.R coloc_${i}_sumstats_output.txt
done
