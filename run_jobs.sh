#!/bin/bash

# We assume running this from the script directory
job_directory=$PWD/.job
file="paramsSmall.txt"

while IFS= read -r line
do
	echo $line
	if [ ${#line} -gt 1 ]
	then
		job_file="${job_directory}/${line// /_}.job"
		echo "#!/bin/bash
			#SBATCH -o .out/${line// /_}.out
			#SBATCH -e .out/${line// /_}.err
			#SBATCH -J ${line// /_}.job
			#SBATCH -p icb_cpu
			#SBATCH -c 1
			#SBATCH --mem=10g
			#SBATCH -t 00-10:00:00
			#SBATCH --nice=10000  # adjusts scheduling priority

			Rscript parameterTuning.r $line" > $job_file
		#sbatch $job_file
	fi
done < $file