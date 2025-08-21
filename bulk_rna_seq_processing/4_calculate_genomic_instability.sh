#!/bin/bash

source /home/m7huang/miniconda3/etc/profile.d/conda.sh

cd /expanse/lustre/projects/csd728/m7huang/coverage_depth_all

#### bam_names = name of all the files that need to have genome coverage calculated ####

#ls /expanse/lustre/projects/csd728/kfisch/kfisch/horii/final/bams/*.bam > bam_names.txt
#ls /expanse/projects/qstore/csd728/raw_data/placenta_bam/*.bam > bam_names.txt


### determine the genome coverage for each of the samples ###
conda activate bedtools
while IFS= read -r line;
do
    result=$(echo "$line" | sed 's:.*/::')
    echo "$result"
    
    genomeCoverageBed -ibam $line -split > /expanse/lustre/projects/csd728/m7huang/coverage_depth_all/${result}_coverage.txt
    
done < bam_names.txt

### calculate the coverage for depth > 50 ###
#### file_list = name of all the coverage calculation files ####

ls | grep -v 'file_list.txt' | grep -v 'bam_names.txt' > file_list.txt

echo -e "File Name \t Coverage" >> output.txt

while IFS= read -r input_file;
do
    num_bases=$(awk 'BEGIN {
        split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM", chroms)
        for (i in chroms) {
            allowed[chroms[i]] = 1
        }
    }
    {
        # Check if the depth is less than 51
        if (allowed[$1] && $2 < 51) {
            sum += $3
        }
    }
    END {
        print sum
    }' "$input_file")
    
    genome_total=$(awk 'BEGIN {
        split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM", chroms)
        for (i in chroms) {
            allowed[chroms[i]] = 1
        }
    }
    {
        # Extract the genome total value for each chromosome
        if ($1 in allowed && !seen[$1]++) {
            sum += $4
        }
    }
    END {
        print sum
    }' "$input_file")

    result=$(awk "BEGIN {print 1 - ($num_bases / $genome_total)}")

    echo -e "$input_file \t $result" >> output.txt
    
done < file_list.txt
 
exit 0
