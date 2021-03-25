#!/bin/bash

# Set prefix as filename 
# Dont need to iterate through sequence file if I can access chromosome column directly e.g. $3

for i in $(ls $1 | grep ".sorted.bam"); do   
    echo 'FILE:::' $i\n
    #for line in text file
    cat ../species_sequences.txt | while read line 
        do
           echo 'SPECIES:::' $line
        done
    done

# test seq = Enterococcus_phage_A2_length_149431

# sed  -e 's/SN:\([0-9XY]*\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/'

# This command works for the header:
# samtools view -H UNC2FT4158_vs_combined.sam.bam | sed -e 's/SN:\([0-9XY]*\)/SN:test_\1/'

# This command works for the file body! Woop 
# samtools view UNC2FT4158_vs_combined.sam.bam | awk -F'\t' -vOFS='\t' '{ $3 = "TEST!!!=" $3 }1' | head -20

 


