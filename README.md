# MicroSNP

To install and run microsnp:

```$ git clone https://github.com/linda5mith/microsnp.git```

```$ cd microsnp```

```$ pip install .```

To view all program options:

```$ microsnp -h```

To view prompts for an individual function:

```$ microsnp mpileup_multi -h```

## SNP Pipeline

MicroSNP aims to automate the finding of SNPs and other variants for every species present in a metagenomic data set. 

The objective of this program is to find SNPs present for all species in a .bam file where a .bam file represents a subject’s metagenomic snapshot at a specific time point. 
It is a python wrapper around the command-line bcftools functions. Bcftools 1.7+ is required. 

When doing SNP calling or even working with .bam (samtools alignment files) there is a different file extension that is added/appended after every step in the pipeline. 
The core of the program is to find files with a specified file_extension and then perform a function on them. 

### Requirements
* Bcftools 1.7+
* reference.fasta - file containing reference genomes of all the known species you wish to find SNPs for.  

### 1. Run mpileup

The first step is to run a mpileup of all the samples for a particular subject to create a read coverage summary vcf. 

```microsnp mpileup_multi <path_to_files> <file_extension> <path_to_reference> --<path_to_output>```

```microsnp mpileup_multi /path/to/snp_project/ .bam /path/to/reference.fasta /path/to/snp_project/mpileups```

For the mpileup_multi function to work the files need to be in the following directory structure:

```
subject_1
      |------- subject1_T1.sorted.bam
               subject1_T2.sorted.bam
               subject1_T3.sorted.bam

subject_2
     |---------subject2_T1.sorted.bam
               subject2_T2.sorted.bam
               ......
```
If you don’t have multiple samples for the same subject you can use:

```microsnp mpileup_single <path_to_file> <file_extension> <path_to_reference> --<path_to_output>``` 

Specifying an output path --path_to_output is optional.

bcftools mpileup parameters are hard-coded and optimized for microbial data. If you have a memory limitation change --max-depth paramater. 

### 2. Filter sequences 

<ins>Filter sequences that passed filter defined in mpileup</ins>

```microsnp filter_qual <path_to_files> <file_extension>```

```microsnp filter_qual /path/to/snp_project/mpileups .flt.vcf```

<ins>Filter by DP (FORMAT field)</ins>

```microsnp filter_dp -h```

```microsnp filter_dp <path_to_files> <file_extension> <DP> <expr>```

```microsnp filter_dp /path/to/snp_project/mpileups/ .fltq.vcf 10 MIN```

<ins>Filter by DP & GQ (FORMAT field)</ins>

```microsnp filter_dp_gq <path_to_files> <file_extension> <DP> <expr> <GQ>```

```microsnp filter_dp_gq /path/to/snp_project/mpileups/ .fltq.vcf 10 AVG 50```

### 3. Check compression status (optional) 

Check to make sure each file is properly compressed e.g. is a bg-zipped file
```microsnp htsfile /path/to/snp_project/mpileups .vcf.gz```

### 4. Index vcfs

```microsnp index_vcf /path/to/snp_project/mpileups .vcf.gz```

### 5. Filter individual species from vcf subject files

```microsnp filter_species /path/to/snp_project/mpileups .vcf.gz```

Will create a new vcf file for each species found in individual subject vcfs e.g.

* allobacillus_halotolerans_subject1.fltqs.vcf.gz
* escherichia_coli_subject1.fltqs.vcf.gz

### 6. Index newly created species files

```microsnp index_vcf /path/to/snp_project/mpileups/ fltqs.vcf.qz```

### 7. Normalize species files for merging

```microsnp bcftools_norm /path/to/snp_project/mpileups/ fltqs.vcf.qz path/to/reference.fasta```

### 8. Merge species files

```microsnp bcftools_merge_norms /path/to/snp_project/mpileups/ .norm.fltqs.vcf.qz```

This generates a vcf file which contains all the SNPs or variants for a given species across all your samples. This can be imported into VcfR for further analysis.

## Troubleshooting installation

```pip install .```

```error: error in 'egg_base' option: 'src' does not exist or is not a directory```

pip version needs to be pip3. Check if pip3 is installed:
```which pip3```

If pip3 is installed edit your ~/.bashrc file with the following:

```nano ~/.bashrc```

```alias pip=pip3```



