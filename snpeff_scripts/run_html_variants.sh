# java -Xmx8g -jar snpEff.jar -v -stats ex1.html GRCh37.75 protocols/ex1.vcf > protocols/ex1.ann.vcf

# 0_0 very important that file name/folder name matches the name of snpeff database entered reference name e.g. 
# /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1/NC_011993.1_Escherichia/mouse_1_NC_011993.vcf.gz
# /data/programs/snpEff/data/Escherichia_coli_LF82

# path to my vcf should be:
# /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1/NC_011993.1_Escherichia/mouse_Escherichia_coli_LF82.vcf.gz

# for each mouse folder in /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/
# find /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_folder/ -name *.gz
# save these all to array or something
# Extract file prefix as this should be the same name as the snpeff database name



for dir in $(ls ./data); do
    filename = split(dir)
    files =`find . -type d -maxdepth 1 -name "${PREFIX}*"`


    for fle in $(ls /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1 | grep "$dir"); do
        echo "$fle"
        echo "$dir"
        # echo "java -jar snpEff.jar -v -stats $dir.html $dir /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/reference_genomes_db/$dir/$dir.vcf \
        # > /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/reference_genomes_db/$dir/$dir.ann.vcf";
    done;
done;


#java -Xmx8g -jar snpEff.jar -v -stats ex1.html GRCh37.75 protocols/ex1.vcf > protocols/ex1.ann.vcf

java -jar snpEff.jar -v -stats mouse_1_Allobacillus_halotolerans.html \
/external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1/Allobacillus_halotolerans_length/mouse_1_Allobacillus_halotolerans_length.vcf.gz \
> /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1/Allobacillus_halotolerans_length/mouse_1_Allobacillus_halotolerans_length.ann.vcf.gz
