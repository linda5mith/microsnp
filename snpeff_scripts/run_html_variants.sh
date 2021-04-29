# java -Xmx8g -jar snpEff.jar -v -stats ex1.html GRCh37.75 protocols/ex1.vcf > protocols/ex1.ann.vcf


# snpEff documentation states that reference used should share the same name across alignment, variant calling and variant annotation


# 0_0 very important that file name/folder name matches the name of snpeff database entered reference name e.g. 
# /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1/NC_011993.1_Escherichia/mouse_1_NC_011993.vcf.gz
# /data/programs/snpEff/data/Escherichia_coli_LF82

# path to my vcf should be:
# /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1/NC_011993.1_Escherichia/mouse_Escherichia_coli_LF82.vcf.gz

# for each mouse folder in /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/
# find /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_folder/ -name *.gz
# save these all to array or something
# Extract file prefix as this should be the same name as the snpeff database name

VCFS=$(find /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/ -name "*vcf.gz")
for fle in $VCFS; do 
    # echo $fle
    # basename=${fle##*/}
    # echo $basename
    species_folder=`dirname "$fle"`
    base_species_folder=`basename "$species_folder"`
    mouse_folder=`dirname "$species_folder"`
    base_mouse_folder=`basename "$mouse_folder"`
    # echo $base_species_folder
    # echo $base_mouse_folder
    echo "----------------------"
    echo "java -jar snpEff.jar -v -stats "$species_folder"/"$base_species_folder".html "$base_species_folder"\
    -v stats $fle > "$species_folder"/"$base_mouse_folder"_"$base_species_folder".html"
done








#java -Xmx8g -jar snpEff.jar -v -stats ex1.html GRCh37.75 protocols/ex1.vcf > protocols/ex1.ann.vcf

# java -jar snpEff.jar -v -stats mouse_1_Allobacillus_halotolerans.html Allobacillus_halotolerans_length_2700297 /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1/Allobacillus_halotolerans_length/mouse_1_Allobacillus_halotolerans_length.vcf.gz > /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1/Allobacillus_halotolerans_length/mouse_1_Allobacillus_halotolerans_length.ann.vcf.gz

# Allobacillus_halotolerans_length_2700297