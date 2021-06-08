# snpEff documentation states that reference used should share the same name across alignment, variant calling and variant annotation

# 0_0 very important that file name/folder name matches the name of snpeff database entered reference name e.g. 
# /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1/NC_011993.1_Escherichia/mouse_1_NC_011993.vcf.gz
# /data/programs/snpEff/data/Escherichia_coli_LF82

# for each mouse folder in /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/
# find /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_folder/ -name *.gz
# save these all to variable VCFS
# Extract file prefix as this should be the same name as the snpeff database name


VCFS=$(find /external_HDD4/linda/unc_mouse_trial/test_snp_pipeline/snp_take2/mouse_1 -name "*vcf.gz" -type f)
for fle in $VCFS; do 
    basename=${fle##*/}
    basename_no_ext=$(echo "$basename" |cut -d '.' -f 1)
    path_to_species_folder=`dirname "$fle"`
    base_species_folder=`basename "$path_to_species_folder"`
    path_to_mouse_folder=`dirname "$path_to_species_folder"`
    base_mouse_folder=`basename "$path_to_mouse_folder"`
    outpath_html="${path_to_species_folder}/${base_mouse_folder}_${base_species_folder}.html"
    outpath_vcf="${path_to_species_folder}/${base_species_folder}/${basename_no_ext}.ann.vcf.gz"
    echo "java -jar snpEff.jar -v -stats "$outpath_html" "$base_species_folder" -v stats "$fle" > "$outpath_vcf"" >> run_html_variants_cmds.sh;    
done
