#!/bin/bash

VCFS=$(find /external_HDD4/linda/unc_mouse_trial/snp_pipeline/ -name "*.fltqs.vcf" -type f)
for i in $(ls /data/programs/snpEff/data); do 
    for fle in $VCFS; do 
        basename=${fle##*/}
        basename_no_ext=$(echo "$basename" |cut -d '.' -f 1)
        species_str=$(echo "$basename_no_ext" |cut -d '_' -f 4-)
        path_to_mouse_folder=`dirname "$fle"`
        outpath_html="${path_to_mouse_folder}/${basename_no_ext}.html"
        outpath_vcf="${path_to_mouse_folder}/${basename_no_ext}.fltqs.ann.vcf.gz"
        # Check if species file name is a substring of snpeff database name
        if [[ "$i" == *"$species_str"* ]]; then
            #echo "It's there."
            echo $i
            echo $species_str
            echo "--------------------"
            echo "java -Xmx8g -jar snpEff.jar -v -stats "$outpath_html" "$i" "$fle" > "$outpath_vcf"" >> run_html_variants_cmds_2.sh;
        fi
    done 
done

