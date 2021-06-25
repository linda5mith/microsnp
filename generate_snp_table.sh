paste <(bcftools view $1 |\
    awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
      !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
    \
  <(bcftools query -f '[\t%SAMPLE=%GT]\n' $1 |\
    awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
    \
  <(bcftools query -f '[\t%SAMPLE=%GT]\n' $1 |\
    awk 'BEGIN {print "nHomAlt"} {print gsub(/1\|1|1\/1/, "")}') \
    \
  <(bcftools query -f '[\t%SAMPLE=%GT]\n' $1 |\
    awk 'BEGIN {print "nHomRef"} {print gsub(/0\|0|0\/0/, "")}') \
    \
  <(bcftools view $1 | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
    \
  <(bcftools view $1 | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\|1|1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
    \
  <(bcftools view $1 | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|0|0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
    \
  | sed 's/,\t/\t/g' | sed 's/,$//g'

  # output to file 


bcftools view /external_HDD4/linda/unc_mouse_trial/snp_pipeline/mouse_1/Allobacillus_halotolerans_length_2700297/mouse_1_Allobacillus_halotolerans_length_merged.norm.fltq.vcf.gz | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}'