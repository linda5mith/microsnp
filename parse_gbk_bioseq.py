from Bio import SeqIO
import pandas as pd

#/external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants
#/external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants/Lactobacillus_plantarum.gbk
# /external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants/Faecalibacterium_prausnitzii.gbk
# /external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants/Bacteroides_vulgatus.gbk
# /external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants/Enterococcus_faecalis.gbk
# /external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants/Escherichia_coli.gbk
# /external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants/A.Halo.gbk
# /external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants/Bifidobacterium_longum.gbk
# /external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants/I.Halo.gbk

# Initialise empty lists to hold gbk data
start=[]
stop=[]
protein_id=[]
gene=[]
product=[]
strand=[]

# Create empty dataframe
df=pd.DataFrame({'Start':start,'End':stop,'Protein_ID':protein_id,'Gene':gene,'Product':product,'Strand':strand})

gb_file=open('/external_HDD4/Tom/S.A.3_MouseTrial/Structural_Variants/I.Halo.gbk','r')

for index, record in enumerate(SeqIO.parse(gb_file, "genbank")):
    print(f'Parsing {record.id}')
    for feature in record.features:
        if feature.type=='CDS':
            #print(dir(feature)) get total attributes of feature method
            #print(dir(feature.location))
            start.append(feature.location.start)
            print(feature.location.end)
            stop.append(feature.location.end)
            print(feature.extract)
            print(feature.id)
            print(feature.qualifiers)
            print(feature.qualifiers.get('protein_id'))
            protein_id.append(feature.qualifiers.get('protein_id'))
            print(feature.qualifiers.get('gene'))
            gene.append(feature.qualifiers.get('gene'))
            print(feature.qualifiers.get('product'))
            product.append(feature.qualifiers.get('product'))
            print(feature.strand)
            strand.append(feature.strand)
            print('\n')

df['Start']=start
df['End']=stop
df['Protein_ID']=protein_id
df['Gene']=gene
df['Product']=product
df['Strand']=strand

print(df)
df.to_csv('tom_gbk_parsed/I.Halo_bioseq.csv')