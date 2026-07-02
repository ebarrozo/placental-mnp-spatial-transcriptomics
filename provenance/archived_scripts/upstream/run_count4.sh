# run_count4.sh # server3 # 7.17.25
## visiumv1_spaceranger_v1.sh

cd
export PATH=/path/to/local_home/cellranger-9.0.0:$PATH
export PATH=/path/to/local_home/spaceranger-3.1.3:$PATH
cd /path/to/local_home/visiumhd_mnp/results/customref

###################################################### GE Slide 5, change fastq folder and sample names. 

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S20-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134368 \
                   --sample=KA26197 \
                   --image=/path/to/local_home/visium/docs/20.tif \
                   --slide=V10N16-022 \
                   --area=A1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10N16-022-A1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300
################################################# 

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S24-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134368 \
                   --sample=KA26198 \
                   --image=/path/to/local_home/visium/docs/24.tif \
                   --slide=V10N16-022 \
                   --area=B1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10N16-022-B1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300
################################################# 

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S25-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134368 \
                   --sample=KA26199 \
                   --image=/path/to/local_home/visium/docs/25.tif \
                   --slide=V10N16-022 \
                   --area=C1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10N16-022-C1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300
################################################# 

spaceranger count --id=S26-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134368 \
                   --sample=KA26200 \
                   --image=/path/to/local_home/visium/docs/26.tif \
                   --slide=V10N16-022 \
                   --area=D1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10N16-022-D1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300