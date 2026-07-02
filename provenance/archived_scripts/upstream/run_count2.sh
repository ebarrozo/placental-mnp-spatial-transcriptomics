# run_count2.sh # server3 # 7.17.25
## visiumv1_spaceranger_v1.sh




cd
export PATH=/path/to/local_home/cellranger-9.0.0:$PATH
export PATH=/path/to/local_home/spaceranger-3.1.3:$PATH
cd /path/to/local_home/visiumhd_mnp/results/customref

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S03-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_133875/25750_03_HD \
                   --sample=KA_25750 \
                   --image=/path/to/local_home/visium/docs/3.tif \
                   --slide=V10S21-369 \
                   --area=B1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-369-B1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S04-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_133875/25751_04_HCMR \
                   --sample=KA_25751 \
                   --image=/path/to/local_home/visium/docs/4.tif \
                   --slide=V10S21-369 \
                   --area=C1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-369-C1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S15-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_133875/25756_15_HCV \
                   --sample=KA_25756 \
                   --image=/path/to/local_home/visium/docs/15.tif \
                   --slide=V10S21-368 \
                   --area=D1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-368-D1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S16-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134233 \
                   --sample=KA26045 \
                   --image=/path/to/local_home/visium/docs/16.tif \
                   --slide=V10S21-370 \
                   --area=A1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-370-A1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S17-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134233 \
                   --sample=KA26046 \
                   --image=/path/to/local_home/visium/docs/17.tif \
                   --slide=V10S21-370 \
                   --area=B1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-370-B1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300
