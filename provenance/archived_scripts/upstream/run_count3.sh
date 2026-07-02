# run_count3.sh # server3 # 7.17.25
## visiumv1_spaceranger_v1.sh

cd
export PATH=/path/to/local_home/cellranger-9.0.0:$PATH
export PATH=/path/to/local_home/spaceranger-3.1.3:$PATH
cd /path/to/local_home/visiumhd_mnp/results/customref

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S18-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134233 \
                   --sample=KA26047 \
                   --image=/path/to/local_home/visium/docs/18.tif \
                   --slide=V10S21-370 \
                   --area=C1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-370-C1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S19-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134233 \
                   --sample=KA26048 \
                   --image=/path/to/local_home/visium/docs/19.tif \
                   --slide=V10S21-370 \
                   --area=D1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-370-D1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S21-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134233 \
                   --sample=KA26049 \
                   --image=/path/to/local_home/visium/docs/21.tif \
                   --slide=V10S21-371 \
                   --area=A1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-371-A1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S22-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134233 \
                   --sample=KA26050 \
                   --image=/path/to/local_home/visium/docs/22.tif \
                   --slide=V10S21-371 \
                   --area=B1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-371-B1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S23a-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134233 \
                   --sample=KA26051 \
                   --image=/path/to/local_home/visium/docs/23a.tif \
                   --slide=V10S21-371 \
                   --area=C1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-371-C1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300