## visiumv1_spaceranger_v1.sh
	# on server3 # not visiumhd but part of hd project
## 7.16.2025
# Enrico Barrozo, Ph.D. 
# Assistant Professor
# Department of Obstetrics & Gynecology
# Baylor College of Medicine and Texas Children’s Hospital

cd
export PATH=/path/to/local_home/cellranger-9.0.0:$PATH
export PATH=/path/to/local_home/spaceranger-3.1.3:$PATH

## make project directory tree
cd
mkdir visiumhd_mnp
cd visiumhd_mnp
mkdir data
mkdir docs
mkdir results
mkdir scripts

# rm -r cd /path/to/local_home/visiumhd_mnp/results/customref ## remove directory if starting over


cd /path/to/local_home/visiumhd_mnp/results
mkdir customref
# rm -r cd /path/to/local_home/visiumhd_mnp/results/customref ## remove directory if starting over
# mkdir customref
cd /path/to/local_home/visiumhd_mnp/results/customref

################################################# 
spaceranger count --id=S11231-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid3/project/symanski_8117/FastQ_Files \
                   --sample=SYM_69956 \
                   --image=/path/to/local_home/visium/docs/Symanski_manual_alignments/V14F05-312_B1_69956_placenta1-2.tif \
                   --slide=V14F05-312 \
                   --area=B1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/Symanski_manual_alignments/V14F05-312-B1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=true \
                   --localmem=300
      ## make sure are all fastq files are named correctly
                  # Consistent with bcl2fastq/mkfastq, 
                  # e.g. "MySample_S1_L001_R1_001.fastq.gz" "MySample_S1_L001_R2_001.fastq.gz"

spaceranger count --id=S13627-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid3/project/symanski_8117/FastQ_Files \
                   --sample=SYM_69956 \
                   --image=/path/to/local_home/visium/docs/Symanski_manual_alignments/V14F05-312_D1_69957_placenta2-4.tif \
                   --slide=V14F05-312 \
                   --area=D1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/Symanski_manual_alignments/V14F05-312-D1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=true \
                   --localmem=300

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)
spaceranger count --id=S23b-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_134233 \
                   --sample=KA26052 \
                   --image=/path/to/local_home/visium/docs/23b.tif \
                   --slide=V10S21-371 \
                   --area=D1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-371-D1.json \
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

################################################# 
## Run w/  manual alignment in Loupe (labelling tissue manually)

spaceranger count --id=S01-customref \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/GRCh38_and_NC_045512.2_and_EF999921.1_and_JN555585.1_and_LS480640.1_and_NC_001526.4_and_NC_021490.2_and_NZ_CP020432.1 \
                   --fastqs=/path/to/local_storage/raid2/project/aagaard-tillery_133875/25749_01_HCV \
                   --sample=KA_25749 \
                   --image=/path/to/local_home/visium/docs/1.tif \
                   --slide=V10S21-369 \
                   --area=A1 \
                   --loupe-alignment=/path/to/local_home/visium/docs/manual_image_alignments/V10S21-369-A1.json \
                   --localcores=48 \
                   --nosecondary \
                   --create-bam=false \
                   --localmem=300



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


