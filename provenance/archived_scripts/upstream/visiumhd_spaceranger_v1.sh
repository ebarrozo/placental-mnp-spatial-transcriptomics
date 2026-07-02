## visiumhd_spaceranger_v1.R
	## on server1 
# Enrico Barrozo, Ph.D. 
# Assistant Professor
# Department of Obstetrics & Gynecology
# Baylor College of Medicine and Texas Children’s Hospital

## 7.17.2025

# curl -o spaceranger-4.0.1.tar.gz "https://cf.10xgenomics.com/releases/spatial-exp/spaceranger-4.0.1.tar.gz?Expires=1752742217&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=cTxJCWcqW-DBlMDxVWRbZGKc-uexH9X8nqIb-Q0vJG3JoFHnECywwJ4keTXW9Z3KfTSZC~erGQILKvDKbEWmIpZBPzpoIO9EsHSHKnHJ~NP8j--dJr5ffeKZoU2Wx7sL9CSJHk4LxKFib2monQeNlanrt~0BFw0gN7e3wz1Es3H0ADL1tqfpRDO3Nlm2GL6I8Yp6ePJ6GuOjLIK1kgRIYuleY17v2UkBMC~QowA-naVDTjI2BHM9DI7IockXCUPaei2LloKOXYCglD4xRKn0I6vchS9yH-HEJ3Q~sQtw9YzJv1aIsDpel6jvlfK8aMJID~hdBns8TyuNA1RJqRyRHQ__"


tar -xzvf spaceranger-4.0.1.tar.gz
cd
## make spaceranger commands executable 
export PATH=/path/to/local_home/spaceranger-4.0.1:$PATH
export PATH=/path/to/local_home/cellranger-9.0.1:$PATH



## make project directory tree
cd
mkdir visiumhd_mnp
cd visiumhd_mnp
mkdir data
mkdir docs
mkdir results
mkdir scripts



/path/to/local_home/visiumhd_mnp/scripts
/path/to/local_home/visiumhd_mnp/results
/path/to/local_home/visiumhd_mnp/data/reprep_slide2_H1-X9J8XBC_D1_74312_13164.json

cd /path/to/local_home/visiumhd_mnp/results
mkdir customref
# rm -r cd /path/to/local_home/visiumhd_mnp/results/customref ## remove directory if starting over
# mkdir customref
cd /path/to/local_home/visiumhd_mnp/results/customref
	## not using customref because probebased ##
################################################# 
spaceranger count --id=S11732 \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/refdata-gex-GRCh38-2024-A \
                   --fastqs=/path/to/local_storage/raid2/project/enrico-barrozo_8432_reseq/FastQ_Files \
                   --probe-set=/path/to/local_home/spaceranger-4.0.1/external/tenx_feature_references/targeted_panels/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
                   --filter-probes=true \
                   --sample=BZO_74305 \
                   --cytaimage=/path/to/local_storage/raid2/project/enrico-barrozo_8432_reseq/CytAssist_Images/reprep_slide1_H1-Q6Z4C2C_A1_74305_11732-2.tif \
                   --slide=H1-Q6Z4C2C \
                   --area=A1 \
                   --loupe-alignment=/path/to/local_home/visiumhd_mnp/data/reprep_slide1_H1-Q6Z4C2C_A1_74305_11732-2.json \
                   --localcores=48 \
                   --create-bam=true \
                   --localmem=300
      ## make sure are all fastq files are named correctly
                  # Consistent with bcl2fastq/mkfastq, 
                  # e.g. "MySample_S1_L001_R1_001.fastq.gz" "MySample_S1_L001_R2_001.fastq.gz"

                  	# reprep_slide2_H1-X9J8XBC_A1_74308_29033.tif
spaceranger count --id=S29033 \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/refdata-gex-GRCh38-2024-A \
                   --fastqs=/path/to/local_storage/raid2/project/enrico-barrozo_8432_reseq/FastQ_Files \
                   --probe-set=/path/to/local_home/spaceranger-4.0.1/external/tenx_feature_references/targeted_panels/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
                   --filter-probes=true \
                   --sample=BZO_74308 \
                   --cytaimage=/path/to/local_storage/raid2/project/enrico-barrozo_8432_reseq/CytAssist_Images/reprep_slide2_H1-X9J8XBC_A1_74308_29033.tif \
                   --slide=H1-X9J8XBC \
                   --area=A1 \
                   --loupe-alignment=/path/to/local_home/visiumhd_mnp/data/reprep_slide2_H1-X9J8XBC_A1_74308_29033.json \
                   --localcores=48 \
                   --create-bam=true \
                   --localmem=300

spaceranger count --id=S3057 \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/refdata-gex-GRCh38-2024-A \
                   --fastqs=/path/to/local_storage/raid2/project/enrico-barrozo_8432_reseq/FastQ_Files \
                   --probe-set=/path/to/local_home/spaceranger-4.0.1/external/tenx_feature_references/targeted_panels/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
                   --filter-probes=true \
                   --sample=BZO_74309 \
                   --cytaimage=/path/to/local_storage/raid2/project/enrico-barrozo_8432_reseq/CytAssist_Images/reprep_slide1_H1-Q6Z4C2C_D1_74309_3057.tif \
                   --slide=H1-Q6Z4C2C \
                   --area=D1 \
                   --loupe-alignment=/path/to/local_home/visiumhd_mnp/data/reprep_slide1_H1-Q6Z4C2C_D1_74309_3057.json \
                   --localcores=48 \
                   --create-bam=true \
                   --localmem=300


spaceranger count --id=S13164 \
                   --transcriptome=/path/to/local_storage/raid2/project/customrefs/refdata-gex-GRCh38-2024-A \
                   --fastqs=/path/to/local_storage/raid2/project/enrico-barrozo_8432_reseq/FastQ_Files \
                   --probe-set=/path/to/local_home/spaceranger-4.0.1/external/tenx_feature_references/targeted_panels/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
                   --filter-probes=true \
                   --sample=BZO_74312 \
                   --cytaimage=/path/to/local_storage/raid2/project/enrico-barrozo_8432_reseq/CytAssist_Images/reprep_slide2_H1-X9J8XBC_D1_74312_13164.tif \
                   --slide=H1-X9J8XBC \
                   --area=D1 \
                   --loupe-alignment=/path/to/local_home/visiumhd_mnp/data/reprep_slide2_H1-X9J8XBC_D1_74312_13164.json \
                   --localcores=48 \
                   --create-bam=true \
                   --localmem=300

                   