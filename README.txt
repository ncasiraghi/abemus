ABEMUS
Adaptive per Base Error Model for Ultra-deep Sequencing data

• Data preparation

1. dbSNP
2. PaCBAM pileup
3. PaCBAM pileup data splitted by chromosome

• Pre-processing

1. Select BED file and prepare target positions

• Processing and Analysis

1. pbem computation
2. AF threshold estimations
3. Call SNVs
4. SNVs annotations



Scripts, Utilities, Data


todo:

[] in callsnvs.R check the parameters chrom.in.parallel and mc.cores
[] change the way the sample info file is read and information are used to create folder and subfloder during the process
[] create repository with stable version of VCF file to annotated bed targets
[] create repository with bed targets for most used kits
[] include functional annotation in step 4 (Oncotator, snpeff, Annovar)
