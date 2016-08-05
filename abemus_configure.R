# abemus configuration file

# [ General parameters ]

# Output directory
outdir = "~/Desktop/temp_working_folder/"

# abemus main directory [ git clone ]
abemusdir = "/Users/nicolacasiraghi/Documents/abemus"

# Sample info file 
filesif = "/Users/nicolacasiraghi/Desktop/temp_working_folder/N250_IPM_samples_with_germlines_info_file.tsv"

# bases in target regions [ previously obtained by using Rcript utility/PrepareTargetRegions.R ]
targetbp = "/elaborazioni/sharedCO/PCF_Project/Analysis/SNVs/DUP/TargetPositions"

# PacBam folder
pacbamfolder = "/scratch/sharedCO/Casiraghi/Plasma_Project/DUP_q10/"

# Number of threads
mc.cores = 2

# abemus working steps
# 1. Compute pbem in target positions exploiting germline samples only;
# 2. Compute AF thresholds across coverage leveles by exploiting germline samples only;
# 3. Call and annotate putative somatic SNVs in case samples;

aws = c(1,2,3)

# [ 1. Estimation of allelic fraction thresold exploiting germline samples only ]

# Bins of coverage into which divide allelic fractions
coverage_binning = 50

# [ 2. Computatation of per-base error model ]

# list of putative somatic SNVs of interest and/or list of known somatic SNVs
snvselection = "/elaborazioni/sharedCO/PCF_Project/Analysis/SNVs/DUP/Results/pmtab_F2_AllSamples.vcf"

# Total number of random positions to build the berr reference random distribution
bprandom = 1000

# compute strand bias at SNP positions
SNP = FALSE

# [ 3. Call SNVs in case samples ]

# Control Data folder
controls = "/elaborazioni/sharedCO/PCF_Project/Analysis/SNVs/DUP/Controls_V7_bq10_covbin50/"

## Parameters to call somatic SNVs

# Min allelic fraction computed by coverage interval 
# AFbycov = TRUE  : af cut-off computed by coverage intervals
# AFbycov = FALSE : af cut-off computed across all coverage bins
# AFbycov = 0.05 : custom af cut-off, i.e. 0.05

AFbycov = 0.02

# Desired specificity to set AF threshold
spec = 0.995

# Minimum locus coverage in plasma/tumor sample
mincov = 10

# Minimum number of reads supporting the alternative allele in plasma/tumor sample
minalt = 1

# Minimum locus coverage in germline sample
mincovgerm = 10

# Maximum allelic fraction in germline samples
maxafgerm = 0.2  
