# abemus configuration file

# [ General parameters ]

# Output directory
outdir = "/scratch/sharedCO/Casiraghi/Abemus_data/IPM_N250/"

# abemus main directory [ git clone ]
abemusdir = "/elaborazioni/sharedCO/Home_casiraghi/Prog/abemus/"

# Sample info file 
filesif = "/elaborazioni/sharedCO/Abemus_data_analysis/Samples_info_files/N250_PMI_samples_only_prad_info_file.tsv"

# bases in target regions [ previously obtained by using Rcript utility/PrepareTargetRegions.R ]
targetbp = "/scratch/sharedCO/Casiraghi/Abemus_data/IPM_N250/TargetPositions"

# PacBam folder
pacbamfolder = "/CIBIO/sharedCO/PaCBAM/N250/data_N250_strands_mod4_byChromosome/"

# Controls folder

# controls_dir = '/scratch/sharedCO/Casiraghi/Abemus_data/IPM_N250/Controls_9950/'
# controls_dir = '/scratch/sharedCO/Casiraghi/Abemus_data/IPM_N250/Controls_9990/'
# controls_dir = '/scratch/sharedCO/Casiraghi/Abemus_data/IPM_N250/Controls_9999/'
controls_dir = '/scratch/sharedCO/Casiraghi/Abemus_data/IPM_N250/Controls_9990_10/'

# Base Error Model folder

pbem_dir = '/scratch/sharedCO/Casiraghi/Abemus_data/IPM_N250/BaseErrorModel_covbin10/'


# Number of threads
mc.cores = 30

# abemus working steps
# 1. Compute pbem in target positions exploiting germline samples only;
# 2. Compute AF thresholds across coverage leveles by exploiting germline samples only;
# 3. Call and annotate putative somatic SNVs in case samples;

aws = c(3)

# [ 1. Computatation of per-base error model ]

# Bins of coverage into which divide allelic fractions
coverage_binning = 10

# [ 2. Estimation of allelic fraction thresold exploiting germline samples only ]

# [ 3. Call SNVs in case samples ]

## Parameters to call somatic SNVs

# Min allelic fraction computed by coverage interval 
# AFbycov = TRUE  : af cut-off computed by coverage intervals
# AFbycov = FALSE : af cut-off computed across all coverage bins
# AFbycov = 0.05 : custom af cut-off, i.e. 0.05

AFbycov = TRUE

# Desired specificity to set AF threshold
spec = 0.9990

# Minimum locus coverage in plasma/tumor sample
mincov = 1

# Minimum number of reads supporting the alternative allele in plasma/tumor sample
minalt = 1

# Minimum locus coverage in germline sample
mincovgerm = 1

# Maximum allelic fraction in germline samples
maxafgerm = 0.2
