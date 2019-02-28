# abemus configuration file

# [ General parameters ]

# Output directory
outdir = "/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/ABEMUS"

# abemus main directory [ git clone ]
abemusdir = "/elaborazioni/sharedCO/Home_casiraghi/Prog/abemus/"

# Sample info file 
filesif = "/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_InSilicoDilution_STM2014/SIF_dilution.tsv"

# bases in target regions [ previously obtained by using Rcript utility/PrepareTargetRegions.R ]
targetbp = "/scratch/sharedCO/Casiraghi/Abemus_data/AmpliSeq_AR_STM2015/TargetPositions/"

# PacBam folder
pacbamfolder = "/scratch/sharedCO/Casiraghi/InSilicoDilution/PaCBAM_byChrom/"

# Controls folder
controls_dir = '/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_MultipleSamples_STM2014/ABEMUS/Controls'

# Base Error Model folder
pbem_dir = '/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_MultipleSamples_STM2014/ABEMUS/BaseErrorModel'

# Number of threads
mc.cores = 30

# abemus working steps
# 1. Compute pbem in target positions exploiting germline samples only;
# 2. Compute AF thresholds across coverage leveles by exploiting germline samples only;
# 3. Call and annotate putative somatic SNVs in case samples;

aws = c(3)

# [ 1. Computatation of per-base error model ]

# Bins of coverage into which divide allelic fractions
coverage_binning = 50

# To compute pbem, consider only positions with AF <= af_max_to_compute_pbem  
af_max_to_compute_pbem = 0.2

# To compute pbem, consider only positions with coverage >= coverage_min_to_compute_pbem 
coverage_min_to_compute_pbem = 10

# When compute pbem, count in how many germline samples the position has an AF >= n_pos_af_th
n_pos_af_th = 0.2

# [ 2. Estimation of allelic fraction thresold exploiting germline samples only ]

# To compute AF threshold, consider only positions with AF <= af_max_to_compute_thresholds  
af_max_to_compute_thresholds = 0.2

# To compute AF threshold, consider only positions with coverage >= coverage_min_to_compute_thresholds 
coverage_min_to_compute_thresholds = 10

# [ 2. Estimation of allelic fraction thresold exploiting germline samples only ]

# [ 3. Call SNVs in case samples ]

## Parameters to call somatic SNVs

# Min allelic fraction computed by coverage interval 
# AFbycov = TRUE  : af cut-off computed by coverage intervals
# AFbycov = FALSE : af cut-off computed across all coverage bins
# AFbycov = 0.05 : custom af cut-off, i.e. 0.05

AFbycov = TRUE

# Desired specificity to set AF threshold
spec = 0.9950

# Minimum locus coverage in plasma/tumor sample
mincov = 1

# Minimum number of reads supporting the alternative allele in plasma/tumor sample
minalt = 1

# Minimum locus coverage in germline sample
mincovgerm = 1

# Maximum allelic fraction in germline samples
maxafgerm = 0.01

# minaf corrected table
minaf_corrected_table = "/elaborazioni/sharedCO/Abemus_data_analysis/Analysis_MultipleSamples_STM2014/ABEMUS/Controls/minaf_cov_corrected_10r_0.995_b.RData"

