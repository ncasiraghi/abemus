# abemus configuration file

# [ General parameters ]

# Output directory
outdir = "/scratch/sharedCO/Casiraghi/Abemus_data/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/PLASMA_DEDUP/"

# abemus main directory [ git clone ]
abemusdir = "/elaborazioni/sharedCO/Home_casiraghi/Prog/abemus/"

# Sample info file 
filesif = "/elaborazioni/sharedCO/Abemus_data_analysis/Samples_info_files/IPM_PLASMA/SelectionGermlinesPbem/samples_for_pbem.tsv"

# bases in target regions [ previously obtained by using Rcript utility/PrepareTargetRegions.R ]
targetbp = "/scratch/sharedCO/Casiraghi/Abemus_data/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/TargetPositions/"

# PacBam folder
pacbamfolder = "/SPICE-NEW/Casiraghi/PaCBAM/Plasma_Project/DEDUP_mode4_strands/mbq20_mrq20_dedup_bychrom/"

# Controls folder
controls_dir = '/scratch/sharedCO/Casiraghi/Abemus_data/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/PLASMA_DEDUP/Controls/'

# Base Error Model folder
pbem_dir = '/scratch/sharedCO/Casiraghi/Abemus_data/NimbleGen_SeqCapEZ_Exome_v3_Plasma_IPM/PLASMA_DEDUP/BaseErrorModel/'

# Number of threads
mc.cores = 40

# abemus working steps
# 1. Compute pbem in target positions exploiting germline samples only;
# 2. Compute AF thresholds across coverage leveles by exploiting germline samples only;
# 3. Call and annotate putative somatic SNVs in case samples;

aws = c(3)

# [ 1. Computatation of per-base error model ]

# Bins of coverage into which divide allelic fractions
coverage_binning = 50

# To compute pbem, consider only positions with AF <= af_max_to_compute_pbem  
af_max_to_compute_pbem = 1.0

# To compute pbem, consider only positions with coverage >= coverage_min_to_compute_pbem 
coverage_min_to_compute_pbem = 10

# When compute pbem, count in how many germline samples the position has an AF >= n_pos_af_th
n_pos_af_th = 0.15

# [ 2. Estimation of allelic fraction thresold exploiting germline samples only ]

# To compute AF threshold, consider only positions with AF <= af_max_to_compute_thresholds  
af_max_to_compute_thresholds = 0.1

# To compute AF threshold, consider only positions with coverage >= coverage_min_to_compute_thresholds 
coverage_min_to_compute_thresholds = 10

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
