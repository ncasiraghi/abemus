# abemus configuration file

# [ General parameters ]

# Output directory
outdir = "/scratch/sharedCO/Casiraghi/N250_insilico/bgprob_0001_rfac_75"

# abemus main directory [ git clone ]
abemusdir = "/elaborazioni/sharedCO/Home_casiraghi/Prog/abemus/"

# Sample info file 
filesif = "/scratch/sharedCO/Casiraghi/N250_insilico/bgprob_0001_rfac_75/sample_info_file.tsv"

# bases in target regions [ previously obtained by using Rcript utility/PrepareTargetRegions.R ]
targetbp = "/scratch/sharedCO/Casiraghi/N250_insilico/TargetPositions"

# PacBam folder
pacbamfolder = "/CIBIO/sharedCO/PaCBAM/N250/data_N250_insilico_bgprob_0001_rfac_75_byChromosome/"

# Number of threads
mc.cores = 20 

# abemus working steps
# 1. Compute pbem in target positions exploiting germline samples only;
# 2. Compute AF thresholds across coverage leveles by exploiting germline samples only;
# 3. Call and annotate putative somatic SNVs in case samples;

aws = c(1)

# [ 1. Computatation of per-base error model ]

# Bins of coverage into which divide allelic fractions
coverage_binning = 50

# [ 2. Estimation of allelic fraction thresold exploiting germline samples only ]

# [ 3. Call SNVs in case samples ]

# Control Data folder
controls = ""

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
