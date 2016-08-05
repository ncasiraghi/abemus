[0] Prepare dbSNP by chromosome and remove SNPs listed in clinvar.

/elaborazioni/sharedCO/ISMB_ECCB_2015/code/Main/utility/dbSNPPER.R

[1a] Get full information from target regions.

/elaborazioni/sharedCO/ISMB_ECCB_2015/code/Main/utility/PrepareTargetPositions_v2.R

[1b] Prepare PacBam pileup and snvs by chromosome

/elaborazioni/sharedCO/ISMB_ECCB_2015/code/Main/utility/PaCBAM_by_chromosome.R

> Create folder in main working directory to store configure file (ConfigureFiles)

[2] Consider the set of germline samples in order to compute allelic fraction threshold to apply in somatic calls.

/elaborazioni/sharedCO/ISMB_ECCB_2015/code/Main/PM_controls.R
/elaborazioni/sharedCO/ISMB_ECCB_2015/code/Main/PM_controls_configure.R

[3a] Compute look-up table to define AF threshold to apply to call somatic PMs
[3b] Select and annotate positions with AF >= AF cut-off estimated

/elaborazioni/sharedCO/ISMB_ECCB_2015/code/Main/PM_analysis.R

[4] Functional annotation of putative somatic point mutations

[5] Compute Base error model for retained postions and/or a list of known SNVs

/elaborazioni/sharedCO/ISMB_ECCB_2015/code/Main/PM_bem.R
/elaborazioni/sharedCO/ISMB_ECCB_2015/code/Main/PM_bem_configure.R

