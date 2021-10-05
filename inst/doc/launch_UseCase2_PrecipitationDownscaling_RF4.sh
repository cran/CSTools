#!/bin/bash
#BSUB -W 6:00
#BSUB -n 16 
#BSUB -M 7000
#BSUB -J RainFARM_Downsc
#BSUB -o /esarchive/scratch/nperez/CSTools_manuscript/v20210603/lsf-%J.out
#BSUB -e /esarchive/scratch/nperez/CSTools_manuscript/v20210603/lsf-%J.err

module load R
Rscript /esarchive/scratch/nperez/git/cstools/inst/doc/UseCase2_PrecipitationDownscaling_RainFARM_RF4.R
