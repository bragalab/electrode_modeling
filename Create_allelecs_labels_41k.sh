#!/bin/bash
#
#SBATCH --account=b1134                	# Our account/allocation
#SBATCH --partition=buyin      		# 'Buyin' submits to our node qhimem0018
#SBATCH --mem=8GB
#SBATCH -t 0:30:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name allelecs
#SBATCH -o /projects/b1134/analysis/elec2roi/logs/allelecs_41k_%a_%A.out
#SBATCH -e /projects/b1134/analysis/elec2roi/logs/allelecs_41k_%a_%A.err

#Written by C. Cyr March 2023

# From surface projections of monopolar electrode spheres in 41k space, 
#merge all projections into one allelecs file
#Usage: 
# echo "sbatch /projects/b1134/tools/electrode_modeling/Create_allelecs_labels_41k.sh ATHUAT 3"

module load matlab/r2020b
module load freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh 
module load connectome_workbench/1.5.0
module load fsl
. ${FSLDIR}/etc/fslconf/fsl.sh

SUB=$1
radius=$2

surf=41k
BASEDIR=/projects/b1134
OUTDIR=$BASEDIR/analysis/elec2roi/${SUB}

SUBJECTS_DIR=$BASEDIR/processed/fs/$SUB

LDIR=$BASEDIR/processed/fs/$SUB/$SUB/elec_recon
ELECDIRsurf=$OUTDIR/elecs_surf_${radius}mm_${surf}
GROUPDIRsurf=$OUTDIR/grouped_${radius}mm_${surf}
mkdir -p $ELECDIRsurf
mkdir -p $GROUPDIRsurf

#create allelecs dlabel file
matlab -batch "addpath('$BASEDIR/tools/electrode_modeling'); Create_allelecs_labels_41k('$SUB', '$radius')"

#create border files

#dlabel to label
wb_command -cifti-separate $GROUPDIRsurf/allelecs_${radius}mm_${surf}.dlabel.nii COLUMN -label CORTEX_LEFT $GROUPDIRsurf/allelecs_${radius}mm_${surf}_lh.label.gii -label CORTEX_RIGHT $GROUPDIRsurf/allelecs_${radius}mm_${surf}_rh.label.gii

#label to border
wb_command -label-to-border $BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/lh.pial_infl2.surf.gii $GROUPDIRsurf/allelecs_${radius}mm_${surf}_lh.label.gii $GROUPDIRsurf/allelecs_${radius}mm_${surf}_lh.border -placement 0.5

wb_command -label-to-border $BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/rh.pial_infl2.surf.gii $GROUPDIRsurf/allelecs_${radius}mm_${surf}_rh.label.gii $GROUPDIRsurf/allelecs_${radius}mm_${surf}_rh.border -placement 0.5




