#!/bin/bash
# Master Script for creating electrode ROIs, following electrode localization
#
# USAGE: sh /projects/b1134/tools/electrode_modeling/Create_elecrois.sh SubjectID
# If this patient has grid or strip electrodes, this must be run on a compute node
# for parallel processing requirements of the dkystra script (line 22)
#
# BEFORE YOU RUN THIS SCRIPT
# Create_Mgrid.sh must be run and electrodes coordinates must be touched up in BioImageSuite
# 
# Created by Chris Cyr, Braga Lab, February 2024
################################################
module purge
module load matlab/r2020b
module load freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
SubjectID=$1
SUBJECTS_DIR=/projects/b1134/processed/fs/$SubjectID

#correct for brainshift
matlab -batch "addpath(genpath('/projects/b1134/tools/iELVis'));makeIniLocTxtFile('$SubjectID')"
matlab -batch "addpath(genpath('/projects/b1134/tools/iELVis'));dykstraElecPjct('$SubjectID')"
matlab -batch "addpath(genpath('/projects/b1134/tools/iELVis'));get_depth_coords_wshift_zplus1_leptovox('$SubjectID')"

#create final coordinate file used for electrode modeling
sh /projects/b1134/tools/iELVis/Add_eleclabels_to_coords.sh $SubjectID

#create 3mm monopolar electrode ROIs in volume and surface
jobid=($(sbatch /projects/b1134/tools/electrode_modeling/elec2rois_gm_41k_parallel_jobarray.sh $SubjectID 3))
echo "Creating electrode rois: log files at /projects/b1134/analysis/elec2roi/logs/elec2roi*${jobid[-1]}*"

#group ROI files by electrode shaft
jobid1=($(sbatch --dependency=afterany:${jobid[-1]} /projects/b1134/tools/electrode_modeling/Combine_electrode_labels_41k.sh $SubjectID 3))
echo "Grouping electrode rois by shaft/grid/strip: log files at /projects/b1134/analysis/elec2roi/logs/eleccombine*${jobid1[-1]}*"

#create QC document
jobid2=($(sbatch --dependency=afterany:${jobid[-1]} /projects/b1134/tools/electrode_visualization/elec_master.sh $SubjectID))
echo "Creating electrode visualization QC document: log files at /projects/b1134/processed/elec_zoom/logs/eleczoom*${jobid2[-1]}*"
