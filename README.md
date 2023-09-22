# electrode_modeling
For modeling implanted electrodes (iEEG) onto brain surfaces/images following electrode localization.

Creating volumetric spheres surrounding each electrode coordinate (nifti files), and additionally projecting these spheres to an individual subject's brain surface (connectome workbench files).

elec2rois_gm_41k_parallel_jobarray.sh creates binary spheres (all 1s inside sphere, all 0s outside of sphere)
Combine_electrode_labels_41k.sh groups binary spheres by depth/strip/grid groupings to allow for easier viewing in Connectome Workbench
elec2rois_gm_41k_gaussian_parallel_jobarray.sh creates gaussian spheres (1 at center of sphere, 0 at circumference of sphere, following a gaussian curve)

# Requirements
This package requires Freesurfer, Connectome Workbench, FSL, and matlab/r2020b to be downloaded, and is meant to be run using a SLURM workload manager to utilize for parallel computing resources.

All scripts depend on a text file with voxel-based electrode coordinates located at at /processed/fs/SUB/SUB/elec_recon/brainmask_0_w_labels.txt with the format of

A1 173 172 172

A2 173 171 171

. . . .

. . . .

# Usage
To run via SLURM workload manager:
sbatch elec2rois_gm_41k_parallel_jobarray.sh SubjectID SphereRadius
sbatch Combine_electrode_labels_41k.sh SubjectID SphereRadius
sbatch elec2rois_gm_41k_gaussian_parallel_jobarray.sh SubjectID SphereFullWidthHalfMaximum

To run locally:
sh elec2rois_gm_41k_parallel_jobarray.sh SubjectID SphereRadius
sh Combine_electrode_labels_41k.sh SubjectID SphereRadius
sh elec2rois_gm_41k_gaussian_parallel_jobarray.sh SubjectID SphereFullWidthHalfMaximum

