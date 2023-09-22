#!/bin/bash
#
#SBATCH --array=0-225 ## number of jobs to run "in parallel"
#SBATCH --account=b1134                	# Our account/allocation
#SBATCH --partition=buyin      		# 'Buyin' submits to our node qhimem0018
#SBATCH --mem=1GB
#SBATCH -t 0:5:00
#SBATCH --job-name makespheres
#SBATCH -o /projects/b1134/analysis/elec2roi/logs/elec2roi_%A_%a.out
#SBATCH -e /projects/b1134/analysis/elec2roi/logs/elec2roi_%A_%a.err

#Written by R.Braga May 2019, parallelized by C. Cyr July 2022

# Create spheres around each electrode, to model sampling volume, for projection to the surface.
# Needs brainmask_coords_0.txt to be in freesurfer elec_recon directory (output of get_depth_coords.m script by A. Kucyi & R. Braga)

#Usage: 
# echo "sbatch /projects/b1134/tools/electrode_modeling/elec2rois_gm_41k_parallel_jobarray.sh ATHUAT 1"

#load modules
module load afni
module load matlab/r2020b
module load freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh 
module load fsl
module load connectome_workbench/1.5.0

#load paths
SUB=$1
radius=$2
radiussq=$((radius * radius))

surf=41k
BASEDIR=/projects/b1134
OUTDIR=$BASEDIR/analysis/elec2roi/${SUB}

SUBJECTS_DIR=$BASEDIR/processed/fs/$SUB

LDIR=$BASEDIR/processed/fs/$SUB/$SUB/elec_recon
ELECDIR=$OUTDIR/elecs_vol_${radius}mm
ELECDIRsurf=$OUTDIR/elecs_surf_${radius}mm_${surf}

SUBJECTS_DIR=$BASEDIR/processed/fs/$SUB
lutt=$BASEDIR/tools/workbench_tools/Label_color_codes/MultipleNetworkColors.txt

mkdir -p $ELECDIR
mkdir -p $ELECDIRsurf

#create blank mask
if [ ! -e $OUTDIR/brainmask_blank.nii.gz ]; then
	fslmaths $LDIR/brainmask.nii.gz -mul 0 $OUTDIR/brainmask_blank
	ln -s $LDIR/brainmask.nii.gz $OUTDIR/brainmask.nii.gz
	ln -s $LDIR $OUTDIR/elec_recon_link
fi

#read in coordinates for all electrodes
OLDIFS=$IFS
IFS=$'\n' GLOBIGNORE='*' command eval  'electrode_list=($(cat $LDIR/brainmask_coords_0_wlabels.txt))'

if (( $SLURM_ARRAY_TASK_ID > ${#electrode_list[@]} )); then
    echo "This array ID is unused because it exceeds the number of electrodes."
    exit
fi

#use one electrode per job
IFS=' ' read -r -a input_args <<< "${electrode_list[$SLURM_ARRAY_TASK_ID]}"
IFS=$OLDIFS

label=${input_args[0]}
x=${input_args[1]}
y=${input_args[2]}
z=${input_args[3]}

echo $label $x $y $z

#Create a volumetric sphere centered on xyz coordinates
echo "Creating volumetric sphere: ${label}_sphere_${x}_${y}_${z}.nii.gz" 
if (( radius < 6 )); then
fslmaths $OUTDIR/brainmask_blank -mul 0 -add 1 -roi $x 1 $y 1 $z 1 0 1 $ELECDIR/${label}_pt_${x}_${y}_${z} -odt float
fslmaths $ELECDIR/${label}_pt_${x}_${y}_${z} -kernel sphere $radius -fmean $ELECDIR/${label}_sphere_${x}_${y}_${z} -odt float
fslmaths $ELECDIR/${label}_sphere_${x}_${y}_${z} -bin $ELECDIR/${label}_sphere_${x}_${y}_${z}
rm $ELECDIR/${label}_pt_${x}_${y}_${z}.nii.gz
chmod 775 $ELECDIR/${label}_sphere_${x}_${y}_${z}.nii.gz
else
3dcalc -a $OUTDIR/brainmask.nii.gz -expr "step(${radiussq}-(i-${x})*(i-${x})-(j-${y})*(j-${y})-(k-${z})*(k-${z}))" -prefix $ELECDIR/${label}_sphere_${x}_${y}_${z}_ball+orig
3dAFNItoNIFTI -prefix $ELECDIR/${label}_sphere_${x}_${y}_${z}.nii.gz $ELECDIR/${label}_sphere_${x}_${y}_${z}_ball+orig
rm $ELECDIR/${label}_sphere_${x}_${y}_${z}_ball+orig*
fi

#project sphere to 41k surface
echo "Projecting sphere to surface ${label}_sphere_${x}_${y}_${z}.nii.gz" 
mri_vol2surf --mov $ELECDIR/${label}_sphere_${x}_${y}_${z}.nii.gz --regheader ${SUB}_41k --hemi lh --projfrac-avg 0 1 0.2 --trgsubject ${SUB}_41k --o $ELECDIRsurf/lh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths.nii --reshape --interp trilinear
mri_vol2surf --mov $ELECDIR/${label}_sphere_${x}_${y}_${z}.nii.gz --regheader ${SUB}_41k --hemi rh --projfrac-avg 0 1 0.2 --trgsubject ${SUB}_41k --o $ELECDIRsurf/rh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths.nii --reshape --interp trilinear

chmod 775 $ELECDIRsurf/*${label}_sphere_${x}_${y}_${z}_${surf}*

mxl=`fslstats $ELECDIRsurf/lh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths.nii -M`
mxr=`fslstats $ELECDIRsurf/rh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths.nii -M`


# Skip electrodes that do not overlap with gray matter
if [ $mxl == 0.000000 ] && [ $mxr == 0.000000 ]; then

echo "# # # # #"
echo "Found non-gray-matter sphere (empty file): $ELECDIRsurf/?h.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths.nii"
echo "SKIPPING this electrode"
echo "# # # # #"
rm $ELECDIRsurf/?h.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths.nii
else

if [ $mxl != 0.000000 ] && [ $mxr == 0.000000 ]; then hemi=lh; fi
if [ $mxl == 0.000000 ] && [ $mxr != 0.000000 ]; then hemi=rh; fi
if [ $mxl != 0.000000 ] && [ $mxr != 0.000000 ]; then hemi=both; fi

fslmaths $ELECDIRsurf/lh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths.nii -bin $ELECDIRsurf/lh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.nii.gz
fslmaths $ELECDIRsurf/rh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths.nii -bin $ELECDIRsurf/rh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.nii.gz

gunzip $ELECDIRsurf/lh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.nii.gz
gunzip $ELECDIRsurf/rh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.nii.gz

matlab -nodisplay -nodesktop -nojvm -r "addpath('$BASEDIR/tools/iProc_archive/iProc/iProc_analysis'); Surf_TaskMap_nii2wb_natsurf_CC('$ELECDIRsurf/lh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.nii','$ELECDIRsurf/rh.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.nii','$ELECDIRsurf','${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.dscalar.nii','$BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/${SUB}_41k_cifti_template.dscalar.nii'); exit"

wb_command -cifti-label-import $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.dscalar.nii $lutt $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.dlabel.nii

wb_command -cifti-separate $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.dlabel.nii COLUMN -label CORTEX_LEFT $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_lh.label.gii -label CORTEX_RIGHT $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_rh.label.gii

if [ $hemi == both ]; then 
wb_command -label-to-border $BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/lh.pial_infl2.surf.gii $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_lh.label.gii $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_lh.border -placement 0.5
wb_command -label-to-border $BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/rh.pial_infl2.surf.gii $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_rh.label.gii $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_rh.border -placement 0.5

elif [ $hemi == lh ]; then 
wb_command -label-to-border $BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/lh.pial_infl2.surf.gii $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_lh.label.gii $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_lh.border -placement 0.5

rm $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_rh.label.gii

elif [ $hemi == rh ]; then 
wb_command -label-to-border $BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/rh.pial_infl2.surf.gii $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_rh.label.gii $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_rh.border -placement 0.5

rm $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin_rh.label.gii
fi

# Cleanup
rm $ELECDIRsurf/?h.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths.nii
rm $ELECDIRsurf/?h.${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.nii

fi



