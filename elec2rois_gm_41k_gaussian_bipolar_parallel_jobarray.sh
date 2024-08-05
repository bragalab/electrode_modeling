#!/bin/bash
#SBATCH --array=0-225 ## number of jobs to run "in parallel"
#SBATCH --account=b1134                	# Our account/allocation
#SBATCH --partition=buyin      		# 'Buyin' submits to our node qhimem0018
#SBATCH --mem=1GB
#SBATCH -t 0:5:00
#SBATCH --job-name makespheres
#SBATCH -o /projects/b1134/analysis/elec2roi/logs/elec2roi_gauss_%A_%a.out
#SBATCH -e /projects/b1134/analysis/elec2roi/logs/elec2roi_gauss_%A_%a.err

#Written by R.Braga May 2019, adapted by C. Cyr March 2023

# Create gaussian spheres around each electrode, to model sampling volume, and project to the 41k surface.
# Needs brainmask_coords_0.txt to be in freesurfer elec_recon directory (output of get_depth_coords.m script by A. Kucyi & R. Braga)

#Usage: 
# echo "sbatch /projects/b1134/tools/electrode_modeling/elec2rois_gm_41k_gaussian_bipolar_parallel_jobarray.sh ATHUAT 1"

#load modules
module load matlab/r2020b
module load freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh 
module load fsl
module load connectome_workbench/1.5.0

#load paths
SUB=$1
FWHM=$2
FWHM_to_sigma=$(bc -l <<< "${FWHM}*10000/23548")
surf=41k
BASEDIR=/projects/b1134
OUTDIR=$BASEDIR/analysis/elec2roi/${SUB}

SUBJECTS_DIR=$BASEDIR/processed/fs/$SUB

LDIR=$BASEDIR/processed/fs/$SUB/$SUB/elec_recon
ELECDIR=$OUTDIR/elecs_gauss_${FWHM}mm_fwhm_bipolar
ELECDIRsurf=$OUTDIR/elecs_gauss_surf_${FWHM}mm_fwhm_${surf}_bipolar

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

#read in monopolar coordinates for all electrodes
OLDIFS=$IFS
IFS=$'\n' GLOBIGNORE='*' command eval  'coordinate_list=($(cat $LDIR/brainmask_coords_0_wlabels.txt))'

#read in labels for bipolar electrode pairs
IFS=$'\n' GLOBIGNORE='*' command eval  'electrode_list=($(cat $LDIR/${SUB}_bipolarelectrodeNames.txt))'
IFS=$OLDIFS

#exit if the job is unneeded
if (( $SLURM_ARRAY_TASK_ID > ${#electrode_list[@]} )); then
    echo "This array ID is unused because it exceeds the number of electrodes."
    exit
fi

#use one bipolar electrode pair per job
IFS=' ' read -r -a label <<< "${electrode_list[$SLURM_ARRAY_TASK_ID]}"
IFS='-' read -r -a contacts <<< "$label"
IFS=$OLDIFS
label1=${contacts[0]}
label2=${contacts[1]}

#calculate bipolar coordinates
for i in "${!coordinate_list[@]}"; do
   if [[ "${coordinate_list[$i]}" = "${label1} "* ]]; then
	IFS=' ' read -r -a monopolar_coordinates1 <<< "${coordinate_list[$i]}"
	IFS=$OLDIFS
	x1=${monopolar_coordinates1[1]}
	y1=${monopolar_coordinates1[2]}
	z1=${monopolar_coordinates1[3]}
   elif [[ "${coordinate_list[$i]}" = "${label2} "* ]]; then
	IFS=' ' read -r -a monopolar_coordinates2 <<< "${coordinate_list[$i]}"
	IFS=$OLDIFS
	x2=${monopolar_coordinates2[1]}
	y2=${monopolar_coordinates2[2]}
	z2=${monopolar_coordinates2[3]}
   fi
done
x=$(bc -l <<< "scale=2; (${x1}+${x2})/2")
y=$(bc -l <<< "scale=2; (${y1}+${y2})/2")
z=$(bc -l <<< "scale=2; (${z1}+${z2})/2")

echo $label $x $y $z

#Create a gaussian sphere centered on xyz coordinates
echo "Creating gaussian sphere: ${label}_gauss.nii.gz" 
fslmaths $OUTDIR/brainmask_blank -mul 0 -add 1000 -roi $x 1 $y 1 $z 1 0 1 $ELECDIR/${label}_pt -odt float #create target point
fslmaths $ELECDIR/${label}_pt -kernel gauss $FWHM_to_sigma -fmean $ELECDIR/${label}_gauss -odt float #create gaussian sphere around target
minandmax=$(fslstats $ELECDIR/${label}_gauss -R)
IFS=' ' read -ra max <<< "$minandmax"
IFS=$OLDIFS
fslmaths $ELECDIR/${label}_gauss -div ${max[1]} $ELECDIR/${label}_gauss #scale values to max of 1
rm $ELECDIR/${label}_pt.nii.gz
chmod 775 $ELECDIR/${label}_gauss.nii.gz

#project sphere to 41k surface
echo "Projecting gaussian sphere to surface ${label}_gauss.nii.gz" 
mri_vol2surf --mov $ELECDIR/${label}_gauss.nii.gz --regheader ${SUB}_41k --hemi lh --projfrac-avg 0 1 0.2 --trgsubject ${SUB}_41k --o $ELECDIRsurf/lh.${label}_gauss_${surf}_alldepths.nii --reshape --interp trilinear
mri_vol2surf --mov $ELECDIR/${label}_gauss.nii.gz --regheader ${SUB}_41k --hemi rh --projfrac-avg 0 1 0.2 --trgsubject ${SUB}_41k --o $ELECDIRsurf/rh.${label}_gauss_${surf}_alldepths.nii --reshape --interp trilinear
chmod 775 $ELECDIRsurf/*${label}_gauss_${surf}*

mxl=`fslstats $ELECDIRsurf/lh.${label}_gauss_${surf}_alldepths.nii -M`
mxr=`fslstats $ELECDIRsurf/rh.${label}_gauss_${surf}_alldepths.nii -M`

# Skip electrodes that do not overlap with gray matter
if [ $mxl == 0.000000 ] && [ $mxr == 0.000000 ]; then

echo "# # # # #"
echo "Found non-gray-matter sphere (empty file): $ELECDIRsurf/?h.${label}_gauss_${surf}_alldepths.nii"
echo "SKIPPING this electrode"
echo "# # # # #"
rm $ELECDIRsurf/?h.${label}_gauss_${surf}_alldepths.nii
else

matlab -nodisplay -nodesktop -nojvm -r "addpath('$BASEDIR/tools/iProc_archive/iProc/iProc_analysis'); Surf_TaskMap_nii2wb_natsurf_CC('$ELECDIRsurf/lh.${label}_gauss_${surf}_alldepths.nii','$ELECDIRsurf/rh.${label}_gauss_${surf}_alldepths.nii','$ELECDIRsurf','${label}_gauss_${surf}_alldepths.dscalar.nii','$BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/${SUB}_41k_cifti_template.dscalar.nii'); exit"

# Cleanup
rm $ELECDIRsurf/?h.${label}_gauss_${surf}_alldepths.nii

fi


