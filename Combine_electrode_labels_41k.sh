#!/bin/sh

#SBATCH --account=b1134                	# Our account/allocation
#SBATCH --partition=buyin      		# 'Buyin' submits to our node qhimem0018
#SBATCH --mem=4GB
#SBATCH -t 00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name eleccombine
#SBATCH -o /projects/b1134/analysis/elec2roi/logs/eleccombine_%a_%A.out
#SBATCH -e /projects/b1134/analysis/elec2roi/logs/eleccombine_%a_%A.err

#Usage: 
# echo "sbatch /projects/b1134/tools/electrode_modeling/Combine_electrode_labels_41k.sh SSYQZJ 2"

module load matlab/r2018a

module load freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh 

SUB=$1
radius=$2

surf=41k

BASEDIR=/projects/b1134
OUTDIR=$BASEDIR/analysis/elec2roi/${SUB}

LDIR=$BASEDIR/processed/fs/$SUB/$SUB/elec_recon
ELECDIR=$OUTDIR/elecs_vol_${radius}mm
ELECDIRsurf=$OUTDIR/elecs_surf_${radius}mm_${surf}

SUBJECTS_DIR=$BASEDIR/processed/fs/$SUB

COORDIR=$SUBJECTS_DIR/$SUB/elec_recon

#coords=$COORDIR/brainmask_coords_pialvox_0.txt # 166
coords=$COORDIR/brainmask_coords_0_wlabels.txt # 166

lutt=$BASEDIR/tools/workbench_tools/Label_color_codes/MultipleNetworkColors.txt

GDIR=$OUTDIR/grouped_${radius}mm_$surf
mkdir -p $GDIR
cd $GDIR

ln -s $coords $GDIR/Elec_coords.txt

# # #
#Loop over each line in the brainmask_coord_0_wlabels
OLDIFS=$IFS

IFS=$'\n'
set -f

count=0

for i in $(cat $GDIR/Elec_coords.txt)
do

IFS=$OLDIFS #reset internal field separator
set +f

x=$(echo $i | awk '{print $2}')
y=$(echo $i | awk '{print $3}')
z=$(echo $i | awk '{print $4}')

label=$(echo $i | awk '{print $1}')
electrode=$(echo $label | tr -d '0123456789')

echo $x $y $z

count=$((count+1))

mkdir -p $GDIR/${electrode}

# Skip empty (non-gray-matter) files
if [ -e $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.dlabel.nii ]; then

wb_command -cifti-label-import $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.dlabel.nii $lutt $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.dlabel.nii -drop-unused-labels

wb_command -set-map-names $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.dlabel.nii -map 1 $label

wb_command -cifti-separate $ELECDIRsurf/${label}_sphere_${x}_${y}_${z}_${surf}_alldepths_bin.dlabel.nii COLUMN -label CORTEX_LEFT $GDIR/${electrode}/${label}_sphere_${x}_${y}_${z}_${surf}_lh.label.gii -label CORTEX_RIGHT $GDIR/${electrode}/${label}_sphere_${x}_${y}_${z}_${surf}_rh.label.gii

else

echo "File not found - Skipping"

fi
done

cd $GDIR

for lon in `ls --ignore=*.txt --ignore=*.gii --ignore=*.border`;
do

cd $GDIR/$lon

pwd

touch $GDIR/${lon}/lh_${lon}_contact_list.txt
ls $GDIR/${lon}/${lon}?_sphere_*_*_*_${surf}_lh.label.gii > $GDIR/${lon}/lh_${lon}_contact_list.txt
ls $GDIR/${lon}/${lon}??_sphere_*_*_*_${surf}_lh.label.gii >> $GDIR/${lon}/lh_${lon}_contact_list.txt

touch $GDIR/${lon}/rh_${lon}_contact_list.txt
ls $GDIR/${lon}/${lon}?_sphere_*_*_*_${surf}_rh.label.gii > $GDIR/${lon}/rh_${lon}_contact_list.txt
ls $GDIR/${lon}/${lon}??_sphere_*_*_*_${surf}_rh.label.gii >> $GDIR/${lon}/rh_${lon}_contact_list.txt

rline=$(cat lh_${lon}_contact_list.txt)
rray=($rline)
for i in ${rray[@]} # NEED TO GROUP IN ORDER
do
echo " -label " >> tmp_$lon
echo $i >> tmp_$lon
done
wb_command -label-merge $GDIR/${lon}_electrode_lh.label.gii `cat tmp_$lon | tr -d '\n'`
wb_command -label-to-border $BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/lh.pial_infl2.surf.gii $GDIR/${lon}_electrode_lh.label.gii $GDIR/${lon}_electrode_lh.border -placement 0.5

rm -r tmp_$lon

rline=$(cat rh_${lon}_contact_list.txt)
rray=($rline)
for i in ${rray[@]} # NEED TO GROUP IN ORDER
do
echo " -label " >> tmp_$lon
echo $i >> tmp_$lon
done

wb_command -label-merge $GDIR/${lon}_electrode_rh.label.gii `cat tmp_$lon | tr -d '\n'`
wb_command -label-to-border $BASEDIR/processed/fs/${SUB}/${SUB}_41k/surf/rh.pial_infl2.surf.gii $GDIR/${lon}_electrode_rh.label.gii $GDIR/${lon}_electrode_rh.border -placement 0.5

rm -r tmp_$lon

done
