#!/bin/bash
#SBATCH --array=2-11 ## number of jobs to run "in parallel"
#SBATCH --account=b1134                	# Our account/allocation
#SBATCH --partition=buyin      		# 'Buyin' submits to our node qhimem0018
#SBATCH --mem=10GB
#SBATCH -t 12:00:00
#SBATCH --job-name BuildTable
#SBATCH -o /projects/b1134/analysis/ccyr/logs/BuildTable_%A_%a.out
#SBATCH -e /projects/b1134/analysis/ccyr/logs/BuildTable_%A_%a.err
# add information about distances to each network to SiteInfo excel spreadsheet
#USAGE: sbatch /projects/b1134/tools/electrode_modeling/Append_SiteInfoTable.sh
declare -a pathlist=("/projects/b1134/processed/fs/DQTAWH/DQTAWH/elec_recon" "/projects/b1134/processed/fs/SSYQZJ/SSYQZJ/elec_recon_Surgery1" "/projects/b1134/processed/fs/SSYQZJ/SSYQZJ/elec_recon_Surgery2" "/projects/b1134/processed/fs/PHPKQJ/PHPKQJ/elec_recon" "/projects/b1134/processed/fs/XVFXFI/XVFXFI/elec_recon" "/projects/b1134/processed/fs/KKYNWL/KKYNWL/elec_recon" "/projects/b1134/processed/fs/ZWLWDL/ZWLWDL/elec_recon" "/projects/b1134/processed/fs/TTHMMI/TTHMMI/elec_recon" "/projects/b1134/processed/fs/YKBYHS/YKBYHS/elec_recon" "/projects/b1134/processed/fs/CEWLLT/CEWLLT/elec_recon" "/projects/b1134/processed/fs/XBSGST/XBSGST/elec_recon" "/projects/b1134/processed/fs/DZAEWN/DZAEWN/elec_recon" "/projects/b1134/processed/fs/UYUKJF/UYUKJF/elec_recon" "/projects/b1134/processed/fs/PLLBNH/PLLBNH/elec_recon_Surgery1" "/projects/b1134/processed/fs/DVYZVK/DVYZVK/elec_recon" "/projects/b1134/processed/fs/VPWMYH/VPWMYH/elec_recon" "/projects/b1134/processed/fs/ATHUAT/ATHUAT/elec_recon" "/projects/b1134/processed/fs/S1242/S1242/elec_recon" "/projects/b1134/processed/fs/S19145/S19145/elec_recon" "/projects/b1134/processed/fs/PLLBNH/PLLBNH/elec_recon_Surgery2")

elec_recon_path=${pathlist[$SLURM_ARRAY_TASK_ID]}

module purge
module load R/4.1.1
Rscript /projects/b1134/tools/electrode_modeling/Append_SiteInfoTable.R $elec_recon_path
