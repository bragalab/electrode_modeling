function Create_allelecs_labels_41k(SubjectID, radius)
% Creates 41k sphace allelecs sphere file in connectome workbench from 
% monopolar sphere files.
%
% This function only creates .dlabel.nii files, after running this script,
% you must run Create_allelecs_labels_fs6.sh to create the
% accompanying .border files
%% add paths
%set matlab environment
currentdir = pwd;
cd('/projects/b1134/tools/iProc_archive/iProc/iProc_analysis'); 
matlab_setenv
cd(currentdir)

%add paths for gifti and workbench tools
addpath('/projects/b1134/tools/gifti-master');
wbpath = '/projects/b1134/tools/workbench/bin_rh_linux64/wb_command';
%Set electrode ROI path
Lpath = sprintf('/projects/b1134/analysis/elec2roi/%s/elecs_surf_%smm_41k', SubjectID, radius);
Lfiles = dir([Lpath,'/','*','.dlabel.nii']);

%find electrodes to exclude
path_channel_file = dir(fullfile(['/projects/b1134/processed/eegproc/BNI/', SubjectID, '/**', '/Pathologic_Channels.txt']));
path_channels = readcell([path_channel_file(end).folder, '/', path_channel_file(end).name]);
for i = 1:height(path_channels)
    path_channels{i} = [path_channels{i}, '_'];
end
%% Loop through electrode list
for i = 1:length(Lfiles) %for each channel  
    if sum(contains(Lfiles(i).name, path_channels)) == 0
        %open electrode file
        g = ciftiopen([Lpath,'/', Lfiles(i).name],wbpath,1);
        maskids = g.cdata>0;

        if ~exist('allrois', 'var')
            allrois = g.cdata; %initialize all electrode all matrix
        else
            allrois(maskids) = g.cdata(maskids); %add new electrode data
        end    
    end
end
allrois(allrois==0) = nan;

g.cdata = allrois;
OUTDIR = sprintf('/projects/b1134/analysis/elec2roi/%s/grouped_%smm_41k', SubjectID, radius);
mkdir(OUTDIR)
outname = sprintf('%s/allelecs_%smm_41k.dlabel.nii', OUTDIR, radius);
ciftisavereset(g,outname,wbpath)





end


