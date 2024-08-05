function Create_Mgrid(SubjectID)
% Convert electrode coordinate given to us from Dr. Kokkinos' automated
% electrode reconstruction pipeline from RAS to voxel, to then be used for
% Braga Lab electrode modeling pipeline
%
%% add paths
addpath('/projects/b1134/tools/electrode_modeling')
addpath(genpath('/projects/b1134/tools/iELVis/'))
CoordDirectory = ['/projects/b1134/processed/fs/',SubjectID,'/', SubjectID,'/elec_recon'];
Freesurfer_CTFile = [CoordDirectory, '/postInPre.nii.gz'];

%% load coordinates and volumes
Brainstorm_RAS_coords = readcell(sprintf('%s/%s_brainstorm_RAS.txt', CoordDirectory, SubjectID));
Freesurfer_voxel_coords = readcell(sprintf('%s/%s_freesurfer_VOX.txt', CoordDirectory, SubjectID));
Freesurfer_voxel_coords(1,:) = [];
Freesurfer_voxel_coords(:,[2,4,6]) = [];
Freesurfer_CT = MRIread(Freesurfer_CTFile);

%% convert from Voxels to RAS
Freesurfer_RAS_coords = zeros(size(Freesurfer_voxel_coords,1),3);
for i = 1:height(Freesurfer_voxel_coords)
    %RAS space to voxel space
    converted_coords = Freesurfer_CT.vox2ras * [cell2mat(Freesurfer_voxel_coords(i,:)) 1]';
    Freesurfer_RAS_coords(i,:) = converted_coords(1:3);  
end

%% save out mgrid file
elec.elecpos = cell2mat(Freesurfer_voxel_coords);
elec.elecpos(:,3) =  -1 * elec.elecpos(:,3) + size(Freesurfer_CT.vol,3);
for i = 1:height(Freesurfer_RAS_coords)
    %determine electrode name
    shaftinfo = Brainstorm_RAS_coords{i,1}(isstrprop(Brainstorm_RAS_coords{i,1},'alpha'));
    %if there's R and L hemisphere depth electrodes with this name
    if sum(contains(Brainstorm_RAS_coords(:,1), ['L' shaftinfo(2:end)])) > 0 && ...
            sum(contains(Brainstorm_RAS_coords(:,1), ['R' shaftinfo(2:end)]) ) > 0
        
        elec.label{i} = [Brainstorm_RAS_coords{i,1}(1:2), '_', Brainstorm_RAS_coords{i,1}(1), Brainstorm_RAS_coords{i,1}(3:end)];
    else
        elec.label{i} = [Brainstorm_RAS_coords{i,1}(1:2), '_', Brainstorm_RAS_coords{i,1}(3:end)];        
    end
end
write_bioimage_mgrid_CC(sprintf('%s/%s.mgrid',CoordDirectory,SubjectID), elec)

end
