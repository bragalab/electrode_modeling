function Convert_BrainStorm_Coords(SubjectID)
% Convert electrode coordinate given to us from Dr. Kokkinos' automated
% electrode reconstruction pipeline from RAS to voxel, to then be used for
% Braga Lab electrode modeling pipeline
%
CoordDirectory = ['/projects/b1134/processed/fs/',SubjectID,'/', SubjectID,'/elec_recon'];
BrainStorm_CTfile = [CoordDirectory, '/postimpRaw.nii.gz'];

%% load coordinates and volumes
Brainstorm_RAS_coords = readcell(sprintf('%s/%s_brainstorm_RAS.txt', CoordDirectory, SubjectID));
Brainstorm_CT = MRIread(BrainStorm_CTfile);

%% convert from RAS to Voxels
Brainstorm_Voxel_coords = zeros(size(Brainstorm_RAS_coords,1),3);
for i = 1:height(Brainstorm_Voxel_coords)
    %RAS space to voxel space
    converted_coords = Brainstorm_CT.vox2ras \ [cell2mat(Brainstorm_RAS_coords(i,2:4)) 1]';
    Brainstorm_Voxel_coords(i,:) = converted_coords(1:3);  
end

%% save out data
writecell(num2cell(Brainstorm_Voxel_coords), sprintf('%s/%s_brainstorm_VOX.txt', CoordDirectory, SubjectID),...
    'Delimiter', ' ');

end
