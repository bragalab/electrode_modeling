function Build_SiteInfoTable(RASDirectory)
% This code needs the ab2v function for the angles and the
% freesurfer_read_surf function for reading surf files. It can also use the
% read_ply function to read ply files but it may not be necessary.

%% add paths
fileinfo = split(RASDirectory,'/');
SubjectID = fileinfo{end-2};
SurgeryInfo = split(fileinfo{end},'_');
if length(SurgeryInfo) == 3
    SessionID = ['_' SurgeryInfo{end}];
else
    SessionID = [];
end
addpath('/projects/b1134/tools/electrode_modeling')
addpath('/projects/b1134/tools/electrode_modeling/CashLab')
SurfDirectory = ['/projects/b1134/processed/fs/',SubjectID,'/', SubjectID,'/surf'];
MRIDirectory = ['/projects/b1134/processed/fs/',SubjectID,'/', SubjectID,'/mri'];

%% load surfaces and volumes
% Reads the left and right hemisphere white matter and pial surf files which allow the subsequent calculations
[pialverticesrh, ~] = freesurfer_read_surf([SurfDirectory,'/rh.pial']);
[whiteverticesrh, whitefacesrh] = freesurfer_read_surf([SurfDirectory,'/rh.white']);

[pialverticeslh, ~] = freesurfer_read_surf([SurfDirectory,'/lh.pial']);
[whiteverticeslh, whitefaceslh] = freesurfer_read_surf([SurfDirectory,'/lh.white']);

%Merges the hemispheres
PialVertices=[pialverticesrh;pialverticeslh];
WhiteVertices=[whiteverticesrh;whiteverticeslh];
WhiteFaces=[whitefacesrh;(whitefaceslh+height(whiteverticesrh))];

%reads subcortical/white matter/gray matter labels for each voxel
SCLabels = MRIread([MRIDirectory,'/aseg.mgz']);
CortexLabels = MRIread([MRIDirectory,'/aparc.a2009s+aseg_edited.mgz']);


%electrode sphere files
volpath = sprintf('/projects/b1134/analysis/elec2roi/%s%s/elecs_vol_2mm', SubjectID, SessionID);
volfiles = dir([volpath, '/', '*sphere*.nii.gz']);
volfiles = arrayfun(@(x) x.name, volfiles, 'UniformOutput', false);

%% load electrode coordinates and labels
% load the RAS coordinates for each monopolar electrode
fid = fopen(sprintf('%s/%s.PIAL', RASDirectory, SubjectID));
coordinfo = textscan(fid, '%s %s %s');
fclose(fid);
coordinfo = [coordinfo{1}, coordinfo{2}, coordinfo{3}];
RAS_coords = str2double(coordinfo(4:end,:));

% load the labels for each monopolar electrode
fid = fopen(sprintf('%s/%s.electrodeNames', RASDirectory, SubjectID));
labelinfo = textscan(fid, '%s %s %s');
fclose(fid);
monopolar_labels = labelinfo{:,1};
monopolar_labels(1:2) = [];

%read in bipolar labels
bipolar_file = sprintf('%s/%s_bipolarelectrodeNames.txt', RASDirectory, SubjectID);
bipolar_labels = readcell(bipolar_file);

%%    calculate bipolar electrode metrics  
Electrode_info = cell2table(cell(height(bipolar_labels), 15));
Electrode_info.Properties.VariableNames = {'ChannelID', 'Coords_x', 'Coords_y',...
    'Coords_z','CoordsBipol_1_x','CoordsBipol_1_y','CoordsBipol_1_z',...
    'CoordsBipol_2_x','CoordsBipol_2_y','CoordsBipol_2_z','InterElectrodeDistance',...
    'DisttoWMBoundary','Orientation', 'TissueType', 'BrainRegion'};

for i = 1:length(bipolar_labels)
    %determine bipolar coordinate
    Electrode_info.ChannelID{i} = bipolar_labels{i};
    contactinfo = split(bipolar_labels{i},'-');
    index1 = matches(monopolar_labels, contactinfo{1});
    index2 = matches(monopolar_labels, contactinfo{2}); 
    Coords = nanmean([RAS_coords(index1,:);RAS_coords(index2,:)]);
    Electrode_info.Coords_x{i} = Coords(1);
    Electrode_info.Coords_y{i} = Coords(2);
    Electrode_info.Coords_z{i} = Coords(3);
    CoordsBipol = [RAS_coords(index1,:);RAS_coords(index2,:)];
    Electrode_info.CoordsBipol_1_x{i} = CoordsBipol(1,1);
    Electrode_info.CoordsBipol_1_y{i} = CoordsBipol(1,2);    
    Electrode_info.CoordsBipol_1_z{i} = CoordsBipol(1,3);    
    Electrode_info.CoordsBipol_2_x{i} = CoordsBipol(2,1);    
    Electrode_info.CoordsBipol_2_y{i} = CoordsBipol(2,2);   
    Electrode_info.CoordsBipol_2_z{i} = CoordsBipol(2,3);    
    
    %calculate inter-electrode distance
    Electrode_info.InterElectrodeDistance{i} = sqrt(sum(bsxfun(@minus, RAS_coords(index1,:), RAS_coords(index2,:)).^2,2));

    % Determine Site Tissue Type
    index1 = ~cellfun(@isempty,regexp(volfiles, ['^',contactinfo{1},'_sphere_', ]));
    elec1 = MRIread([volpath,'/',volfiles{index1}]);
    index2 = ~cellfun(@isempty,regexp(volfiles, ['^',contactinfo{2},'_sphere_', ]));
    elec2 = MRIread([volpath,'/',volfiles{index2}]);
    elecvol = elec1.vol > 0 | elec2.vol > 0; %bipolar elecrode volume
    
    %determine tissue type
    if sum(SCLabels.vol(elec1.vol > 0)) == 0 || sum(SCLabels.vol(elec2.vol > 0)) == 0
        Electrode_info.TissueType{i} = 'Out of Brain';
    else
        current_SClabels = SCLabels.vol(logical(elecvol)); %parcellations that overlap with bipolar electrode volume
        Electrode_info.TissueType{i} = DetermineTissueType(current_SClabels);
    end
    
    %determine brain region
    current_CortexLabels = CortexLabels.vol(elecvol);
    Electrode_info.BrainRegion{i} = DetermineBrainRegion(current_CortexLabels);
    
    %determine distance to gm/wm boundary 
    distancesWM = sqrt(sum(bsxfun(@minus, WhiteVertices, Coords).^2,2)); %determine distance to white each matter boundary vertex
    closestWMvertex_index = find(distancesWM==min(distancesWM)); %white matter boundary vertex closest to electrode
    closestWMfaces = [find(ismember(WhiteFaces(:,1),closestWMvertex_index));... %faces that inclue this white matter boundary vertex
                        find(ismember(WhiteFaces(:,2),closestWMvertex_index));...
                        find(ismember(WhiteFaces(:,3),closestWMvertex_index))];
    distancetofaces = NaN(size(closestWMfaces));
    for face = 1:length(closestWMfaces) %for each face
        facecoordinates = WhiteVertices(WhiteFaces(closestWMfaces(face),:),:);
        F = scatteredInterpolant([facecoordinates(:,1),facecoordinates(:,2)],facecoordinates(:,3));
        F.Method = 'natural';
        F.ExtrapolationMethod = 'none' ;
        [xq,yq] = meshgrid(linspace(min(facecoordinates(:,1)),max(facecoordinates(:,1)),100),...
            linspace(min(facecoordinates(:,2)),max(facecoordinates(:,2)),100));
        zq = F(xq,yq); %interpolate 10000 points onto the face
        distances = [sqrt((xq(:)-Coords(1)).^2 + (yq(:)-Coords(2)).^2 + (zq(:)-Coords(3)).^2);min(distancesWM)]; %calculate distance to each of these points
        distancetofaces(face) = min(distances); %determine nearest distance to this face
    end
    if strcmp(Electrode_info.TissueType(i), 'White Matter') 
        Electrode_info.DisttoWMBoundary{i} = -min(distancetofaces); %nearest distance to any of the faces involved in the nearest white matter voundary vertex
    else
        Electrode_info.DisttoWMBoundary{i} = min(distancetofaces);
    end
    
    %Compute orientation to gray matter column
    closestWM = WhiteVertices(closestWMvertex_index);
    distancesWMGM = sqrt(sum(bsxfun(@minus, PialVertices, closestWM).^2,2));
    closestWMGM = PialVertices(distancesWMGM==min(distancesWMGM),:);
    v1 = CoordsBipol-repmat(CoordsBipol(1,:),2,1);
    v2 = closestWMGM-repmat(closestWM,2,1);
    [~, angleWMGM] = ab2v(v1(2,:), v2(1,:));
    Electrode_info.Orientation{i} = angleWMGM(7);
end
% save out table
outfile = sprintf('%s/%s_SiteInfoTable_bipolar.xlsx', RASDirectory, SubjectID);
if exist(outfile, 'file')
    delete(outfile)
end
writetable(Electrode_info, outfile);

%%    calculate monopolar electrode metrics
Electrode_info = cell2table(cell(height(monopolar_labels), 7));
Electrode_info.Properties.VariableNames = {'ChannelID', 'Coords_x', 'Coords_y',...
    'Coords_z', 'DisttoWMBoundary','TissueType', 'BrainRegion'};

for i = 1:length(monopolar_labels)
    %determine monopolar coordinate
    Electrode_info.ChannelID{i} = monopolar_labels{i};
    Electrode_info.Coords_x{i} = RAS_coords(i,1);
    Electrode_info.Coords_y{i} = RAS_coords(i,2);
    Electrode_info.Coords_z{i} = RAS_coords(i,3);
      
     % Determine Site Tissue Type
    index = ~cellfun(@isempty,regexp(volfiles, ['^',monopolar_labels{i},'_sphere_', ]));
    elec = MRIread([volpath,'/',volfiles{index}]);
    elecvol = elec.vol > 0; %monopolar elecrode volume
    
    %determine tissue type
    current_SClabels = SCLabels.vol(elecvol); %parcellations that overlap with bipolar electrode volume
    Electrode_info.TissueType{i} = DetermineTissueType(current_SClabels);
    
    %determine brain region
    current_CortexLabels = CortexLabels.vol(elecvol);
    Electrode_info.BrainRegion{i} = DetermineBrainRegion(current_CortexLabels);
    
    %determine distance to gm/wm boundary 
    distancesWM = sqrt(sum(bsxfun(@minus, WhiteVertices, Coords).^2,2)); %determine distance to white each matter boundary vertex
    closestWMvertex_index = find(distancesWM==min(distancesWM)); %white matter boundary vertex closest to electrode
    closestWMfaces = [find(ismember(WhiteFaces(:,1),closestWMvertex_index));... %faces that inclue this white matter boundary vertex
                        find(ismember(WhiteFaces(:,2),closestWMvertex_index));...
                        find(ismember(WhiteFaces(:,3),closestWMvertex_index))];
    distancetofaces = NaN(size(closestWMfaces));
    for face = 1:length(closestWMfaces) %for each face
        facecoordinates = WhiteVertices(WhiteFaces(closestWMfaces(face),:),:);
        F = scatteredInterpolant([facecoordinates(:,1),facecoordinates(:,2)],facecoordinates(:,3));
        F.Method = 'natural';
        F.ExtrapolationMethod = 'none' ;
        [xq,yq] = meshgrid(linspace(min(facecoordinates(:,1)),max(facecoordinates(:,1)),100),...
            linspace(min(facecoordinates(:,2)),max(facecoordinates(:,2)),100));
        zq = F(xq,yq); %interpolate 10000 points onto the face
        distances = [sqrt((xq(:)-Coords(1)).^2 + (yq(:)-Coords(2)).^2 + (zq(:)-Coords(3)).^2);min(distancesWM)]; %calculate distance to each of these points
        distancetofaces(face) = min(distances); %determine nearest distance to this face
    end
    if strcmp(Electrode_info.TissueType(i), 'White Matter') 
        Electrode_info.DisttoWMBoundary{i} = -min(distancetofaces); %nearest distance to any of the faces involved in the nearest white matter voundary vertex
    else
        Electrode_info.DisttoWMBoundary{i} = min(distancetofaces);
    end
end
% save out table
outfile = sprintf('%s/%s_SiteInfoTable_monopolar.xlsx', RASDirectory, SubjectID);
if exist(outfile, 'file')
    delete(outfile)
end
writetable(Electrode_info, outfile);

end
