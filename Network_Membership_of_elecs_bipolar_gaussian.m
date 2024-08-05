% Created by C. Cyr March 2022
%
% Creates a network membership table, with one row for each bipolar electrode
% Column 1: Bipolar Channel Label
% Column 2: Total number of gray matter vertices
% Column 3: Average tSNR
% Column 4:2:end: Network IDs sorted in order of ranked membership
% Column 5:2:end: Number of vertices in each of those networks
%
function Network_Membership_of_elecs_bipolar_gaussian(SubjectID, spheresize)
%% Add paths
addpath('/projects/b1134/tools/workbench/bin_rh_linux64')
addpath('/projects/b1134/tools/HCP_WB_Tutorial_1.0/fsaverage6')
addpath('/projects/b1134/tools/gifti-master')
wbpath='/projects/b1134/tools/workbench/bin_rh_linux64/wb_command';
currdirr = pwd;
cd('/projects/b1134/tools/iProc_archive/iProc/iProc_analysis'); 
matlab_setenv
cd(currdirr)

if strcmp(SubjectID, 'DQTAWH')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/DQTAWH/REST_test_41k/2mm/vonmises_parcellations/17';
elseif strcmp(SubjectID, 'SSYQZJ')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/SSYQZJ/REST_41k_test/2mm/vonmises_parcellations/14';
elseif strcmp(SubjectID, 'PHPKQJ')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/PHPKQJ/REST/2mm/vonmises_parcellations/18';
elseif strcmp(SubjectID, 'YKBYHS')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/YKBYHS/REST_41k_best/2mm/vonmises_parcellations/20';
elseif strcmp(SubjectID, 'XVFXFI')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/XVFXFI/REST/2mm/vonmises_parcellations/14';    
elseif strcmp(SubjectID, 'CEWLLT')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/CEWLLT/REST/2mm/vonmises_parcellations/14';   
elseif strcmp(SubjectID, 'TTHMMI')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/TTHMMI/REST/2mm/vonmises_parcellations/13';   
elseif strcmp(SubjectID, 'XBSGST')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/XBSGST/REST_41k/2mm/vonmises_parcellations/17';   
elseif strcmp(SubjectID, 'ZWLWDL')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/ZWLWDL/REST/2mm/vonmises_parcellations/17';   
elseif strcmp(SubjectID, 'KKYNWL')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/KKYNWL/REST/2mm/vonmises_parcellations/17';
elseif strcmp(SubjectID, 'DZAEWN')
    vonmises_path = '/projects/b1134/analysis/surfFC/BNI/DZAEWN/REST/2mm/vonmises_parcellations/11';  
end

%% Load Data
%individual level parcellation
fileinfo = split(vonmises_path, '/');
k_number = str2double(fileinfo{end});
vonmises_file = sprintf('%s/%s_k%i_init100_roifsaverage3_vonmises_parcellation.dlabel.nii',...
    vonmises_path, SubjectID, k_number);
vonmises_data = ciftiopen(vonmises_file,wbpath,1);
vonmises_ids = importdata(sprintf('%s/Network_ids.txt', vonmises_path)); %network IDs of each k cluster
vonmises_ids =  cellfun(@(x) strsplit(x, ' '), vonmises_ids, 'UniformOutput', false);
vonmises_ids = vertcat(vonmises_ids{:});

%individual level parcellation confidence
confidence_file = sprintf('%s/%s_k%i_init100_roifsaverage3_vonmises_parcellation_silhouette.dscalar.nii',...
    vonmises_path, SubjectID, k_number); 
confidence_data = ciftiopen(confidence_file,wbpath,1);

%group level parcellation
yeo_file = sprintf('/projects/b1134/analysis/surfFC/BNI/%s/YeoAtlas/Yeo2011_17Networks_N1000.dlabel.nii',SubjectID);
yeo_data = ciftiopen(yeo_file,wbpath,1);

%group level parcellation confidence
yeo_confidence_file = sprintf('/projects/b1134/analysis/surfFC/BNI/%s/YeoAtlas/Yeo2011_17NetworksConfidence_N1000.dscalar.nii',SubjectID);
yeo_confidence_data = ciftiopen(yeo_confidence_file,wbpath,1);

%Set electrode ROI path
Lpath = sprintf('/projects/b1134/analysis/elec2roi/%s/elecs_gauss_surf_%imm_fwhm_41k_bipolar', SubjectID, spheresize);

%Load list of electrode names
Lnamesfile = sprintf('/projects/b1134/processed/fs/%s/%s/elec_recon/brainmask_coords_0_wlabels.txt', ...
    SubjectID, SubjectID);
Lnames=importdata(Lnamesfile);
Lnames = Lnames.textdata;

%create bipolar electrode names
if exist(sprintf('/projects/b1134/processed/fs/%s/%s/elec_recon/%s_bipolarelectrodeNames.txt', ...
        SubjectID, SubjectID, SubjectID), 'file')
    bipolar_Lnames = readcell(sprintf('/projects/b1134/processed/fs/%s/%s/elec_recon/%s_bipolarelectrodeNames.txt', ...
        SubjectID, SubjectID, SubjectID));
else
    shafts = cell(length(Lnames),1);
    for i = 1:length(Lnames)
        shafts{i} = Lnames{i}(1:find(isletter(Lnames{i}), 1, 'last'));
    end 
    bipolar_Lnames = cell(length(Lnames),1);
    for i = length(bipolar_Lnames):-1:2
        if strcmp(shafts{i}, shafts{i-1})
            bipolar_Lnames{i} = sprintf('%s-%s', Lnames{i}, Lnames{i-1});
        end    
    end
    empty_indices = cellfun(@isempty, bipolar_Lnames);
    bipolar_Lnames(empty_indices) = [];
end

%load file names of individual bipolar spheres
Allfiles = arrayfun(@(x) x.name, dir(Lpath), 'UniformOutput', false);
Lfilesindex = contains(Allfiles, '.dscalar.nii');
Lfiles = Allfiles(Lfilesindex);


%% Loop through electrode list
structheight = length(bipolar_Lnames);
Network_info = struct('ChannelID', cell(structheight,1), 'TotalGaussian', cell(structheight,1),...
    'TotalVertices', cell(structheight,1),'Confidence', cell(structheight,1),...
    'DNA_Percentage', cell(structheight,1),'DNB_Percentage', cell(structheight,1),...
    'FPNA_Percentage', cell(structheight,1),'FPNB_Percentage', cell(structheight,1),...
    'dATNA_Percentage', cell(structheight,1),'dATNB_Percentage', cell(structheight,1),...
    'SAL_Percentage', cell(structheight,1),'LANG_Percentage', cell(structheight,1),...
    'UNI_Percentage', cell(structheight,1),'YeoConfidence', cell(structheight,1),...
    'Yeo5_Percentage', cell(structheight,1), 'Yeo6_Percentage', cell(structheight,1),...
    'Yeo7_Percentage', cell(structheight,1), 'Yeo8_Percentage', cell(structheight,1),...
    'Yeo11_Percentage', cell(structheight,1), 'Yeo12_Percentage', cell(structheight,1),...
    'Yeo13_Percentage', cell(structheight,1), 'Yeo14_Percentage', cell(structheight,1),...
    'Yeo15_Percentage', cell(structheight,1), 'Yeo16_Percentage', cell(structheight,1),...
    'Yeo17_Percentage', cell(structheight,1), 'YeoUNI_Percentage', cell(structheight,1));
networks = {'DNA','DNB','FPNA','FPNB','dATNA','dATNB','SAL','LANG','UNI'};
yeo_networks = {'5','6','7','8','11','12','13','14','15','16','17','UNI'};
for i = 1:numel(bipolar_Lnames)
    Network_info(i).ChannelID = bipolar_Lnames{i};
    
    %Load electrode channel label file
    if sum(contains(Lfiles, [bipolar_Lnames{i},'_gauss_'])) == 1
        elecroifile = Lfiles(contains(Lfiles, [bipolar_Lnames{i},'_gauss_']));
    
        %find electrode percentage each network, gaussian weighted by distance
        g = ciftiopen([Lpath,'/', elecroifile{1}],wbpath,1);
        Network_info(i).TotalVertices = sum(g.cdata > 0);
        Network_info(i).TotalGaussian = sum(g.cdata);
        Network_info(i).Confidence = sum(g.cdata.*confidence_data.cdata)/Network_info(i).TotalGaussian;
        Network_info(i).YeoConfidence = sum(g.cdata.*yeo_confidence_data.cdata)/Network_info(i).TotalGaussian;
        network_percentages = zeros(length(networks),1);     
        for j = 1:length(networks)
            if strcmp(networks{j},'UNI')
                network_index = find(matches(vonmises_ids(:,2),{'VIS','VISA','VISB','VISC','VISD','VISE','SMOT','SMOTA','SMOTB','SMOTC','AUD','AUDA','AUDB','AUDC'}));
            else
                network_index = find(matches(vonmises_ids(:,2),networks{j}));
            end
            maskids = ismember(vonmises_data.cdata,network_index);
            network_percentages(j) = sum(g.cdata(maskids))/Network_info(i).TotalGaussian;
        end   
        Network_info(i).DNA_Percentage = network_percentages(1);
        Network_info(i).DNB_Percentage = network_percentages(2);
        Network_info(i).FPNA_Percentage = network_percentages(3);        
        Network_info(i).FPNB_Percentage = network_percentages(4);        
        Network_info(i).dATNA_Percentage = network_percentages(5);        
        Network_info(i).dATNB_Percentage = network_percentages(6);        
        Network_info(i).SAL_Percentage = network_percentages(7);        
        Network_info(i).LANG_Percentage = network_percentages(8);        
        Network_info(i).UNI_Percentage = network_percentages(9);        
        
        yeo_network_percentages = zeros(length(yeo_networks),1);     
        for j = 1:length(yeo_networks)
            if strcmp(yeo_networks{j},'UNI')
                yeo_network_index = 1:4;
            else
                yeo_network_index = str2double(yeo_networks{j});
            end
            yeo_maskids = ismember(yeo_data.cdata,yeo_network_index);
            yeo_network_percentages(j) = sum(g.cdata(yeo_maskids))/Network_info(i).TotalGaussian;
        end   
        Network_info(i).Yeo5_Percentage = yeo_network_percentages(1);
        Network_info(i).Yeo6_Percentage = yeo_network_percentages(2);        
        Network_info(i).Yeo7_Percentage = yeo_network_percentages(3);        
        Network_info(i).Yeo8_Percentage = yeo_network_percentages(4);        
        Network_info(i).Yeo11_Percentage = yeo_network_percentages(5);        
        Network_info(i).Yeo12_Percentage = yeo_network_percentages(6);        
        Network_info(i).Yeo13_Percentage = yeo_network_percentages(7);        
        Network_info(i).Yeo14_Percentage = yeo_network_percentages(8);
        Network_info(i).Yeo15_Percentage = yeo_network_percentages(9);        
        Network_info(i).Yeo16_Percentage = yeo_network_percentages(10);        
        Network_info(i).Yeo17_Percentage = yeo_network_percentages(11);               
        Network_info(i).YeoUNI_Percentage = yeo_network_percentages(12);  
        
    else
        Network_info(i).TotalVertices = NaN;
        Network_info(i).TotalGaussian = NaN;
        Network_info(i).Confidence = NaN;
        Network_info(i).YeoConfidence = NaN;
        Network_info(i).DNA_Percentage = NaN;
        Network_info(i).DNB_Percentage = NaN;
        Network_info(i).FPNA_Percentage = NaN;        
        Network_info(i).FPNB_Percentage = NaN;        
        Network_info(i).dATNA_Percentage = NaN;        
        Network_info(i).dATNB_Percentage = NaN;       
        Network_info(i).SAL_Percentage = NaN;        
        Network_info(i).LANG_Percentage = NaN;       
        Network_info(i).UNI_Percentage = NaN;  
        Network_info(i).Yeo5_Percentage = NaN;
        Network_info(i).Yeo6_Percentage = NaN;        
        Network_info(i).Yeo7_Percentage = NaN;       
        Network_info(i).Yeo8_Percentage = NaN;        
        Network_info(i).Yeo11_Percentage = NaN;        
        Network_info(i).Yeo12_Percentage = NaN;        
        Network_info(i).Yeo13_Percentage = NaN;        
        Network_info(i).Yeo14_Percentage = NaN;
        Network_info(i).Yeo15_Percentage = NaN;        
        Network_info(i).Yeo16_Percentage = NaN;        
        Network_info(i).Yeo17_Percentage = NaN;              
        Network_info(i).YeoUNI_Percentage = NaN;
    end
end

%% Save Results

T = struct2table(Network_info);
outfile=sprintf('%s/Bipolar_gauss_%imm_FWHM_Elec_Network_Membership.csv', vonmises_path, spheresize);
if exist(outfile, 'file')
    delete(outfile)
end
writetable(T, outfile)
end


