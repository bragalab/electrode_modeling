function BOLD_RSFC_ElectrodeROI_Matrix(RASDirectory, iTaskDirectory)
%% Add toolboxes
addpath('/projects/b1134/tools/workbench/bin_rh_linux64')
addpath('/projects/b1134/tools/HCP_WB_Tutorial_1.0/fsaverage6')
addpath('/projects/b1134/tools/gifti-master')
wbpath='/projects/b1134/tools/workbench/bin_rh_linux64/wb_command';
currdirr = pwd;
cd('/projects/b1134/tools/iProc_archive/iProc/iProc_analysis'); 
matlab_setenv
cd(currdirr)

%% load fMRI file paths
lh_surflist_file = [iTaskDirectory, '/lh_surflist_row.txt'];
fid = fopen(lh_surflist_file, 'r');
fileinfo = textscan(fid, '%s');
lh_surflist = fileinfo{1};

rh_surflist_file = [iTaskDirectory, '/rh_surflist_row.txt'];
fid = fopen(rh_surflist_file, 'r');
fileinfo = textscan(fid, '%s');
rh_surflist = fileinfo{1};

%% load electrode ROI file paths
%Set electrode ROI path
pathinfo = split(RASDirectory,'/');
SubjectID = pathinfo{end-1};
sessioninfo = split(pathinfo{end},'_');
if length(sessioninfo) == 3
    SessionID = sessioninfo{end};
    Lpath = sprintf('/projects/b1134/analysis/elec2roi/%s_%s/elecs_gauss_surf_10mm_fwhm_41k_bipolar', SubjectID, SessionID);
else
    SessionID = [];
    Lpath = sprintf('/projects/b1134/analysis/elec2roi/%s/elecs_gauss_surf_10mm_fwhm_41k_bipolar', SubjectID);
end

%Load list of electrode names
Lnamesfile = sprintf('%s/brainmask_coords_0_wlabels.txt', RASDirectory);
Lnames=importdata(Lnamesfile);
Lnames = Lnames.textdata;

%create bipolar electrode names
if exist(sprintf('%s/%s_bipolarelectrodeNames.txt', RASDirectory, SubjectID), 'file')
    bipolar_Lnames = readcell(sprintf('%s/%s_bipolarelectrodeNames.txt', ...
        RASDirectory, SubjectID));
else
    shafts = cell(length(Lnames),1);
    for run = 1:length(Lnames)
        shafts{run} = Lnames{run}(1:find(isletter(Lnames{run}), 1, 'last'));
    end 
    bipolar_Lnames = cell(length(Lnames),1);
    for run = length(bipolar_Lnames):-1:2
        if strcmp(shafts{run}, shafts{run-1})
            bipolar_Lnames{run} = sprintf('%s-%s', Lnames{run}, Lnames{run-1});
        end    
    end
    empty_indices = cellfun(@isempty, bipolar_Lnames);
    bipolar_Lnames(empty_indices) = [];
end

%load file names of individual bipolar spheres
Allfiles = arrayfun(@(x) x.name, dir(Lpath), 'UniformOutput', false);
Lfilesindex = contains(Allfiles, 'dscalar.nii');
Lfiles = Allfiles(Lfilesindex);

%% Create Correlation Matrix
for run = 1:length(lh_surflist)
    %load vertex-wise bold timeseries
    input = erase(lh_surflist{run},'_sm2');
    input_series = MRIread(input);
    lh_vertex_series = single(transpose(reshape(input_series.vol,...
        size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4))));

    input = erase(rh_surflist{run},'_sm2');
    input_series = MRIread(input);
    rh_vertex_series = single(transpose(reshape(input_series.vol,...
        size(input_series.vol, 1) * size(input_series.vol, 2) * size(input_series.vol, 3), size(input_series.vol, 4))));

    wb_vertex_series = [lh_vertex_series, rh_vertex_series];
    
    %create electrode ROI-wise bold timesieres
    elecROI_series = NaN(height(wb_vertex_series),length(bipolar_Lnames));
    for channel = 1:length(bipolar_Lnames)
        %Load electrode channel label file
        if sum(contains(Lfiles, [bipolar_Lnames{channel},'_gauss_'])) == 1
            elecroifile = Lfiles(contains(Lfiles, [bipolar_Lnames{channel},'_gauss_']));
            g = ciftiopen([Lpath,'/', elecroifile{1}],wbpath,1);
            %normalize weights and transpose vector
            g.cdata = (g.cdata/sum(g.cdata))';
            %created weighted average timeseries
            elecROI_series(:,channel) = sum(wb_vertex_series .* repmat(g.cdata,height(wb_vertex_series),1),2);
        end
    end
    
    % normalize series (note that series are now of dimensions: T x N)
    elecROI_series = bsxfun(@minus, elecROI_series, mean(elecROI_series, 1));
    elecROI_series = bsxfun(@times, elecROI_series, 1./sqrt(sum(elecROI_series.^2, 1)));

    sbj_corr_mat = elecROI_series' * elecROI_series;
    if(run == 1)
        sbj_z_mat = StableAtanh(sbj_corr_mat); % fisher-z transform
    else
        sbj_z_mat = sbj_z_mat + StableAtanh(sbj_corr_mat);
    end
end
sbj_z_mat = sbj_z_mat/length(lh_surflist);
corr_mat = tanh(sbj_z_mat);
OUTPATH = strrep(iTaskDirectory,'2mm','0mm');
mkdir(OUTPATH)
save(sprintf('%s/corr_mat_elecROI.mat', OUTPATH),'corr_mat','bipolar_Lnames')

%% Create Dconn File
%g = ciftiopen(sprintf('/projects/b1134/processed/fs/%s/%s_41k/surf/%s_41k_cifti_template.dscalar.nii',...
%    SubjectID,SubjectID,SubjectID),wbpath,1);
%g.cdata = zeros(81924,length(bipolar_Lnames));
%weight_matrix = zeros(size(g.cdata));
%for channel = 1:length(bipolar_Lnames)
%    if sum(contains(Lfiles, [bipolar_Lnames{channel},'_gauss_'])) == 1
%        elecroifile = Lfiles(contains(Lfiles, [bipolar_Lnames{channel},'_gauss_']));
%        elecroi = ciftiopen([Lpath,'/', elecroifile{1}],wbpath,1);
%        weight_matrix(:,channel) = elecroi.cdata/sum(elecroi.cdata);
%    end
%end
%weight_matrix = weight_matrix ./ repmat(sum(weight_matrix,2),1,length(bipolar_Lnames));

%for channel1 = 1:length(bipolar_Lnames)
%    channel_values = zeros(81924,length(bipolar_Lnames));
%    for channel2 = 1:length(bipolar_Lnames)
%        maskids = weight_matrix(:,channel2) > 0;
%        channel_values(:,channel2) = maskids * corr_mat(channel1,channel2);
%    end
%    g.cdata(:,channel1) = sum(weight_matrix .* channel_values,2);
%end

%ciftisavereset(g,sprintf('%s/corr_mat_elecROI.dscalar.nii', OUTPATH),...
%    '/projects/b1134/tools/workbench/bin_rh_linux64/wb_command')

