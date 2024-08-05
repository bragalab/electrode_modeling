#load packages
if(!require("readxl")) {install.packages("readxl"); require("readxl")}
if(!require("writexl")) {install.packages("writexl"); require("writexl")}
if(!require("stringr")) {install.packages("stringr"); require("stringr")}
if(!require("tidyverse")) {install.packages("tidyverse"); require("tidyverse")}
if(!require("dplyr")) {install.packages("dplyr"); require("dplyr")}
if(!require("devtools")) {install.packages(c("devtools")); require("devtools")}
if(!require("ssbrain")) {install_github("nlanderson9/ssbrain"); require("ssbrain")}
set_wbpath("/projects/b1134/tools/workbench/bin_rh_linux64/wb_command")
source('/projects/b1134/tools/eeganalysis/FUNC_MAPPING/Retrieve_SubjectFiles.R')
############################################## load data
elec_recon_dir <- toString(commandArgs(trailingOnly = TRUE))
SubjectID <- unlist(str_split(elec_recon_dir,'/'))[6]
SurgeryInfo <- unlist(str_split(elec_recon_dir,'/'))[8]
if (length(unlist(str_split(SurgeryInfo,'_')))==3){
  SessionID <- paste0('_',unlist(str_split(SurgeryInfo,'_')))[3]
} else {
  SessionID <- ''
}
Electrode_Info <- read_excel(paste0(elec_recon_dir,'/',SubjectID,'_','SiteInfoTable_bipolar.xlsx'))
Electrode_Info$DisttoDNA <- NA
Electrode_Info$DisttoDNB <- NA
Electrode_Info$DisttoFPNA <- NA
Electrode_Info$DisttoFPNB <- NA
Electrode_Info$DisttodATNA <- NA
Electrode_Info$DisttodATNB <- NA
Electrode_Info$DisttoSAL <- NA
Electrode_Info$DisttoLANG <- NA
Electrode_Info$DisttoUNI <- NA
Electrode_Info$DisttoLowTask <- NA
Electrode_Info$DisttoMidTask <- NA
Electrode_Info$DisttoHighTask <- NA
ParcellationNames <- c('DNA','DNB','FPNA','FPNB','dATNA','dATNB','SAL','LANG','UNI')
ParcellationLists <- list(c('DNA'),c('DNB'),c('FPNA'),c('FPNB'),c('dATNA'),c('dATNB'),c('SAL'),c('LANG'),
                     c('VIS','VISA','VISB','VISC','VISD','VISE','SMOT','SMOTA','SMOTB','SMOTC','AUD','AUDA','AUDB','AUDC'))
TaskNames <- c('LowTask','MidTask','HighTask')
TaskThresholds <- c(0.5, 1, 1.5)
#set path for electrode rois
elecroipath <- paste0('/projects/b1134/analysis/elec2roi/',SubjectID,SessionID)

#set paths to subject's surfaces
surfL <- paste0('/projects/b1134/processed/fs/',SubjectID,'/',SubjectID,'_41k/surf/lh.pial_infl3.surf.gii')
surfR <- paste0('/projects/b1134/processed/fs/',SubjectID,'/',SubjectID,'_41k/surf/rh.pial_infl3.surf.gii')

#load parcellation info
subjectinfo <- Retrieve_SubjectFiles(SubjectID)
labelinfo <- read.table(subjectinfo$labelinfo_file)
my_brain = ss_surf(surfL = surfL,surfR = surfR) + ss_dlabel(filename = subjectinfo$kmeans_file, labels=5, colors='red')
parcellation_info = c(my_brain$dlabel_info$dlabel_info$left,my_brain$dlabel_info$dlabel_info$right)
#load task info
if (subjectinfo$task_file != ''){
  my_brain = ss_surf(surfL = surfL,surfR = surfR) + ss_dscalar(filename = subjectinfo$task_file)
  task_info = c(my_brain$dscalar_info$dscalar_data$left,my_brain$dscalar_info$dscalar_data$right)
} else {
  Electrode_Info[,c('DisttoLowTask','DisttoMidTask','DisttoHighTask')] <- 21
}

for (site in 1:dim(Electrode_Info)[1]){#iterate through each site using increasing ROI sizes until overlapped with a significant parcellation cluster
  contact1 <- str_split(Electrode_Info$ChannelID[site],'-')[[1]][1]
  contact2 <- str_split(Electrode_Info$ChannelID[site],'-')[[1]][2]
  ROIsize <- 1
  while (ROIsize <= 20){#loop through ROI sizes to determine distance to each parcellation   
    
    #load ROI
    roifile1 = list.files(path=paste0(elecroipath,'/elecs_surf_',ROIsize,'mm_41k'), pattern=paste0('^', contact1,'_sphere_','.*','_41k_alldepths_bin.dlabel.nii'), full.names=TRUE)
    if (length(roifile1) > 0){
      roi1 = ss_surf(surfL = surfL,surfR = surfR) +
        ss_dlabel(filename = roifile1, labels=5, colors='red')
      roi1info <- data.frame(location = c(roi1$dlabel_info$dlabel_info$left,roi1$dlabel_info$dlabel_info$right))
    } else {roi1info <- data.frame(location = rep(0,81924))}
    roifile2 = list.files(path=paste0(elecroipath,'/elecs_surf_',ROIsize,'mm_41k'), pattern=paste0('^', contact2,'_sphere_','.*','_41k_alldepths_bin.dlabel.nii'), full.names=TRUE)
    if (length(roifile2) > 0){
      roi2 = ss_surf(surfL = surfL,surfR = surfR) +
        ss_dlabel(filename = roifile2, labels=5, colors='red')
      roi2info <- data.frame(location = c(roi2$dlabel_info$dlabel_info$left,roi2$dlabel_info$dlabel_info$right))
    } else {roi2info <- data.frame(location = rep(0,81924))}      
    roi <- roi1info$location > 0 | roi2info$location > 0
    
    #check which parcellations the ROI overlaps with
    for (parcellation in 1:length(ParcellationNames)){
      if (is.na(Electrode_Info[site,paste0('Distto',ParcellationNames[parcellation])])){ #if distance to parcellation has not already been established
        if (any(parcellation_info[roi] %in% which(labelinfo[,2] %in% ParcellationLists[[parcellation]]))){ #check if this ROI overlaps with parcellation of interest
          Electrode_Info[site,paste0('Distto',ParcellationNames[parcellation])] <- ROIsize
          print(paste('Subject: ', SubjectID, 'Session: ', SessionID,'Site: ', Electrode_Info$ChannelID[site], 'ROI Size: ',ROIsize, 'parcellation: ', ParcellationNames[parcellation]))
        }
      }
    }
    #check which task maps the ROI overlaps with
    for (task in 1:length(TaskNames)){
      if (is.na(Electrode_Info[site,paste0('Distto',TaskNames[task])])){ #if distance to task map has not already been established
        if (any(task_info[roi] > TaskThresholds[task])){ #check if this ROI overlaps with task map
          Electrode_Info[site,paste0('Distto',TaskNames[task])] <- ROIsize
          print(paste('Subject: ', SubjectID, 'Session: ', SessionID,'Site: ', Electrode_Info$ChannelID[site], 'ROI Size: ',ROIsize, 'Task: ', TaskNames[task]))
        }
      }
    }
    if (all(!is.na(Electrode_Info[site,c('DisttoDNA','DisttoDNB','DisttoFPNA','DisttoFPNB',
                                    'DisttodATNA','DisttodATNB','DisttoSAL','DisttoLANG',
                                    'DisttoUNI','DisttoLowTask','DisttoMidTask','DisttoHighTask')]))){
      break
    }
    ROIsize <- ROIsize + 1
  }
}
  
#reformat distances to each parcellation
Electrode_Info$DisttoDNA[is.na(Electrode_Info$DisttoDNA)] <- 21 
Electrode_Info$DisttoDNB[is.na(Electrode_Info$DisttoDNB)] <- 21 
Electrode_Info$DisttoFPNA[is.na(Electrode_Info$DisttoFPNA)] <- 21 
Electrode_Info$DisttoFPNB[is.na(Electrode_Info$DisttoFPNB)] <- 21 
Electrode_Info$DisttodATNA[is.na(Electrode_Info$DisttodATNA)] <- 21 
Electrode_Info$DisttodATNB[is.na(Electrode_Info$DisttodATNB)] <- 21 
Electrode_Info$DisttoSAL[is.na(Electrode_Info$DisttoSAL)] <- 21 
Electrode_Info$DisttoLANG[is.na(Electrode_Info$DisttoLANG)] <- 21 
Electrode_Info$DisttoUNI[is.na(Electrode_Info$DisttoUNI)] <- 21 
Electrode_Info$DisttoLowTask[is.na(Electrode_Info$DisttoLowTask)] <- 21 
Electrode_Info$DisttoMidTask[is.na(Electrode_Info$DisttoMidTask)] <- 21 
Electrode_Info$DisttoHighTask[is.na(Electrode_Info$DisttoHighTask)] <- 21 
write_xlsx(Electrode_Info, paste0(elec_recon_dir,'/',SubjectID,'_','SiteInfoTable_bipolar_appended.xlsx'))

