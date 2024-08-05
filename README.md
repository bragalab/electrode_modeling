# electrode_modeling
For localizing implanted electrodes and modeling them on brain surfaces/images. Creating volumetric spheres surrounding each electrode coordinate (nifti files), and additionally projecting these spheres to an individual subject's brain surface (connectome workbench files).

Create_Mgrid.sh is a master script that can be used for electrode localization, otherwise the patient must be manually localized in Bioimage Suite.

Create_elecrois.sh is a master script for creating electrode spheres and an electrode visualization document, to be run after the completion of electrode localization.

Network_Membership_of_elecs_bipolar_gaussian.m, Build_SiteInfoTable.m, and Append_SiteInfoTable.R/.sh can then be used to build spreadsheets that characterize each electrode (brain region, distance to white matter, proximity to networks, etc.)
