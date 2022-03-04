################################################################################
################################################################################
#####GEOMETRIC MORPHOMETRIC ANALYSIS OF Southern cassowary CASQUE OUTLINES######
################LANDMARK-BOOKSTEIN ALIGNMENT with SINGLE ORIGIN#################
#############################(c)2021 David Ian Kay##############################
############################written in R3.6.3 64-bit############################
################################################################################
#command for clearing out all read in data, objects, attached data, etc.
rm(list=ls())

#PACKAGES
#Read in Momocs, a morphometrics package. We am using v1.3.0
library(Momocs)
#Read in tibble to coerce dataframes to tibble type for momocs to use v3.0.3
library(tibble)
#Read in ggplot2 library, v3.3.2
library(ggplot2)
#Packages necessary for the linear discriminant analysis (LDA) with the principal coordinate analysis (PCO) results
library(MASS) #v7.3-51.5
library(vegan) #v2.5-6
#Command to display the loaded packages, versions, and other attached information
sessionInfo()

################################################################################
################################################################################
########################RIGHT/LEFT SPECIFICITY ANALYSIS#########################
################################################################################
################################################################################

##############################Initial Data read-in##############################
#Initial data read-in of the Out object with landmarks assigned created in the data_landmarking_.R script
load("FILE_PATH/cc_test_ldk_bk_slid.RData")
#Check the objects from the loaded data file
ls()

#######################Elliptical Fourier Transformation########################
#Calibrate harmonics needed to capture 99.9% shape variance (or N-1 harmonics)
cc_test_ldk_bk_slid_harm<-calibrate_harmonicpower_efourier(cc_test_ldk_bk_slid, nb.h=20, plot=T)
cc_test_ldk_bk_slid_harm ##16 harmonics capture 99.9% shape variance
#Altered graph to show more than the default 10 harmonics
cc_test_ldk_bk_slid_harm$gg+theme_minimal()+
coord_cartesian(xlim=c(0.5,17),ylim=c(0,100))+
ggtitle('left/right landmark bk slid test Harmonic calibration')
#Calibrate the chosen number of reconstructions
calibrate_reconstructions_efourier(cc_test_ldk_bk_slid,range=1:16)
#Elliptical Fourier transformation with the number of harmonics previously chosen from the calibration
cc_test_ldk_bk_slid_e_trans<-efourier(cc_test_ldk_bk_slid,nb.h=16,norm=F)

#########################Principal coordinate analysis##########################
#Using the MDS function in Momocs
cc_test_ldk_bk_slid_MDS_MOMOCS<-MDS(cc_test_ldk_bk_slid_e_trans,k=39)
#Examine the output for the MDS
str(cc_test_ldk_bk_slid_MDS_MOMOCS)
summary(cc_test_ldk_bk_slid_MDS_MOMOCS)
cc_test_ldk_bk_slid_MDS_MOMOCS

#Plot out results (first two principal coordinates only)
plot_MDS(cc_test_ldk_bk_slid_MDS_MOMOCS)
plot_MDS(cc_test_ldk_bk_slid_MDS_MOMOCS, ~side)
plot_MDS(cc_test_ldk_bk_slid_MDS_MOMOCS, ~side,chullfilled=T)

#using the capland function in vegan package
#to set up a PCO, the "formula" portion of the function call is #the data~1
#We chose euclidean distance to construct the dissimilarity #matrix
cc_test_ldk_bk_slid_PCO_MOMOCS<-capscale(cc_test_ldk_bk_slid_e_trans~1,distance="euclidean")
#look at output to determine number of axes to explain 99% of variance
summary(cc_test_ldk_bk_slid_PCO_MOMOCS)

cc_test_ldk_bk_slid_slid_MDS_MOMOCS_scores<-scores(cc_test_ldk_bk_slid_slid_MDS_MOMOCS, choices=c(1,2,3,4,5,6,7,8),display=c('sites'))
#Write results to csv to be used in the MANOVA and LDA and to have a separate copy
write.csv(cc_test_ldk_bk_slid_PCO_MOMOCS_scores,file='FILE_PATH/cc_test_ldk_bk_slid_PCO_MOMOCS.csv')

#read in the data for the MANOVA and LDA
cc_test_ldk_bk_slid_PCO_MOMOCS_dat<-read.csv('FILE_PATH/cc_test_ldk_bk_slid_PCO_MOMOCS.csv',header=T,row.names=1)
cc_test_ldk_bk_slid_PCO_MOMOCS_dat

#read in the factor data for sides
cass_info<-read.csv("FILE_PATH/test_info.csv",header=T,row.names=1)

#create a dataframe of the principal coordinates and the side factor data
cc_test_ldk_bk_slid_PCO_MOMOCS_dat_fac<-cbind(cc_test_ldk_bk_slid_PCO_MOMOCS_dat,cass_info$fac.side)
cc_test_ldk_bk_slid_PCO_MOMOCS_dat_fac

#####################################MANOVA#####################################
#Run a multiple analysis of variance on the principal coordinates
cc_test_ldk_bk_slid_PCO_MANOVA<-aov(MDS1+MDS2+MDS3+MDS4+MDS5+MDS6+MDS7+MDS8~fac.side, data=cc_test_ldk_bk_slid_PCO_MOMOCS_dat_fac)
#Look at the output
summary(cc_test_ldk_bk_slid_PCO_MANOVA)
cc_test_ldk_bk_slid_PCO_MANOVA

######################################LDA#######################################
#Use the data from the 8 prinipal coordinates produced
cc_test_ldk_bk_slid_PCO_MOMOCS_lda<-lda(fac.side~MDS1+MDS2+MDS3+MDS4+MDS5+MDS6+MDS7+MDS8,data=cc_test_ldk_bk_slid_PCO_MOMOCS_dat_fac)
#Use the predict function to test the LDA, but establish the principal coordinate data as a dataframe
cc_test_ldk_bk_slid_PCO_MOMOCS_lda_predict<-predict(cc_test_ldk_bk_slid_PCO_MOMOCS_lda,newdata=as.data.frame(cc_test_ldk_bk_slid_PCO_MOMOCS_dat))
#Check the predicted portion for a % correct
cc_test_ldk_bk_slid_PCO_MOMOCS_lda_predict$class
#Build a CV table
CV.fac_all<-cc_test_ldk_bk_slid_PCO_MOMOCS_lda_predict$class
CV.tab_all<-table(cc_test_ldk_bk_slid_PCO_MOMOCS_dat[,9],CV.fac_all)
names(dimnames(CV.tab_all))<-c('actual','classified')
CV.correct_all<-sum(diag(CV.tab_all))/sum(CV.tab_all)
tab_all <- CV.tab_all
  ce_all <- sapply(seq_along(1:nrow(tab_all)),
               function(i) 1-(sum(tab_all[i, -i])/sum(tab_all[i, ])))
  names(ce_all) <- rownames(tab_all)
#Correct classification rate
ce_all
#Classification table
tab_all

######################################PCA#######################################
#Running a principal components analysis as well, mostly to confirm that the results would show a "horseshoe" shape due to the outlines being auto-correlated
cc_test_ldk_bk_slid_PCA<-PCA(cc_test_ldk_bk_slid_e_trans)
#Output of PCA
cc_test_ldk_bk_slid_PCA
summary(cc_test_ldk_bk_slid_PCA)
#Plot the PCA results
plot_PCA(cc_test_ldk_bk_slid_PCA)
plot_PCA(cc_test_ldk_bk_slid_PCA,'side')
plot_PCA(cc_test_ldk_bk_slid_PCA,'side',chullfilled=T)

#####################################MANOVA#####################################
#MANOVA of the principal coordinates
cc_test_ldk_bk_slid_MANOVA<-MANOVA(cc_test_ldk_bk_slid_PCA,'side')
cc_test_ldk_bk_slid_MANOVA
#Pairwise MANOVA
cc_test_ldk_bk_slid_MANOVA_PW<-MANOVA_PW(cc_test_ldk_bk_slid_PCA,'side')
cc_test_ldk_bk_slid_MANOVA_PW

######################################LDA#######################################
#LDA of the principal components
cc_test_ldk_bk_slid_LDA<-LDA(cc_test_ldk_bk_slid_PCA,'side')
cc_test_ldk_bk_slid_LDA

plot_CV(cc_test_ldk_bk_slid_LDA)

##############################Workspace save/load###############################
#Save workspace, commented to prevent erroneous execution
# save.image("FILE_PATH/cass_test_ldk_bk_slid.RData")
#load workspace, commented to prevent erroneous execution
# load("FILE_PATH/cass_test_ldk_bk_slid.RData")
#Check the objects in the R session
ls()


################################################################################
################################################################################
##############################ALL CASQUE ANALYSIS###############################
################################################################################
################################################################################
################################################################################
###############################CC LATERAL ASPECT################################

##############################Initial Data read-in##############################
#Initial data read-in of the Out object with landmarks assigned created in the data_landmarking_2.R script
load("FILE_PATH/cc_lat_all_ldk_bk_slid.RData")
ls()

#######################Elliptical Fourier Transformation########################
#Calibrate harmonics needed to capture 99.9% shape variance (or a maximum of N-1 harmonics)
cal_cc_lat_all_ldk_bk_slid<-calibrate_harmonicpower_efourier(cc_lat_all_ldk_bk_slid,nb.h=20,plot=T)
cal_cc_lat_all_ldk_bk_slid
#Altered graph to show more than the default 10 harmonics
cal_cc_lat_all_ldk_bk_slid$gg+theme_minimal()+
coord_cartesian(xlim=c(0.5,17),ylim=c(0,100))+
ggtitle('Southern Cassowary lateral landmark bookstein slid Harmonic calibration')
#Visualize the reonstruction at various harmonics
cal_cc_lat_all_ldk_bk_slid_recon<-calibrate_reconstructions_efourier(cc_lat_all_ldk_bk_slid,range=1:16)
cal_cc_lat_all_ldk_bk_slid_recon
#Elliptical Fourier transformation with the number of harmonics previously chosen from the calibration.
cc_lat_all_ldk_bk_slid_ef<-efourier(cc_lat_all_ldk_bk_slid,nb.h=16,norm=F)

######################################PCA#######################################
#Ordinate the data with a Principal Components Analysis
cc_lat_all_ldk_bk_slid_PCA<-PCA(cc_lat_all_ldk_bk_slid_ef)
cc_lat_all_ldk_bk_slid_PCA
summary(cc_lat_all_ldk_bk_slid_PCA)
#Plot out the results, export to an svg
plot_PCA(cc_lat_all_ldk_bk_slid_PCA)
svg(filename="FILE_PATH/cc_lat_all_ldk_bk_slid_PCA_001.svg")
plot_PCA(cc_lat_all_ldk_bk_slid_PCA)
dev.off()
#PCA plot with point labels, export to an svg
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,labelpoints=T)
svg(filename="FILE_PATH/cc_lat_all_ldk_bk_slid_PCA_002.svg")
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,labelpoints=T)
dev.off()
#Plotting sex results
#Plot with point labels, export to an svg
plot(cc_lat_all_ldk_bk_slid_PCA,'sex',ellipses=T)
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'sex',labelpoints=T)
svg(filename='FILE_PATH/cc_lat_all_ldk_bk_slid_PCA_003.svg')
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'sex',labelpoints=T)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'sex')
svg(filename='FILE_PATH/cc_lat_all_ldk_bk_slid_PCA_004.svg')
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'sex')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'sex',chullfilled=T)
svg(filename='FILE_PATH/cc_lat_all_ldk_bk_slid_PCA_005.svg')
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'sex',chullfilled=T)
dev.off()
#Plotting geography results
#Plot with point labels, export to an svg
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'geo',labelpoints=T)
svg(filename='FILE_PATH/cc_lat_all_ldk_bk_slid_PCA_006.svg')
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'geo',labelpoints=T)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'geo')
svg(filename='FILE_PATH/cc_lat_all_ldk_bk_slid_PCA_007.svg')
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'geo')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'geo',chullfilled=T)
svg(filename='FILE_PATH/cc_lat_all_ldk_bk_slid_PCA_008.svg')
plot_PCA(cc_lat_all_ldk_bk_slid_PCA,'geo',chullfilled=T)
dev.off()

#############################Outlier Identification#############################
#Identify potential outliers using the "Which_out" function
which_out(cc_lat_all_ldk_bk_slid_PCA)

##############################Workspace save/load###############################
#Save workspace, commented to prevent erroneous execution
# save.image("FILE_PATH/cass_lat_all_ldk_bk_slid.RData")
#load workspace, commented to prevent erroneous execution
# load("FILE_PATH/cass_lat_sex_ldk_bk_slid.RData")
#Check the objects in the R session
ls()

################################################################################
###############################CC ANTERIOR ASPECT###############################

##############################Initial Data read-in##############################
#Initial data read-in of the Out object with landmarks assigned created in the data_landmarking_2.R script
load("FILE_PATH/cc_ant_all_ldk_bk_slid.RData")
#Check the objects from the loaded data file
ls()

#######################Elliptical Fourier Transformation########################
#Calibrate harmonics needed to capture 99.9% shape variance (or a maximum of N-1 harmonics)
cal_cc_ant_all_ldk_bk_slid<-calibrate_harmonicpower_efourier(cc_ant_all_ldk_bk_slid,nb.h=20,plot=T)
cal_cc_ant_all_ldk_bk_slid ##14 harmonics capture 99.9% of shape
#Altered graph to show more than the default 10 harmonics
cal_cc_ant_all_ldk_bk_slid$gg+theme_minimal()+
coord_cartesian(xlim=c(0.5,15),ylim=c(0,100))+
ggtitle('Southern Cassowary anterior all landmark bookstein-slid Harmonic calibration')
#Calibrate the chosen number of reconstructions
cal_cc_ant_all_ldk_bk_slid_recon<-calibrate_reconstructions_efourier(cc_ant_all_ldk_bk_slid,range=1:14)
cal_cc_ant_all_ldk_bk_slid_recon
#Elliptical Fourier transformation with the number of harmonics previously chosen from the calibration
cc_ant_all_ldk_bk_slid_ef<-efourier(cc_ant_all_ldk_bk_slid,nb.h=14,norm=F)

######################################PCA#######################################
##Ordinate the data with a Principal Components Analysis
cc_ant_all_ldk_bk_slid_PCA<-PCA(cc_ant_all_ldk_bk_slid_ef)
#PCA output
cc_ant_all_ldk_bk_slid_PCA
summary(cc_ant_all_ldk_bk_slid_PCA)
#Plot out the results, export to an svg
plot_PCA(cc_ant_all_ldk_bk_slid_PCA)
svg(filename="FILE_PATH/cc_ant_all_ldk_bk_slid_PCA_001.svg")
plot_PCA(cc_ant_all_ldk_bk_slid_PCA)
dev.off()
#PCA plot with point labels, export to an svg
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,labelpoints=T)
svg(filename="FILE_PATH/cc_ant_all_ldk_bk_slid_PCA_002.svg")
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,labelpoints=T)
dev.off()
#Plotting sex results
#Plot with point labels, export to an svg
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'sex',labelpoints=T)
svg(filename='FILE_PATH/cc_ant_all_ldk_bk_slid_PCA_003.svg')
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'sex',labelpoints=T)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'sex')
svg(filename='FILE_PATH/cc_ant_all_ldk_bk_slid_PCA_004.svg')
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'sex')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'sex',chullfilled=T)
svg(filename='FILE_PATH/cc_ant_all_ldk_bk_slid_PCA_005.svg')
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'sex',chullfilled=T)
dev.off()

#Plotting geography results
##Plot with point labels, export to an svg
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'geo',labelpoints=T)
svg(filename='FILE_PATH/cc_ant_all_ldk_bk_slid_PCA_006.svg')
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'geo',labelpoints=T)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'geo')
svg(filename='FILE_PATH/cc_ant_all_ldk_bk_slid_PCA_007.svg')
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'geo')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'geo',chullfilled=T)
svg(filename='FILE_PATH/cc_ant_all_ldk_bk_slid_PCA_008.svg')
plot_PCA(cc_ant_all_ldk_bk_slid_PCA,'geo',chullfilled=T)
dev.off()

#############################Outlier Identification#############################
#Identify potential outliers using the "Which_out" function
which_out(cc_ant_all_ldk_bk_slid_PCA)

##############################Workspace save/load###############################
#save the workspace for future access
# save.image("FILE_PATH/cass_ant_all_ldk_bk_slid.RData")
#command to load the workspace
# load("FILE_PATH/cass_ant_all_ldk_bk_slid.RData")
#Check the objects in the R session
ls()


################################################################################
################################################################################
#############################SEX CASQUE ANALYSIS################################
################################################################################
################################################################################
################################################################################
#############################CC SEX LATERAL ASPECT##############################

##############################Initial Data read-in##############################
#Initial data read-in of the Out object with landmarks assigned created in the data_landmarking_2.R script
# load("FILE_PATH/cc_lat_sex_ldk_bk_slid.RData")
ls()
stack(cc_lat_sex_ldk_bk_slid)

#######################Elliptical Fourier Transformation########################
#Calibrate harmonics needed to capture 99.9% shape variance (or N-1 harmonics)
cal_cc_lat_sex_ldk_bk_slid<-calibrate_harmonicpower_efourier(cc_lat_sex_ldk_bk_slid,nb.h=20,plot=T)
cal_cc_lat_sex_ldk_bk_slid
#Altered graph to show more than the default 10 harmonics
cal_cc_lat_sex_ldk_bk_slid$gg+theme_minimal()+
coord_cartesian(xlim=c(0.5,17),ylim=c(0,100))+
ggtitle('Southern Cassowary known-sex lateral landmark bookstein slid Harmonic calibration')
#Calibrate the chosen number of reconstructions
cal_cc_lat_sex_ldk_bk_slid_recon<-calibrate_reconstructions_efourier(cc_lat_sex_ldk_bk_slid,range=1:16)
cal_cc_lat_sex_ldk_bk_slid_recon
#Elliptical Fourier transformation with the number of harmonics previously chosen from the calibration
cc_lat_sex_ldk_bk_slid_ef<-efourier(cc_lat_sex_ldk_bk_slid,nb.h=16,norm=F)

######################################PCA#######################################
##Ordinate the data with a Principal Components Analysis
cc_lat_sex_ldk_bk_slid_PCA<-PCA(cc_lat_sex_ldk_bk_slid_ef,fac='sex')
cc_lat_sex_ldk_bk_slid_PCA
summary(cc_lat_sex_ldk_bk_slid_PCA)
#Plot out the results
plot_PCA(cc_lat_sex_ldk_bk_slid_PCA)
svg(filename="FILE_PATH/cc_lat_sex_ldk_bk_slid_PCA_001.svg")
plot_PCA(cc_lat_sex_ldk_bk_slid_PCA)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_lat_sex_ldk_bk_slid_PCA,'sex')
svg(filename='FILE_PATH/cc_lat_sex_ldk_bk_slid_PCA_002.svg')
plot_PCA(cc_lat_sex_ldk_bk_slid_PCA,'sex')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_lat_sex_ldk_bk_slid_PCA,'sex',chullfilled=T)
svg(filename='FILE_PATH/cc_lat_sex_ldk_bk_slid_PCA_003.svg')
plot_PCA(cc_lat_sex_ldk_bk_slid_PCA,'sex',chullfilled=T)
dev.off()

######################################LDA#######################################
#Run an LDA on the PCA results, retaining the number of axes capturing 99% variance or supported by sample size
cc_lat_sex_ldk_bk_slid_LDA<-LDA(cc_lat_sex_ldk_bk_slid_PCA,'sex',retain=9)
#LDA output
cc_lat_sex_ldk_bk_slid_LDA
#Create an LDA plot, and an LDA cross-validation plot
plot_LDA(cc_lat_sex_ldk_bk_slid_LDA,'sex')
plot_CV(cc_lat_sex_ldk_bk_slid_LDA)

#####################################MANOVA#####################################
#MANOVA of the PCA ordinated results, retaining the number of axes capturing 99% variance or supported by sample size
cc_lat_sex_ldk_bk_slid_MANOVA<-MANOVA(cc_lat_sex_ldk_bk_slid_PCA,fac='sex',retain=9)
#MANOVA output
cc_lat_sex_ldk_bk_slid_MANOVA
#Pairwise MANOVA
cc_lat_sex_ldk_bk_slid_PW_MANOVA<-MANOVA_PW(cc_lat_sex_ldk_bk_slid_PCA,fac='sex')
#Pairwise MANOVA output
cc_lat_sex_ldk_bk_slid_PW_MANOVA

##############################Workspace save/load###############################
#save the workspace for future access, commented to prevent erroneous execution
# save.image("FILE_PATH/cass_lat_sex_ldk_bk_slid.RData")
#command to load the workspace, commented to prevent erroneous execution
# load("FILE_PATH/cass_lat_sex_ldk_bk_slid.RData")


################################################################################
#############################CC SEX ANTERIOR ASPECT#############################

##############################Initial Data read-in##############################
#Initial data read-in of the Out object with landmarks assigned created in the data_landmarking_2.R script
# load("FILE_PATH/cc_ant_sex_ldk_bk_slid.RData")
ls()

#######################Elliptical Fourier Transformation########################
#Calibrate harmonics needed to capture 99.9% shape variance (or N-1 harmonics)
cal_cc_ant_ldk_bk_slid<-calibrate_harmonicpower_efourier(cc_ant_sex_ldk_bk_slid,nb.h=20,plot=T)
cal_cc_ant_ldk_bk_slid ##14 harmonicscapture 99.9% of the shape
#Altered graph to show more than the default 10 harmonics
cal_cc_ant_ldk_bk_slid$gg+theme_minimal()+
coord_cartesian(xlim=c(0.5,15),ylim=c(0,100))+
ggtitle('Southern Cassowary known sex anterior landmark bookstein-slid Harmonic calibration')
#Calibrate the chosen number of reconstructions
cal_cc_ant_sex_ldk_bk_slid_recon<-calibrate_reconstructions_efourier(cc_ant_sex_ldk_bk_slid,range=1:14)
cal_cc_ant_sex_ldk_bk_slid_recon
#Elliptical Fourier transformation with the number of harmonics previously chosen from the calibration
cc_ant_sex_ldk_bk_slid_ef<-efourier(cc_ant_sex_ldk_bk_slid,nb.h=14,norm=F)

######################################PCA#######################################
##Ordinate the data with a Principal Components Analysis
cc_ant_sex_ldk_bk_slid_PCA<-PCA(cc_ant_sex_ldk_bk_slid_ef,fac='sex')
#PCA output
cc_ant_sex_ldk_bk_slid_PCA
summary(cc_ant_sex_ldk_bk_slid_PCA)
#Plot out the results, export to an svg
plot_PCA(cc_ant_sex_ldk_bk_slid_PCA)
svg(filename="FILE_PATH/cc_ant_sex_ldk_bk_slid_PCA_001.svg")
plot_PCA(cc_ant_sex_ldk_bk_slid_PCA)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_ant_sex_ldk_bk_slid_PCA,'sex')
svg(filename='FILE_PATH/cc_ant_sex_ldk_bk_slid_PCA_002.svg')
plot_PCA(cc_ant_sex_ldk_bk_slid_PCA,'sex')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_ant_sex_ldk_bk_slid_PCA,'sex',chullfilled=T)
svg(filename='FILE_PATH/cc_ant_sex_ldk_bk_slid_PCA_003.svg')
plot_PCA(cc_ant_sex_ldk_bk_slid_PCA,'sex',chullfilled=T)
dev.off()

######################################LDA#######################################
#Run an LDA on the PCA results, retaining the number of axes capturing 99% variance or supported by sample size
cc_ant_sex_ldk_bk_slid_LDA<-LDA(cc_ant_sex_ldk_bk_slid_PCA,'sex',retain=9)
cc_ant_sex_ldk_bk_slid_LDA
#Create an LDA plot, and an LDA cross-validation plot
plot_LDA(cc_ant_sex_ldk_bk_slid_LDA,'sex')
plot_CV(cc_ant_sex_ldk_bk_slid_LDA)

#####################################MANOVA#####################################
#MANOVA of the PCA ordinated results, retaining the number of axes capturing 99% variance or supported by sample size
cc_ant_sex_ldk_bk_slid_MANOVA<-MANOVA(cc_ant_sex_ldk_bk_slid_PCA,fac='sex',retain=9)
#MANOVA output
cc_ant_sex_ldk_bk_slid_MANOVA
#Pairwise MANOVA
cc_ant_sex_ldk_bk_slid_PW_MANOVA<-MANOVA_PW(cc_ant_sex_ldk_bk_slid_PCA,fac='sex')
#Pairwise MANOVA output
cc_ant_sex_ldk_bk_slid_PW_MANOVA

##############################Workspace save/load###############################
#save the workspace for future access, commented to prevent erroneous execution
# save.image("FILE_PATH/cass_ant_sex_ldk_bk_slid.RData")
#command to load the workspace, commented to prevent erroneous execution
# load("FILE_PATH/cass_ant_sex_ldk_bk_slid.RData")


################################################################################
################################################################################
##########################GEOGRAPHY CASQUE ANALYSIS#############################
################################################################################
################################################################################
################################################################################
##########################CC GEOGRAPHY LATERAL ASPECT###########################

##############################Initial Data read-in##############################
#Initial data read-in of the Out object with landmarks assigned created in the data_landmarking_2.R script
# load("FILE_PATH/cc_lat_geo_ldk_bk_slid.RData")
ls()

#######################Elliptical Fourier Transformation########################
#Calibrate harmonics needed to capture 99.9% shape variance (or N-1 harmonics)
cal_cc_lat_geo_ldk_bk_slid<-calibrate_harmonicpower_efourier(cc_lat_geo_ldk_bk_slid,nb.h=20,plot=T)
cal_cc_lat_geo_ldk_bk_slid ##15 harmonics capture 99.9% of the shape
#Altered graph to show more than the default 10 harmonics
cal_cc_lat_geo_ldk_bk_slid$gg+theme_minimal()+
coord_cartesian(xlim=c(0.5,16),ylim=c(0,100))+
ggtitle('Southern Cassowary lateral aspect landmark bookstein slid geography Harmonic calibration')
#Calibrate the chosen number of reconstructions
cal_cc_lat_geo_ldk_bk_slid_recon<-calibrate_reconstructions_efourier(cc_lat_geo_ldk_bk_slid,range=1:15)
cal_cc_lat_geo_ldk_bk_slid_recon
#Elliptical Fourier transformation with the number of harmonics previously chosen from the calibration
cc_lat_geo_ldk_bk_slid_ef<-efourier(cc_lat_geo_ldk_bk_slid,nb.h=15,norm=F)

######################################PCA#######################################
##Ordinate the data with a Principal Components Analysis
cc_lat_geo_ldk_bk_slid_PCA<-PCA(cc_lat_geo_ldk_bk_slid_ef,fac='geo')
#PCA output
cc_lat_geo_ldk_bk_slid_PCA
summary(cc_lat_geo_ldk_bk_slid_PCA)
#Plot out the results
plot_PCA(cc_lat_geo_ldk_bk_slid_PCA)
svg(filename="FILE_PATH/cc_lat_geo_ldk_bk_slid_PCA_001.svg")
plot_PCA(cc_lat_geo_ldk_bk_slid_PCA)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_lat_geo_ldk_bk_slid_PCA,'geo')
svg(filename='FILE_PATH/cc_lat_geo_ldk_bk_slid_PCA_002.svg')
plot_PCA(cc_lat_geo_ldk_bk_slid_PCA,'geo')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_lat_geo_ldk_bk_slid_PCA,'geo',chullfilled=T)
svg(filename='FILE_PATH/cc_lat_geo_ldk_bk_slid_PCA_003.svg')
plot_PCA(cc_lat_geo_ldk_bk_slid_PCA,'geo',chullfilled=T)
dev.off()

######################################LDA#######################################
#Run an LDA on the PCA results, retaining the number of axes capturing 99% variance or supported by sample size
cc_lat_geo_ldk_bk_slid_LDA<-LDA(cc_lat_geo_ldk_bk_slid_PCA,'geo',retain=8)
#LDA output
cc_lat_geo_ldk_bk_slid_LDA
#Create an LDA plot, and an LDA cross-validation plot
plot_LDA(cc_lat_geo_ldk_bk_slid_LDA,'geo')
plot_CV(cc_lat_geo_ldk_bk_slid_LDA)

#####################################MANOVA#####################################
#MANOVA of the PCA ordinated results, retaining the number of axes capturing 99% variance or supported by sample size
cc_lat_geo_ldk_bk_slid_MANOVA<-MANOVA(cc_lat_geo_ldk_bk_slid_PCA,fac='geo',retain=8)
#MANOVA output
cc_lat_geo_ldk_bk_slid_MANOVA
#Pairwise MANOVA
cc_lat_geo_ldk_bk_slid_PW_MANOVA<-MANOVA_PW(cc_lat_geo_ldk_bk_slid_PCA,fac='geo')
#Pairwise MANOVA output
cc_lat_geo_ldk_bk_slid_PW_MANOVA

##############################Workspace save/load###############################
#save the workspace for future access, commented to prevent erroneous execution
# save.image("FILE_PATH/bookstein_slid/cass_lat_geo_ldk_bk_slid.RData")
#command to load the workspace, commented to prevent erroneous execution
# load("FILE_PATH/cass_lat_geo_ldk_bk_slid.RData")


################################################################################
##########################CC GEOGRAPHY ANTERIOR ASPECT##########################

##############################Initial Data read-in##############################
#Initial data read-in of the Out object with landmarks assigned created in the data_landmarking_2.R script
# load("FILE_PATH/cc_ant_geo_ldk_bk_slid.RData")
ls()
stack(cc_ant_geo_ldk_bk_slid)

#######################Elliptical Fourier Transformation########################
#Calibrate harmonics needed to capture 99.9% shape variance (or N-1 harmonics)
cal_cc_ant_ldk_bk_slid<-calibrate_harmonicpower_efourier(cc_ant_geo_ldk_bk_slid,nb.h=20,plot=T)
cal_cc_ant_ldk_bk_slid
#Altered graph to show more than the default 10 harmonics
cal_cc_ant_ldk_bk_slid$gg+theme_minimal()+
coord_cartesian(xlim=c(0.5,15),ylim=c(0,100))+
ggtitle('Southern Cassowary known geo anterior landmark bookstein Harmonic calibration')

##Typically choose N-1 harmonics
#calibrate the deviations.
#calibrate_deviations_efourier(cass_pro)
#calibrate the reconstructions, as a visualization for
cal_cc_ant_geo_ldk_bk_slid_recon<-calibrate_reconstructions_efourier(cc_ant_geo_ldk_bk_slid,range=1:14)
cal_cc_ant_geo_ldk_bk_slid_recon
#elliptical Fourier transformation with the number of harmonics previously
#chosen from the calibration. I do not normalize the coefficients because
#of the poor alignments produced by
cc_ant_geo_ldk_bk_slid_ef<-efourier(cc_ant_geo_ldk_bk_slid,nb.h=14,norm=F)
#need to access harmonics data to then send to the principal coordinate
#analysis
#can access harmonic data by calling the $coe part of the Fourier transformed
#data object
cc_ant_geo_ldk_bk_slid_ef
cc_ant_geo_ldk_bk_slid_ef$coe
cc_ant_geo_ldk_bk_slid_ef$fac

######################################PCA#######################################
##Ordinate the data with a Principal Components Analysis
cc_ant_geo_ldk_bk_slid_PCA<-PCA(cc_ant_geo_ldk_bk_slid_ef)
cc_ant_geo_ldk_bk_slid_PCA
summary(cc_ant_geo_ldk_bk_slid_PCA)
#plot out the results, export to an svg
plot_PCA(cc_ant_geo_ldk_bk_slid_PCA)
svg(filename="FILE_PATH/cc_ant_geo_ldk_bk_slid_PCA_001.svg")
plot_PCA(cc_ant_geo_ldk_bk_slid_PCA)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_ant_geo_ldk_bk_slid_PCA,'geo')
svg(filename='FILE_PATH/cc_ant_geo_ldk_bk_slid_PCA_002.svg')
plot_PCA(cc_ant_geo_ldk_bk_slid_PCA,'geo')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_ant_geo_ldk_bk_slid_PCA,'geo',chullfilled=T)
svg(filename='FILE_PATH/cc_ant_geo_ldk_bk_slid_PCA_003.svg')
plot_PCA(cc_ant_geo_ldk_bk_slid_PCA,'geo',chullfilled=T)
dev.off()

######################################LDA#######################################
cc_ant_geo_ldk_bk_slid_LDA<-LDA(cc_ant_geo_ldk_bk_slid_PCA,'geo',retain=6)
cc_ant_geo_ldk_bk_slid_LDA
#Create an LDA plot, and an LDA cross-validation plot
plot_LDA(cc_ant_geo_ldk_bk_slid_LDA,'geo')
plot_CV(cc_ant_geo_ldk_bk_slid_LDA)

#####################################MANOVA#####################################
cc_ant_geo_ldk_bk_slid_MANOVA<-MANOVA(cc_ant_geo_ldk_bk_slid_PCA,fac='geo',retain=6)
cc_ant_geo_ldk_bk_slid_MANOVA
cc_ant_geo_ldk_bk_slid_PW_MANOVA<-MANOVA_PW(cc_ant_geo_ldk_bk_slid_PCA,fac='geo')
cc_ant_geo_ldk_bk_slid_PW_MANOVA

##############################Workspace save/load###############################
#save the workspace for future access, commented to prevent erroneous execution
# save.image("FILE_PATH/cass_ant_geo_ldk_bk_slid.RData")
#command to load the workspace, commented to prevent erroneous execution
# load("FILE_PATH/cass_ant_geo_ldk_bk_slid.RData")


################################################################################
################################################################################
##########################GEOGRAPHY-SEX CASQUE ANALYSIS#########################
################################################################################
################################################################################
################################################################################
###############################CC LATERAL GEO-SEX###############################

##############################Initial Data read-in##############################
#Initial data read-in of the Out object with landmarks assigned created in the data_landmarking_2.R script
# load("FILE_PATH/cc_lat_geo_sex_ldk_bk_slid.RData")
ls()

#######################Elliptical Fourier Transformation########################
#Calibrate harmonics needed to capture 99.9% shape variance (or N-1 harmonics)
cal_cc_lat_geo_sex_ldk_bk_slid<-calibrate_harmonicpower_efourier(cc_lat_geo_sex_ldk_bk_slid,nb.h=20,plot=T)
cal_cc_lat_geo_sex_ldk_bk_slid ##16 harmonics capture 99.9% of the shape
#Altered graph to show more than the default 10 harmonics
cal_cc_lat_geo_sex_ldk_bk_slid$gg+theme_minimal()+
coord_cartesian(xlim=c(0.5,17),ylim=c(0,100))+
ggtitle('Southern Cassowary lateral aspect landmark bookstein slid geography-sex Harmonic calibration')
#Calibrate the chosen number of reconstructions
cal_cc_lat_geo_sex_ldk_bk_slid_recon<-calibrate_reconstructions_efourier(cc_lat_geo_sex_ldk_bk_slid,range=1:16)
#Calibrate the chosen number of reconstructions
cal_cc_lat_geo_sex_ldk_bk_slid_recon
#Elliptical Fourier transformation with the number of harmonics previously chosen from the calibration
cc_lat_geo_sex_ldk_bk_slid_ef<-efourier(cc_lat_geo_sex_ldk_bk_slid,nb.h=16,norm=F)

######################################PCA#######################################
#Ordinate the data with a Principal Components Analysis
cc_lat_geo_sex_ldk_bk_slid_PCA<-PCA(cc_lat_geo_sex_ldk_bk_slid_ef,fac=c('sex','geo'))
#PCA output
cc_lat_geo_sex_ldk_bk_slid_PCA
summary(cc_lat_geo_sex_ldk_bk_slid_PCA)
#Plot out the results
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_PCA)
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_PCA,~sex+geo)
#Plot with convex hulls
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_PCA,~sex+geo,chullfilled=T)
#plot with 95% CI ellipses
plot(cc_lat_geo_sex_ldk_bk_slid_PCA,~sex+geo,ellipses=T)

#####################################MANOVA#####################################
#MANOVA of the PCA ordinated results, retaining the number of axes capturing 99% variance or supported by sample size
cc_lat_geo_sex_ldk_bk_slid_MANOVA<-MANOVA(cc_lat_geo_sex_ldk_bk_slid_PCA,fac='sex',retain=3)
#MANOVA output
cc_lat_geo_sex_ldk_bk_slid_MANOVA
#Pairwise MANOVA
cc_lat_geo_sex_ldk_bk_slid_PW_MANOVA<-MANOVA_PW(cc_lat_geo_sex_ldk_bk_slid_PCA,~sex+geo)
#Pairwise MANOVA output
cc_lat_geo_sex_ldk_bk_slid_PW_MANOVA

########################slice the data to just the SPNG#########################
cc_lat_geo_sex_ldk_bk_slid_SPNG_ef<-slice(cc_lat_geo_sex_ldk_bk_slid_ef,c(1:6,15))
#Check factors for accuracy
cc_lat_geo_sex_ldk_bk_slid_SPNG_ef$fac

######################################PCA#######################################
##Ordinate the data with a Principal Components Analysis
cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA<-PCA(cc_lat_geo_sex_ldk_bk_slid_SPNG_ef,fac='sex')
cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA
summary(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA)
#Plot out the results, export to an svg
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA)
svg(filename="FILE_PATH/cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA_001.svg")
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA,'sex')
svg(filename='FILE_PATH/cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA_002.svg')
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA,'sex')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA,'sex',chullfilled=T)
svg(filename='FILE_PATH/cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA_003.svg')
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA,'sex',chullfilled=T)
dev.off()

######################################LDA#######################################
#Run an LDA on the PCA results, retaining the number of axes capturing 99% variance or supported by sample size
cc_lat_geo_sex_ldk_bk_slid_SPNG_LDA<-LDA(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA,'sex',retain=3)
#LDA output
cc_lat_geo_sex_ldk_bk_slid_SPNG_LDA
#Create an LDA plot, and an LDA cross-validation plot
plot_LDA(cc_lat_geo_sex_ldk_bk_slid_SPNG_LDA,'sex')
plot_CV(cc_lat_geo_sex_ldk_bk_slid_SPNG_LDA)

#####################################MANOVA#####################################
#MANOVA of the PCA ordinated results, retaining the number of axes capturing 99% variance or supported by sample size
cc_lat_geo_sex_ldk_bk_slid_SPNG_MANOVA<-MANOVA(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA,fac='sex',retain=3)
#MANOVA output
cc_lat_geo_sex_ldk_bk_slid_SPNG_MANOVA
#Pairwise MANOVA
cc_lat_geo_sex_ldk_bk_slid_SPNG_PW_MANOVA<-MANOVA_PW(cc_lat_geo_sex_ldk_bk_slid_SPNG_PCA,fac='sex')
#Pairwise MANOVA output
cc_lat_geo_sex_ldk_bk_slid_SPNG_PW_MANOVA

#########################slice the data to just the AUS#########################
cc_lat_geo_sex_ldk_bk_slid_AUS_ef<-slice(cc_lat_geo_sex_ldk_bk_slid_ef,c(7:14,16:21))
#Check factors for accuracy
cc_lat_geo_sex_ldk_bk_slid_AUS_ef$fac

######################################PCA#######################################
##Ordinate the data with a Principal Components Analysis
cc_lat_geo_sex_ldk_bk_slid_AUS_PCA<-PCA(cc_lat_geo_sex_ldk_bk_slid_AUS_ef,fac='sex')
cc_lat_geo_sex_ldk_bk_slid_AUS_PCA
summary(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA)
#Plot out the results, export to an svg
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA)
svg(filename="FILE_PATH/cc_lat_geo_sex_ldk_bk_slid_AUS_PCA_001.svg")
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA,'sex')
svg(filename='FILE_PATH/cc_lat_geo_sex_ldk_bk_slid_AUS_PCA_002.svg')
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA,'sex')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA,'sex',chullfilled=T)
svg(filename='FILE_PATH/cc_lat_geo_sex_ldk_bk_slid_AUS_PCA_003.svg')
plot_PCA(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA,'sex',chullfilled=T)
dev.off()

######################################LDA#######################################
#Run an LDA on the PCA results, retaining the number of axes capturing 99% variance or supported by sample size
cc_lat_geo_sex_ldk_bk_slid_AUS_LDA<-LDA(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA,'sex',retain=6)
#LDA output
cc_lat_geo_sex_ldk_bk_slid_AUS_LDA
#Create an LDA plot, and an LDA cross-validation plot
plot_LDA(cc_lat_geo_sex_ldk_bk_slid_AUS_LDA,'sex')
plot_CV(cc_lat_geo_sex_ldk_bk_slid_AUS_LDA)

#####################################MANOVA#####################################
#MANOVA of the PCA ordinated results, retaining the number of axes capturing 99% variance or supported by sample size
cc_lat_geo_sex_ldk_bk_slid_AUS_MANOVA<-MANOVA(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA,fac='sex',retain=6)
#MANOVA output
cc_lat_geo_sex_ldk_bk_slid_AUS_MANOVA
#Pairwise MANOVA
cc_lat_geo_sex_ldk_bk_slid_AUS_PW_MANOVA<-MANOVA_PW(cc_lat_geo_sex_ldk_bk_slid_AUS_PCA,fac='sex')
#Pairwise MANOVA output
cc_lat_geo_sex_ldk_bk_slid_AUS_PW_MANOVA


##############################Workspace save/load###############################
#save the workspace for future access, commented to prevent erroneous execution
# save.image("FILE_PATH/cass_lat_geo_sex_ldk_bk_slid.RData")
#command to load the workspace, commented to prevent erroneous execution
# load("FILE_PATH/cass_lat_geo_sex_ldk_bk_slid.RData")


################################################################################
##############################CC ANTERIOR GEO-SEX###############################

##############################Initial Data read-in##############################
#Initial data read-in of the Out object with landmarks assigned created in the data_landmarking_2.R script
# load("FILE_PATH/cc_ant_geo_sex_ldk_bk_slid.RData")
ls()
stack(cc_ant_geo_sex_ldk_bk_slid)

#######################Elliptical Fourier Transformation########################
#Calibrate harmonics needed to capture 99.9% shape variance (or N-1 harmonics)
cal_cc_ant_geo_sex_ldk_bk_slid<-calibrate_harmonicpower_efourier(cc_ant_geo_sex_ldk_bk_slid,nb.h=20,plot=T)
cal_cc_ant_geo_sex_ldk_bk_slid ##14 harmonics capture 99.9% of the shape
#Altered graph to show more than the default 10 harmonics
cal_cc_ant_geo_sex_ldk_bk_slid$gg+theme_minimal()+
coord_cartesian(xlim=c(0.5,15),ylim=c(0,100))+
ggtitle('Southern Cassowary anterior aspect landmark bookstein geography-sex Harmonic calibration')
#Calibrate the chosen number of reconstructions
cal_cc_ant_geo_sex_ldk_bk_slid_recon<-calibrate_reconstructions_efourier(cc_ant_geo_sex_ldk_bk_slid,range=1:14)
cal_cc_ant_geo_sex_ldk_bk_slid_recon
#Elliptical Fourier transformation with the number of harmonics previously chosen from the calibration
cc_ant_geo_sex_ldk_bk_slid_ef<-efourier(cc_ant_geo_sex_ldk_bk_slid,nb.h=14,norm=F)

######################################PCA#######################################
##Run a Principal Component Analysis to ordinate the data
cc_ant_geo_sex_ldk_bk_slid_PCA<-PCA(cc_ant_geo_sex_ldk_bk_slid_ef)
#PCA output
cc_ant_geo_sex_ldk_bk_slid_PCA
summary(cc_ant_geo_sex_ldk_bk_slid_PCA)
#Plot out the results, export to an svg
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_PCA)
svg(filename="FILE_PATH/cc_ant_geo_sex_ldk_bk_slid_PCA_001.svg")
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_PCA)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_PCA,'sex')
svg(filename='FILE_PATH/bookstein_slid/cc_ant_geo_sex_ldk_bk_slid_PCA_002.svg')
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_PCA,'sex')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_PCA,'sex',chullfilled=T)
svg(filename='FILE_PATH/cc_ant_geo_sex_ldk_bk_slid_PCA_003.svg')
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_PCA,'sex',chullfilled=T)
dev.off()

#####################################MANOVA#####################################
#MANOVA of the PCA ordinated results, retaining the number of axes capturing 99% variance or supported by sample size
cc_ant_geo_sex_ldk_bk_slid_MANOVA<-MANOVA(cc_ant_geo_sex_ldk_bk_slid_PCA,fac='sex',retain=3)
#MANOVA output
cc_ant_geo_sex_ldk_bk_slid_MANOVA
#Pairwise MANOVA
cc_ant_geo_sex_ldk_bk_slid_PW_MANOVA<-MANOVA_PW(cc_ant_geo_sex_ldk_bk_slid_PCA,~sex+geo)
#Pairwise MANOVA output
cc_ant_geo_sex_ldk_bk_slid_PW_MANOVA

########################slice the data to just the SPNG#########################
cc_ant_geo_sex_ldk_bk_slid_SPNG_ef<-slice(cc_ant_geo_sex_ldk_bk_slid_ef,c(1:6,15))
#Check factors for accuracy
cc_ant_geo_sex_ldk_bk_slid_SPNG_ef$fac

######################################PCA#######################################
##Ordinate the data with a Principal Components Analysis
cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA<-PCA(cc_ant_geo_sex_ldk_bk_slid_SPNG_ef)
cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA
summary(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA)
#plot out the results, export to an svg
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA)
svg(filename="FILE_PATH/cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA_001.svg")
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA,'sex')
svg(filename='FILE_PATH/cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA_002.svg')
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA,'sex')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA,'sex',chullfilled=T)
svg(filename='FILE_PATH/cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA_003.svg')
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA,'sex',chullfilled=T)
dev.off()

######################################LDA#######################################
#Run an LDA on the PCA results, retaining the number of axes capturing 99% variance or supported by sample size
cc_ant_geo_sex_ldk_bk_slid_SPNG_LDA<-LDA(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA,'sex',retain=3)
#LDA output
cc_ant_geo_sex_ldk_bk_slid_SPNG_LDA
#Create an LDA plot, and an LDA cross-validation plot
plot_LDA(cc_ant_geo_sex_ldk_bk_slid_SPNG_LDA,'sex')
plot_CV(cc_ant_geo_sex_ldk_bk_slid_SPNG_LDA)

#####################################MANOVA#####################################
cc_ant_geo_sex_ldk_bk_slid_SPNG_MANOVA<-MANOVA(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA,fac='sex',retain=3)
cc_ant_geo_sex_ldk_bk_slid_SPNG_MANOVA
cc_ant_geo_sex_ldk_bk_slid_SPNG_PW_MANOVA<-MANOVA_PW(cc_ant_geo_sex_ldk_bk_slid_SPNG_PCA,fac='sex')
cc_ant_geo_sex_ldk_bk_slid_SPNG_PW_MANOVA

#########################slice the data to just the AUS#########################
cc_ant_geo_sex_ldk_bk_slid_AUS_ef<-slice(cc_ant_geo_sex_ldk_bk_slid_ef,c(7:14,16))
#Check factors for accuracy
cc_ant_geo_sex_ldk_bk_slid_AUS_ef$fac

######################################PCA#######################################
##Ordinate the data with a Principal Components Analysis
cc_ant_geo_sex_ldk_bk_slid_AUS_PCA<-PCA(cc_ant_geo_sex_ldk_bk_slid_AUS_ef)
cc_ant_geo_sex_ldk_bk_slid_AUS_PCA
summary(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA)
#plot out the results, export to an svg
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA)
svg(filename="FILE_PATH/cc_ant_geo_sex_ldk_bk_slid_AUS_PCA_001.svg")
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA)
dev.off()
#Plot with convex hulls, export to an svg
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA,'sex')
svg(filename='FILE_PATH/cc_ant_geo_sex_ldk_bk_slid_AUS_PCA_002.svg')
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA,'sex')
dev.off()
#Plot with convex hulls filled in, export to an svg
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA,'sex',chullfilled=T)
svg(filename='FILE_PATH/cc_ant_geo_sex_ldk_bk_slid_AUS_PCA_003.svg')
plot_PCA(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA,'sex',chullfilled=T)
dev.off()

######################################LDA#######################################
#Run an LDA on the PCA results, retaining the number of axes capturing 99% variance or supported by sample size
cc_ant_geo_sex_ldk_bk_slid_AUS_LDA<-LDA(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA,'sex',retain=4)
#LDA output
cc_ant_geo_sex_ldk_bk_slid_AUS_LDA
#Create an LDA plot, and an LDA cross-validation plot
plot_LDA(cc_ant_geo_sex_ldk_bk_slid_AUS_LDA,'sex')
plot_CV(cc_ant_geo_sex_ldk_bk_slid_AUS_LDA)

#####################################MANOVA#####################################
#MANOVA of the PCA ordinated results, retaining the number of axes capturing 99% variance or supported by sample size
cc_ant_geo_sex_ldk_bk_slid_AUS_MANOVA<-MANOVA(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA,fac='sex',retain=4)
#MANOVA output
cc_ant_geo_sex_ldk_bk_slid_AUS_MANOVA
#Pairwise MANOVA
cc_ant_geo_sex_ldk_bk_slid_AUS_PW_MANOVA<-MANOVA_PW(cc_ant_geo_sex_ldk_bk_slid_AUS_PCA,fac='sex')
#Pairwise MANOVA output
cc_ant_geo_sex_ldk_bk_slid_AUS_PW_MANOVA

##############################Workspace save/load###############################
#save the workspace for future access, commented to prevent erroneous execution
# save.image("FILE_PATH/cass_ant_geo_sex_ldk_bk_slid.RData")
#command to load the workspace, commented to prevent erroneous execution
# load("FILE_PATH/cass_ant_geo_sex_ldk_bk_slid.RData")
