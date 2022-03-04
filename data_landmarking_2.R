################################################################################
##############################DATA LANDMARKING##################################
#############################(c)2021 David Ian Kay##############################
############################written in R3.6.3 64-bit############################
################################################################################
#PACKAGES
#Load Momocs, a morphometrics package. I am using v1.3.0
library(Momocs)
#Load tibble, a package enabling tibble data structures instead of dataframes for momocs to use v3.0.3
library(tibble)
#Command to check which versions of packages are loaded and attached in the R session
sessionInfo()
#Command to ensure that no extranneous objects exist
rm(list=ls())

################################################################################
##############################LATERAL LANDMARKING###############################
################################################################################

######################Data read-in and outline creation#########################
#Read in the factors list to categorize specimens, initialize as a tibble
cass_info<-as_tibble(read.csv("FILE_PATH/all_cc_lat_fac_ldk.csv", header=T, row.names=1))
#Check the factors are correct
cass_info
#Create an object containing the filenames to reference them for import
lf_lat <- list.files("FOLDER_PATH/caque_pics_oriented_ldk/Cc_LAT", full.names=TRUE)
#Double-check the filenames in the object
lf_lat
#Import the binary mask jpegs using the file list
cass_lat_all<-import_jpg(lf_lat)
#Convert them to outline, simultaneously adding factors to the outlines
cass_lat_all_out<-Out(cass_lat_all, fac=cass_info)
#Interpolate to the average number of points stated in the object information
cass_lat_all_out
cass_lat_all_out_int<-coo_interpolate(cass_lat_all_out,n=5305)
#Double check the point numbers are correct
cass_lat_all_out_int

###################Landmark Assignment and Outline Processing###################
##Make sure to save the object returned by the command
##There does not seem to be a way to go back through the specimens or edit points.
##Landmarks must be assigned in the same order for each specimen (I did anterior first, posterior second)
cass_lat_all_ldk<-def_ldk(cass_lat_all_out_int, 2)
#Inspect to make sure all specimens have 2 landmarks assigned
cass_lat_all_ldk
cass_lat_all_ldk$ldk
#Plot the outlines in their unaligned state, landmarks are red
stack(cc_lat_all_ldk)

#There are only two landmarks on the outlines at the ventral ends of the casques, so a Procrustes alignment cannot work. Bookstein coordinate alignment can, as it works on only two landmarks. These will be used to establish a new baseline calculated from the first and last landmarks (or just first and second if there are only two), which is exactly what is needed.
cass_lat_all_ldk_bk<-coo_bookstein(cass_lat_all_ldk)
#Plot the aligned outlines to make sure the outlines are aligned along the same baseline
stack(cass_lat_all_ldk_bk) #(they are)

#The start points for the outlines are inconsistent. Here I create a single origin for each outline
#Write the coordinates out to a csv to create a vector for the coo_slide function, identifying which point is closest to the origin (0,0) of the aligned outlines using the csv in Microsoft Excel
write.csv(bookstein_test$coo,file='bookstein_coords_lat.csv')
#write a vector for sliding the coordinate start point to the same x-y coordinate for all outlines
#Create the list of points to shift all outline starting points to the same location
slide_vector<-c(200,76,259,277,289,135,147,178,242,159,223,202,339,258,275,286,281,303,384,249,389,181,167,114,273,141,282,68,65,156,134,127,5277,71,159,212,220,261,322,101,222,197,171,190,257,155,277,113,
  171,180,229,215,181,101,141,49,230,88,87,166,156,132,107,107,111,80,131,182,252,221,339,381,96,5247,163,159,267,73,246,91,234,3,29,133,141,257,219,160,90,192,121,145,31,208,92,200,103,89,71,238,199,170,75,138,170,192,237,137,174,269,209,134,224,197,233,144,1,178,289,239,328)
#Slide the start points of all outlines to ~ (0,0) using coo_slide
cass_lat_all_ldk_bk_slid<-coo_slide(cass_lat_all_ldk_bk,id=slide_vector)
#Verify the start points are relatively the same for all specimens
stack(cass_lat_all_ldk_bk_slid)
#Save the data object for redundancy
save(cass_lat_all_ldk_bk_slid,file='cass_lat_all_ldk_bk_slid.RData')

################################################################################
#Slicing the data to different analytical groups
#Left/Right analysis
cc_test_ldk_bk_slid<-slice(cass_lat_all_ldk_bk_slid,c(1,2,7,8,13,14,24,25,26,27,34,36,69,70,71,78,85,86,87,89,102:121))
#Inspect that the correct specimens were sliced into this group
cc_test_ldk_bk_slid
print(cc_test_ldk_bk_slid$fac[,3],n=40)
#generate graph to visually double-check the start points
stack(cc_test_ldk_bk_slid)
#save outline data object
save(cc_test_ldk_bk_slid,file='cc_test_ldk_bk_slid.RData')

#All lateral outlines (excluding the test data)
cc_lat_all_ldk_bk_slid<-slice(cass_lat_all_ldk_bk,c(1:101))
cc_lat_all_ldk_bk_slid_id_vector<-slide_vector[c(1:101)]
cc_lat_all_ldk_bk_slid<-coo_slide(cc_lat_all_ldk_bk_slid, id=cc_lat_all_ldk_bk_slid_id_vector)
#Inspect that the correct specimens were sliced into this group
cc_lat_all_ldk_bk_slid
print(cc_lat_all_ldk_bk_slid$fac[,1],n=100)
#generate graph to visually double-check the start points
stack(cc_lat_all_ldk_bk_slid)
#Save outline data object
save(cc_lat_all_ldk_bk_slid,file="cc_lat_all_ldk_bk_slid.RData")

#lateral sex analysis
cc_lat_sex_ldk_bk_slid<-slice(cass_lat_all_ldk_bk_slid, c(1:4,7,10:23,25,26,27,38,42,43,54:59,61,63,64,65,67,68,69,77,84,87:90,93:101))
#Inspect that the correct specimens were sliced into this group
cc_lat_sex_ldk_bk_slid
#Double check the sex data are accurate
print(cc_lat_sex_ldk_bk_slid$fac[,1],n=53)
#generate graph to visually double-check the start points
stack(cc_lat_sex_ldk_bk_slid)
#Save outline data object
save(cc_lat_sex_ldk_bk_slid,file='cc_lat_sex_ldk_bk_slid.RData')

#Lateral geography analysis
cc_lat_geo_ldk_bk_slid<-slice(cass_lat_all_ldk_bk_slid,c(1:4,6,7,10,11,25,30:33,35,36,38,40,41,43,45,49:52,54:60,62,63,69,76,77,81,84,87:90,95,97,98))
#Inspect that the correct specimens were sliced into this group
cc_lat_geo_ldk_bk_slid
#Double check the geography data are accurate
print(cc_lat_geo_ldk_bk_slid$fac[,2],n=45)
#generate graph to visually double-check the start points
stack(cc_lat_geo_ldk_bk_slid)
#Save outline data object
save(cc_lat_geo_ldk_bk_slid,file='cc_lat_geo_ldk_bk_slid_ed.RData')

#Lateral geography-sex analysis
cc_lat_geo_sex_ldk_bk_slid<-slice(cass_lat_all_ldk_bk_slid,c(1:4,10,11,25,43,54:59,77,84,87:90,95))
#Inspect that the correct specimens were sliced into this group
cc_lat_geo_sex_ldk_bk_slid
#Double check the sex and geography data are accurate
print(cc_lat_geo_sex_ldk_bk_slid$fac,n=21)
#generate graph to visually double-check the start points
stack(cc_lat_geo_sex_ldk_bk_slid)
#Save outline data object
save(cc_lat_geo_sex_ldk_bk_slid,file='cc_lat_geo_sex_ldk_bk_slid.RData')


################################################################################
###############################ANTERIOR LANDMARKING#############################
################################################################################

#######################Data read-in and outline creation########################
#Read in the factors list to categorize specimens
cass_info<-as_tibble(read.csv("FILE_PATH/all_cc_ant_fac_ldk.csv", header=T, row.names=1))
#Check the factors are correct
cass_info
#Create an object containing the filenames to reference them for import
lf_ant <- list.files("FOLDER_PATH/CC_ANT", full.names=TRUE)
#Double-check the filenames in the object
lf_ant
#Import the binary mask jpegs using the file list
cass_ant_all<-import_jpg(lf_ant)
#Convert them to outline, simultaneously adding factors to the outlines
cass_ant_all_out<-Out(cass_ant_all, fac=cass_info)
#Interpolate to the average number of points stated in the object information
cass_ant_all_out
cass_ant_all_out_int<-coo_interpolate(cass_ant_all_out,n=3626)
#Double check the point numbers are correct
cass_ant_all_out_int

###################Landmark Assignment and Outline Processing###################
##Make sure to save the object returned by the command
##There does not seem to be a way to go back through the specimens or edit points.
##Landmarks must be assigned in the same order for each specimen (I did anatomical right first, left second)
cass_ant_all_ldk<-def_ldk(cass_ant_all_out_int, 2)
#Inspect to make sure all specimens have 2 landmarks assigned
cass_ant_all_ldk
cass_ant_all_ldk$ldk
#Plot the outlines in their unaligned state, landmarks are red
stack(cass_ant_all_ldk)

#There are only two landmarks on the outlines at the ventral ends of the casques, so a Procrustes alignment cannot work. Bookstein coordinate alignment can, as it works on only two landmarks. These will be used to establish a new baseline calculated from the first and last landmarks (or just first and second if there are only two), which is exactly what is needed.
cass_ant_all_ldk_bk<-coo_bookstein(cass_ant_all_ldk)
#Plot the aligned outlines to make sure the outlines are aligned along the same baseline
stack(cass_ant_all_ldk_bk)#(they are)

#save this workspace
save.image("C:/Users/dik10/Documents/cass_ant_landmark_workspace.RData")

#The start points for the outlines are inconsistent. Here I create a single origin for each outline
#Write the coordinates out to a csv to create a vector for the coo_slide function, identifying which point is closest to the origin (0,0) of the aligned outlines using the csv in Microsoft Excel
write.csv(cass_ant_all_ldk_bk$coo,file='cc_ant_bookstein_coo.csv')
#write a vector for sliding the coordinate start point to the same x-y coordinate for all outlines
#Create the list of points to shift all outline starting points to the same location
slide_vector<-c(25,3446,3575,3492,3608,58,3616,3624,3608,27,3454,3567,3515,46,46,47,3589,83,51,31,3558,40,3584,3534,3530,46,3623,8,11,3619,3566,3606,3623,3550,3592,3485,3565,3532,3564,2,3552,3541,3612,3600,3567,3597,3617,3526,3623,19,31,3525,3539,3564,3617,3584,3510,3525,3499,3579,38,3532,3618,3585,3529,18,3611,3613,3614,54,3581,31,3573,3571,3610,3527,3603,3589,3586,3509,8,15,5,19,3610,32,3561)
#Slide the start points of all outlines to ~ (0,0) using coo_slide
cass_ant_all_ldk_bk_slid<-coo_slide(cass_ant_all_ldk_bk_slid,id=slide_vector)
#Verify the start points are relatively the same for all specimens
stack(cass_ant_all_ldk_bk_slid)
#Save the data object for redundancy
save(cass_ant_all_ldk_bk_slid,file='cass_ant_all_ldk_bk_slid.RData')

################################################################################
#slicing the data to different analytical populations
#All anterior outlines
cc_ant_all_ldk_bk_slid<-cass_ant_all_ldk_bk_slid
#Inspect that the correct specimens were sliced into this group
cc_ant_all_ldk_bk_slid
#generate graph to visually double-check the start points
stack(cc_ant_all_ldk_bk_slid)
#Save outline data object
save(cc_ant_all_ldk_bk_slid,file="cc_ant_all_ldk_bk_slid.RData")

#anterior sex analysis
cc_ant_sex_ldk_bk_slid<-slice(cass_ant_all_ldk_bk_slid,c(1:4,7,11:26,29,41,53:58,60,61,63:67,69,70,71,76,78,85))
#Inspect that the correct specimens were sliced into this group
cc_ant_sex_ldk_bk_slid
#Double check the sex data are accurate
print(cc_ant_sex_ldk_bk_slid$fac[,1],n=42)
#generate graph to visually double-check the start points
stack(cc_ant_sex_ldk_bk_slid)
#Save outline data object
save(cc_ant_sex_ldk_bk_slid,file='cc_ant_sex_ldk_bk_slid.RData')

#anterior geography analysis
cc_ant_geo_ldk_bk_slid<-slice(cass_ant_all_ldk_bk_slid,c(1:4,6,7,11,12,29:33,35,37,39,40,41,43,47,48,50,51,53:59,61,62,65,71,77,78,85))
#Inspect that the correct specimens were sliced into this group
cc_ant_geo_ldk_bk_slid
#Double check the geography data are accurate
print(cc_ant_geo_ldk_bk_slid$fac[,2],n=37)
#generate graph to visually double-check the start points
stack(cc_ant_geo_ldk_bk_slid)
#Save outline data object
save(cc_ant_geo_ldk_bk_slid,file='cc_ant_geo_ldk_bk_slid.RData')

#anterior geography-sex analysis
cc_ant_geo_sex_ldk_bk_slid<-slice(cass_ant_all_ldk_bk_slid, c(1:4,11,12,29,41,53:58,78,85))
#Inspect that the correct specimens were sliced into this group
cc_ant_geo_sex_ldk_bk_slid
#Double check the sex and geography data are accurate
print(cc_ant_geo_sex_ldk_bk_slid$fac, n=16)
#generate graph to visually double-check the start points
stack(cc_ant_geo_sex_ldk_bk_slid)
#Save outline data object
save(cc_ant_geo_sex_ldk_bk_slid,file='cc_ant_geo_sex_ldk_bk_slid.RData')
