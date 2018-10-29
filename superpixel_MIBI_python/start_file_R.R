## 21st Dec 2016
## BISCUIT R implementation
## Start_file with user inputs
## 
## Code author SP
##
## 6th Feb 2018
## Added code to iterate clustering of pixet1 and Bayesian classification of pixet2.
##
###
###

############## packages required ##############

library(MCMCpack)
library(mvtnorm)
library(ellipse)
library(coda)
library(Matrix)
library(Rtsne)
library(gtools)
library(foreach)
library(doParallel)
library(doSNOW)
library(snow)
library(lattice)
library(MASS)
library(bayesm)
library(robustbase)
library(chron)
library(mnormt)
library(schoolmath)
library(RColorBrewer)
library(R.matlab)
library(parallel)
library(stats)
library(dplyr)
library(tiff)
#install.packages("rapport")
library(rapport)
library(gplots)



#path_t <- "/Users/prabhaks/Documents/Multiplet_Detection/Python_code/MIBI_model"
path_t <- getwd(); #use this in the python call code
working_path <- paste0(path_t,"/superpixel_MIBI_python")

setwd(working_path)

print(working_path)

input_data_tab_delimited <- TRUE; #set to TRUE if the input data is tab-delimited

is_format_genes_cells <-  FALSE; #set to TRUE if input data has rows as genes and columns as cells

num_iter <- 15; #number of iterations, choose based on data size.

num_cores <- detectCores() - 4; #number of cores for parallel processing. Ensure that detectCores() > 1 for parallel processing to work, else set num_cores to 1.

z_true_labels_avl <- FALSE; #set this to TRUE if the true labels of cells are available, else set it to FALSE. If TRUE, ensure to populate 'z_true' with the true labels in 'BISCUIT_process_data.R'

num_cells_batch <- 100; #set this to 1000 if input number of cells is in the 1000s, else set it to 100.

alpha <- 0.1; #DPMM dispersion parameter. A higher value spins more clusters whereas a lower value spins lesser clusters.



#full.data.0 <- readMat(input_file_name);
#print(dim(full .data.0))
#print("hello world")


#build the image tensor
channels_used <- read.csv(paste0(working_path,"/count_data/clusterChannels.csv"),header=FALSE)
lc <- length(channels_used[[1]])
choose_genes <- lc; #comment if you want all the genes to be considered

gene_names <- as.character(unlist(channels_used))

gene_batch <- lc; #number of genes per batch, therefore num_batches = choose_genes (or numgenes)/gene_batch. Max value is 150

iter_pixets <- 1
within_iter_model <- 1;
input_file_name <- paste0(working_path,"/count_data/count_matrix_LM.csv");
output_folder_name <- "/output_6_markers_Leeat_original_py"; #give a name for your output folder.

final_clustering <- 0
#do.L_CM <- TRUE
source("BISCUIT_main.R")

meta_data <- read.table(paste0(working_path,"/count_data/meta_data.txt"))
path_channels <- meta_data[[1]][1]
patient_num <- meta_data[[1]][2]
image_dim <- strtoi(meta_data[[1]][3])
max_MCMC <- strtoi(meta_data[[1]][4])

path_marker <- paste0(path_channels[[1]],"/",channels_used[[1]][1])

image_tensor <- array(rep(0,image_dim*image_dim*lc), dim=c(image_dim,image_dim,lc)) #"countsReshape.mat"

for (i in 1:lc){
    path_marker <- paste0(path_channels[[1]],"/Point",patient_num,"/",channels_used[[1]][i],".tif")
    print(path_marker)
    t=readTIFF(path_marker,as.is=TRUE)
    image_tensor[,,i] = t[1:image_dim,1:image_dim]
}

#read in datasets needed to start clustering and reassigning between superpixet1 and superpixet2
	
patient_file_name <- paste0(working_path,"/count_data/");
input_file_name <- paste0(patient_file_name,"count_matrix_6_updated.csv"); # this is the count matrix created from pixet1 (i.e. nuclear pixels) and for just once through Matlab

valid_WSS_file_name <- paste0(patient_file_name,"valid_WSS.csv");
valid_WSS <- as.matrix(read.csv(valid_WSS_file_name, header=FALSE))
valid_WSS <- valid_WSS+1
#choose_cells <- length(valid_cells) #4876 #4988; #comment if you want all the cells to be considered

valid_cells_file_name <- paste0(patient_file_name,"valid_cells.csv");
valid_cells <- as.matrix(read.csv(valid_cells_file_name, header=FALSE))
valid_cells <- valid_cells + 1

super_pixel_file_name <- paste0(patient_file_name, "count_matrix_super_pixel.csv");
cm_superpixel <- as.matrix(read.csv(super_pixel_file_name,header=FALSE));

choice <- dim(cm_superpixel)[1]

WSS_celllabels_file_name <- paste0(patient_file_name,"WSS_celllab_assgn.csv");
WSS_cell_labels <- as.matrix(read.csv(WSS_celllabels_file_name, header=FALSE))
WSS_cell_labels <- WSS_cell_labels + 1;

neighbours_super_pixel_file_name <- paste0(patient_file_name, "neighbours_super_pixel.csv");
neigh_superpixel <- as.matrix(read.csv(neighbours_super_pixel_file_name,header=FALSE));
neigh_superpixel <- gsub("\\[|\\]","",neigh_superpixel)

#sd <- as.numeric(unlist(strsplit(ns[1]," ")))
#sd <- sd[!is.na(sd)]
# add + 1

pixet_1_struct_file_name <- paste0(patient_file_name, "pixet_1_struct.csv");
pixet_1_struct_py <- as.matrix(read.csv(pixet_1_struct_file_name, header=FALSE));
pixet_1_struct_py <- gsub("\\[|\\]","",pixet_1_struct_py)

rows_pixet_1_str <- nrow(pixet_1_struct_py)
pixet_1_struct <- list()

for (i in 1:rows_pixet_1_str){
    ps <- as.numeric(unlist(strsplit(pixet_1_struct_py[i]," ")))
    ps <- ps[!is.na(ps)]
    ps <- t(matrix(ps,2,length(ps)/2))
    pixet_1_struct[[i]] <- ps + 1
}

pixet_1_struct_cp <- pixet_1_struct;

#read in the stats_revised and stats_WSS_revised
#stats_revised_area
rra_file_name <- paste0(patient_file_name,"regions_revised_area.csv");
rr_area <- as.matrix(read.csv(rra_file_name, header=FALSE))
#stats_revised_centroid
rrcent_file_name <- paste0(patient_file_name, "regions_revised_centroid.csv");
rr_centroid <- as.matrix(read.csv(rrcent_file_name, header=FALSE));
rr_centroid <- rr_centroid + 1
#stats_revised_pixelcoords
rrpcoord_file_name <- paste0(patient_file_name, "regions_revised_pixelcoords.csv");
rr_pcoord_py <- as.matrix(read.csv(rrpcoord_file_name, header=FALSE));
rr_pcoord_py <- gsub("\\[|\\]","",rr_pcoord_py)
#repeat pixet_1_struct circus below

rows_rrp <- nrow(rr_pcoord_py)
rr_pcoord <- list()

for (i in 1:rows_rrp){
    rrp <- rr_pcoord_py[i]
    rrp <- as.numeric(unlist(strsplit(rrp," ")))
    rrp <- rrp[!is.na(rrp)]
    rr_pcoord[[i]] <- t(matrix(rrp,2,length(rrp)/2))
}

rr_pcoord_cp <- rr_pcoord

#cleaned and corrected WSS_revised_area
rWSSra_file_name <- paste0(patient_file_name,"WSS_revised_area.csv");
WSSr_area <- as.matrix(read.csv(rWSSra_file_name, header=FALSE))
#stats_WSS_revised_centroid
rWSSrcent_file_name <- paste0(patient_file_name, "WSS_revised_centroid.csv");
WSSr_centroid <- as.matrix(read.csv(rWSSrcent_file_name, header=FALSE));
WSSr_centroid <- WSSr_centroid + 1
#stats_WSS_revised_pixelcoords
rWSSrpcoord_file_name <- paste0(patient_file_name, "WSS_revised_pixelcoords.csv");
WSSr_pcoord_py <- as.matrix(read.csv(rWSSrpcoord_file_name, header=FALSE));
WSSr_pcoord_py <- gsub("\\[|\\]","",WSSr_pcoord_py)
#repeat pixet_1_struct circus below
rows_WSSrp <- nrow(WSSr_pcoord_py)
WSSr_pcoord <- list()

for (i in 1:rows_WSSrp){
    rrp <- WSSr_pcoord_py[i]
    rrp <- as.numeric(unlist(strsplit(rrp," ")))
    rrp <- rrp[!is.na(rrp)]
    WSSr_pcoord[[i]] <- t(matrix(rrp,2,length(rrp)/2))
}


#do.L_CM <- FALSE; # now we start clustering superpixet1 and iter_pixets is still 1
c_zero <- list(0)
d_WSS_area <- dim(WSSr_area)[1]
reassin_frac <- list(0)

#while(iter_pixets <=max_MCMC){ #10,max_MCMC
while(iter_pixets <=14 && (length(which(WSSr_area==0))/d_WSS_area)*100 <= 95){ #breaks when 95% of superpixet2 are assigned to superpixet1
    print(paste("iteration # ",iter_pixets));
    mcmc_counter <- iter_pixets
    
    reass_frac_val <- (length(which(WSSr_area==0))/d_WSS_area)*100
    print(paste(reass_frac_val,"% of superpixet2 reassigned to superpixet1 "));
    reassin_frac[[iter_pixets]] <- reass_frac_val
    
    alpha <- 0.1;
    output_folder_name <- paste0("output_6_markers_updated_1_",iter_pixets); #give a name for your output folder.
    within_iter_model <- within_iter_model + 1;
    source("BISCUIT_main.R"); #for clustering superpixet1

    if(iter_pixets >1){ #for spatially-aware clustering of superpixet2
        # rm(uq_z);
        # rm(z_inferred_final);
        
        print("Entering MIBI model");
        #if(iter_pixets ==1){
        #   source("MIBI_superpixel_model_2_0.R");
        #}
        source("MIBI_model_main.R"); # add code if you have helper functions
    }
    iter_pixets <- iter_pixets + 1;
}

write.table((iter_pixets-1),paste0(working_path,"/count_data/meta_data.txt"),append=TRUE,row.names=FALSE,col.names=FALSE,sep="")

write.csv(unlist(reassin_frac)[-1],paste0(working_path,"/count_data/reas_frac.txt"),row.names=FALSE)
write.csv(unlist(c_zero)[-1],paste0(working_path,"/count_data/sp2tosp1.txt"),row.names=FALSE)


#final clustering
mcmc_counter <- mcmc_counter + 1;
final_clustering <- 1
output_folder_name <- paste0("output_6_markers_updated_1_",iter_pixets); #give a name for your output folder.
within_iter_model <- within_iter_model + 1;
source("BISCUIT_main.R"); #for clustering superpixet


#plot(1:length(unlist(c_zero)),unlist(c_zero),typ='l')
#pie(unlist(reassin_frac)[-1])
#plot(1:length(unlist(reassin_frac)),unlist(reassin_frac),typ='l')


source("code_snippets_forMIBI_1.R")
