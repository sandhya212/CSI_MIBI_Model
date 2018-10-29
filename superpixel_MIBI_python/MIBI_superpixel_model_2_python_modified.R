
## 31st Jan 2018
## MIBI_model_Bayesian classifier
##
## Code author SP
##
##
## pixet_1 are the set of certain pixels - they are nuclear pixels
## pixet_2 is every other pixel
## MAP per pixel per cluster
##
#### 1st Feb 2018 #####
## MLE per pixel per cluster
##
## 5th Feb 2018
## Euclidean dist based features
##
## 6th Feb 2018
## Updating pixet1 and pixet2 based on NBayes classifier MAP estimates
##
## 7th Feb 2018
## Updated stats_revised with pixet1 and pixet2 movements
## Corrected p logic = pixet_2[p,3].
##
## 15th Feb 2018
## Revised the stats_revised structure for Area and Centroids post pixel assignment
##
## 16th Feb 2018
## accounting for the case when a cluster has just one pixel
##
## 21st Feb 2018
## Added code for rebuilding stats_revised for R Matlab swings
##
## 28th Feb 2018
## Upgrading the prior for clusters to be cell dependent for feature similarity based on neighbours and also capturing inverse dependency of neighbour distances via exponential function
## 16-19th June 2018


# to pull out the correct columns from mu and Sigma per cluster
# dimnames gives us the clusters in order (1....K);

#col_names <- dimnames(mu_final)[[2]]
#col_names <- col_names[!col_names=='dummy']
#col_names

#tabulate gives us the clusters in order (1....K);
#tab_z <- tabulate(z_inferred_final);
#k_prior <- log(tab_z/dim(full.data.3)[1]);
        #k_prior <- tab_z/dim(full.data.3)[1];



# use a list ?
#choice <- dim(cm_superpixel)[1]

#if (iter_pixets!=1){
    col_names <- dimnames(mu_final)[[2]]
    col_names <- col_names[!col_names=='dummy']
    col_names
    
    #tabulate gives us the clusters in order (1....K);
    tab_z <- tabulate(z_inferred_final);
    k_prior <- log(tab_z/dim(full.data.3)[1]);
    #k_prior <- tab_z/dim(full.data.3)[1];
    
    posterior_superpixet_2 <- matrix(0, choice,final_num_K);
    colnames(posterior_superpixet_2) = col_names;
#}
#this module is called only when iter_pixets >=2

llhood_superpixet_2 <- matrix(0, choice,final_num_K);
colnames(llhood_superpixet_2)= col_names;

MAP_superpixet_2 <- matrix(0,choice,1);
MLE_superpixet_2 <- matrix(0,choice,1);
######check if 1 is needed, interesting indeed ###

##Bayesian MAP prediction per pixel per cluster
## TO DO
## https://www.cs.waikato.ac.nz/~eibe/pubs/UAI_200.pdf for weighted NB
## GP classifier

w <- 1e-3 #1e-6 # or less or more?

counter <- 1;
euc_dist_p <- 0;
euc_dist_p1 <- 0;
count_zeroing <- 0;
dist_method <- 'euclidean'; #'minkowski'

# this is done for each 'uncertain' superpixel
for (sp in 1:choice){ ##sp is a counter for superpixels. WWWWWWWto access the actual superpixel number, do pixet_2[p,3]
    
     print(paste("sp is ", sp));
     
     if(WSSr_area[sp]!=0){
        
        rm(euc_dist_p);
        rm(euc_dist_p1);
        
        #neigh_sp = unlist(neigh_superpixel[2,sp,]); #who are sp's neighbours
        neigh_sp <- as.numeric(unlist(strsplit(neigh_superpixel[sp]," ")))
        neigh_sp <- neigh_sp[!is.na(neigh_sp)]
        neigh_sp <- neigh_sp + 1
        print(paste('neighbours are ',neigh_sp));
        
        #WSS_rows <- which(unlist(stats_WS_superpixel[5,,1]) %in% neigh_sp); #which of the neighbours are in superpixet2, find the rows
        WSS_rows <- which(valid_WSS %in% neigh_sp);
        print(paste('rows where sp neighbours are superpixet2 ', WSS_rows));
        
        
        #pixet1_labs <- unlist(stats_WS_superpixel[6,WSS_rows,]); #find those superpixet1 rows (rows here = those in stats_revised) as that is how it is created in Matlab
        pixet1_labs <- WSS_cell_labels[WSS_rows]
        print(paste('true pixet1_labs  ', pixet1_labs));
        
        ## find the Euclidean distance between all of sp's neighbours
        ## irrespective of them being in superpixet1 or not
        
        #### >>>> # compute Euclidean distance of sp to neighbouring superpixels' centroids
        #x <- as.numeric(unlist(stats_WS_superpixel[2,sp,1]));
        x <- as.numeric(WSSr_centroid[sp,])
        euc_dist_p1 <- 0;
        #neigh_rev_ind <- unlist(stats_WS_superpixel[5,WSS_rows,]);#as.numeric(neigh_sp) #matrix(0,length(neigh_sp),1);
        neigh_rev_ind <- neigh_sp;
        
        for (s in 1:length(neigh_sp)){
            #compute Euc distance of sp to cell centroids of neighbourhoods
            #cell_label_temp <- unlist(neigh_sp[s]);
            
            ##
            #neigh_rev_ind[s] <- unlist(stats_WS_superpixel[5,WSS_rows[s],]); #which(stats_revised[5,,]==cell_label_temp)
            ##
            #y <- as.numeric(c(unlist(stats_WS_superpixel[2,WSS_rows[s],1])[1],unlist(stats_WS_superpixel[2,WSS_rows[s],1])[2]));
            #y <- as.numeric(unlist(stats_revised[2,cell_row[j],1]));
            y <- as.numeric(WSSr_centroid[WSS_rows[s],])

            mat <- rbind(x,y)
            euc_temp <- as.numeric(dist(mat,dist_method));
            euc_dist_p1 <- rbind(euc_dist_p1,euc_temp);
        }
        #sp is also in the neighbour set of p
        #y <- as.numeric(c(unlist(stats_WS_superpixel[2,sp,1])[1],unlist(stats_WS_superpixel[2,sp,1])[2]));
        #mat <- rbind(x,y);
        #euc_temp = as.numeric(dist(mat,'euclidean'));
        #euc_dist_p1 <- rbind(euc_dist_p1,euc_temp);
        euc_dist_p1 <- euc_dist_p1[2:length(euc_dist_p1)];
        
        print(euc_dist_p1);
        
        
        
        ######
        sum_pixet1_labs <- sum(pixet1_labs);
        if (sum_pixet1_labs>0){ # get into the processing for current superpixel only if has some pixet1 neighbours
        
            ##find classes of all pixet1 that are sp's neighbours
            px_len <- length(pixet1_labs); #if this is 0, then this is because non of the sp neighbours belong to pixet1
            class_px <- matrix(0,1,px_len);
        
            for (i in 1:px_len){
                # x <- vector(mode="numeric", length=0)
                # length(x) = 0
                print(paste(i,"th neighbour"))
                if (pixet1_labs[i]>0){
                    #cell_label_row <- which(stats_revised[5,,]==pixet1_labs[i])
                    #print(paste("cell_label_row is ", cell_label_row));
                    #class_px[i] <- cnames[z_inferred_final[cell_label_row]];
                    #print(paste("z is", z_inferred_final[cell_label_row]));
                    ###CHECK
                    ####
                    #superpixel1 <- as.numeric(stats_revised[5,pixet1_labs[i],]);
                    
                    #superpixel1 <- as.numeric(which(stats_revised[5,,]==pixet1_labs[i]));
                    superpixel1 <- as.numeric(which(valid_cells==pixet1_labs[i])); #this will be the row of the valid cell id
                    if (length(superpixel1!=0)){
                        class_px[i] <- cnames[z_inferred_final[superpixel1]];
                        print(paste("z is", z_inferred_final[superpixel1]));
                        print(paste("class is" , class_px[i]));
                        
                    }
                }
            
            }
            
            
            #find sp likelihood
            for (k in 1:final_num_K){
                
                ####Likelihood
                superpixet_2_loglkhd <-  dmnorm(cm_superpixel[sp,],mu_final[,col_names[k]],Sigma_final[,,col_names[k]][1]);
                sum_superpixet_2_ll <- sum(superpixet_2_loglkhd);
                llhood_superpixet_2[sp,col_names[k]] <- sum_superpixet_2_ll;
                
                
                
                ####Spatial prior
                ####Cluster k prior p(k) is replaced with a prior taking into account the pixel's neighbours. This spatially-sensitive prior is therefore 'dependent' upon each pixel (or cell). (Earlier p(k) was agnostic to its neighbours, therefore independent). Every pixel will now have a linearly-smoothed contribution from its neighbours and itself. This will also be balanced using distances as in inverse-square law.
                spat_prior <- 0.00001
                exp_spat_prior <- exp(spat_prior); #this alone would be like adding 1 to the variable
                
                ##
                ##
                ##
                ######
                if(length(neigh_sp)!=0){
                   for (nei in 1:length(neigh_sp)){
                       #spat_prior <- spat_prior + (1/euc_dist_p1[nei]) * posterior_superpixet_2[WSS_rows[nei],col_names[k]];
                       #row_WSS <- which(stats_WS_superpixel[5,,]==neigh_sp[nei]);
                       #spat_prior <- spat_prior + (1/euc_dist_p1[nei]) * posterior_superpixet_2[row_WSS,col_names[k]];#31st Aug
                       spat_prior <- spat_prior + (1/euc_dist_p1[nei]);

                    }
                   #spat_prior <- spat_prior + posterior_superpixet_2[sp,col_names[k]];
                    exp_spat_prior <- exp(w/(length(neigh_sp)) * spat_prior);
                }
               
                posterior_superpixet_2[sp,col_names[k]] = sum_superpixet_2_ll + exp_spat_prior;
            } #end of for loop with k
            
            MAP_superpixet_2[sp] <- cnames[which.max(posterior_superpixet_2[sp,])];
            MLE_superpixet_2[sp] <- cnames[which.max(llhood_superpixet_2[sp,])];
            
            #check if sp's current class is present in any of its pixet1 neighbours
            if (length(intersect(class_px,MAP_superpixet_2[sp]))>0){
                
                # compute Euclidean distance of sp to neighbouring superpixels' centroids
                centroid_sp <- x; #as.numeric(unlist(stats_WS_superpixel[2,sp,1]));
                euc_dist_p <- 0;
                
                #find which neighbour it is closest to
                nei_sim_z <- which(class_px %in% MAP_superpixet_2[sp])
                print(paste('neighbours with same class ',neigh_sp[nei_sim_z]));
                
                ln_nei <- length(nei_sim_z);
                cell_row <- matrix(0,1,ln_nei);
                for (j in 1:ln_nei){
                    #cell_row[j] <- which(stats_revised[5,,]==pixet1_labs[nei_sim_z[j]])
                    #cell_row[j] <- which(stats_revised[5,,]==pixet1_labs[nei_sim_z[j]]);#pixet1_labs[nei_sim_z[j]]
                    cell_row[j] <- which(valid_cells==pixet1_labs[nei_sim_z[j]]);#pixet1_labs[nei_sim_z[j]]
                    if(length(cell_row[j])!=0){
                        #y <- as.numeric(unlist(stats_revised[2,cell_row[j],1]));
                        y <- as.numeric(unlist(rr_centroid[cell_row[j],]));
                        mat <- rbind(centroid_sp,y)
                        euc_temp <- as.numeric(dist(mat,dist_method));
                        euc_dist_p <- rbind(euc_dist_p,euc_temp);
                    }
                    
                }
                euc_dist_p <- euc_dist_p[2:length(euc_dist_p)];
                print(paste("Euclidean distance to closest cell is ", euc_dist_p));

                min_row <- which.min(euc_dist_p)
                print(paste("min row is: ",min_row));
                print(paste("closest neighbour : ",neigh_sp[nei_sim_z[min_row]]));
                
                #min_cell_row <- which(stats_revised[5,,]==pixet1_labs[nei_sim_z[min_row]]);
                min_cell_row <- which(valid_cells==pixet1_labs[nei_sim_z[min_row]]);
                
                print(paste("min cell row is: ",min_cell_row));
                
                
                #remove sp from stats_WS_superpixel
                # update sp to stats_revised cell_label to which it is closest to
                if(length(min_cell_row)!=0){ #which cant be since this is the cell label aka the pixet1 label
                    count_zeroing <- count_zeroing + 1;
                    source("Rearrange_MIBI_datastructs.R")
                }

            }
                
        }# end 'if' when you do processing for sp with pixet1 neighbours
     }#end 'if' if area is zero
}# end 'for'




#plots
apply(posterior_superpixet_2[1:choice,],1, order)
apply(llhood_superpixet_2[1:choice,],1, order)

#
table(MAP_superpixet_2)
table(MLE_superpixet_2)


#####


##write the imputed matrix, still in logspace
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/pixet_1.txt");
#write.matrix(pixet_1,file=f,sep="\t")
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/pixet_1.mat");
#writeMat(f,pixet_1 = pixet_1);
##
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/pixet_2.txt");
#write.matrix(pixet_2,file=f,sep="\t")
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/pixet_2.mat");
#writeMat(f,pixet_2 = pixet_2);
##
###################################################
######### Below IS WHAT WE HAD BEFORE ##############
########Commented out on 9th Sep 2018 to write Python-specific csvs
###################################################
#?f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/stats_revised.txt");
#?write.matrix(stats_revised,file=f,sep="\t")
#better
#s_area <- unlist(stats_revised[1,,1]);
#l_s <- length(unlist(stats_revised[2,,1]))/2;
#s_centroid <- t(matrix(unlist(stats_revised[2,,1]),2,l_s));
#s_cl <- unlist(stats_revised[5,,1]);
#s_pixelind <- unlist(stats_revised[3,,1]);
#s_pixelxy <- matrix(0,1,2);
#for (d in 1 :choose_cells){
#    s_pixelxy_temp <- matrix(unlist(stats_revised[4,d,1]),s_area[d],2);
#    s_pixelxy <- rbind(s_pixelxy,s_pixelxy_temp);
#}
#s_pixelxy <- s_pixelxy[-1,]

#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/s_area.mat");
#writeMat(f,Area = s_area);
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/s_centroid.mat");
#writeMat(f,Centroid = s_centroid);
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/s_cl.mat");
#writeMat(f,cell_label = s_cl);
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/s_pixelxy.mat");
#writeMat(f,PixelList = s_pixelxy);
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/s_pixelind.mat");
#writeMat(f,PixelIdxList = s_pixelind);
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/full_data_3.mat");
#writeMat(f,fulldata = full.data.3);
###################################################
###################################################

f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/full_data_3.csv");
#write.table(as.matrix(full.data.3),f,row.names=FALSE,col.names=FALSE);
write.csv(as.matrix(full.data.3),f,row.names=FALSE);

f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/regions_revised_area_updated.csv");
write.table(rr_area,f,row.names=FALSE,col.names=FALSE);

f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/regions_revised_centroid_updated.csv");
#write.table(rr_centroid,f,row.names=FALSE,col.names=FALSE);
write.csv(as.matrix(rr_centroid),f,row.names=FALSE);

s_pixelxy <- matrix(0,1,2);
for (d in 1 :choose_cells){
    s_pixelxy_temp <- matrix(rr_pcoord[[d]],rr_area[d],2);
    s_pixelxy <- rbind(s_pixelxy,s_pixelxy_temp);
}
s_pixelxy <- s_pixelxy[-1,]

f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/regions_revised_pixelcoords_updated.csv");
#write.table(s_pixelxy,f,row.names=FALSE,col.names=FALSE);
write.csv(as.matrix(s_pixelxy),f,row.names=FALSE);

print(paste("Number of superpixet2 moved to superpixet1: ",count_zeroing));
c_zero[[iter_pixets]] <- count_zeroing;
