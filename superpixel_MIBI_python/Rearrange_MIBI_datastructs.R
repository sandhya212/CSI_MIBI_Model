## 28th Feb 2018
## Rearranging the MIBI data structures based on closest neighbour to pixel under consideration.
##
## Code author SP
##

###
###

#pixet_2[p,4] <- closest_neigh;

### update the stats_revised for p
### i.e. remove the pixedindlist and pixellist (px,py) for pth pixel. Area, centroid and cell_label stay the same.

#print(paste('Zero out area of superpixel ',stats_WS_superpixel[5,sp,],' stats_WS_superpixel '));
print(paste('Zero out area of superpixel ',valid_WSS[sp],'from stats of WS_superpixel '));


#find which entry is p in the previous cell
#d1 <- which(as.numeric(unlist(stats_revised[3,stats_rev_ind,1]))==pixet_2[p,3]);

#if(is.matrix(stats_revised[3,stats_rev_ind,1][[1]])){
#   indi <- logical(length=dim(stats_revised[3,stats_rev_ind,1][[1]])[1]);
#}else{
#    indi <- logical(length=length(stats_revised[3,stats_rev_ind,1][[1]]));
#}
#indi = !indi; #so now it is all TRUE
#indi[d1] = FALSE;

#####stats_revised[3,stats_rev_ind,1][[1]][d1] = 0; #set this to 0 or erase it. See which is better



#consider the case when that stats entry is all gone....


#find the corresponding (px,py) for pth pixel and set them to 0 as well.
#if(length(stats_revised[4,stats_rev_ind,1][[1]])==2){ #this can also be checked with stats_revised[3...] meaning this entry will be of length 1, just one pixel
#   d1_px <- stats_revised[4,stats_rev_ind,1][[1]][1]; #px
#   d1_py <- stats_revised[4,stats_rev_ind,1][[1]][2]; #py
    #stats_revised[4,stats_rev_ind,1][[1]][1] <- 0;
    #stats_revised[4,stats_rev_ind,1][[1]][2] <- 0;
    
    #27th Feb 2018 commented !!!! grrr
    
    #reversed the order with first 2 and then 4, so that you at least have the last pixel info as centroid, else it is 0 and that is not right to make the euc distane calculation.
#   stats_revised[2,stats_rev_ind,1][[1]] <- stats_revised[4,stats_rev_ind,1][[1]];
#   stats_revised[4,stats_rev_ind,1][[1]] <- stats_revised[4,stats_rev_ind,1][[1]][indi];
    
    #stats_revised[2,stats_rev_ind,1][[1]][indi]; #stats_revised[4,stats_rev_ind,1][[1]]; # do we keep the pixel here? no we remove it
#   stats_revised[1,stats_rev_ind,1][[1]] <- 0
    
#}else{
#   d1_px <- stats_revised[4,stats_rev_ind,1][[1]][d1,1]; #px
#   d1_py <- stats_revised[4,stats_rev_ind,1][[1]][d1,2]; #py
    #stats_revised[4,stats_rev_ind,1][[1]][d1,1] <- 0;
    #stats_revised[4,stats_rev_ind,1][[1]][d1,2] <- 0;
#   stats_revised[1,stats_rev_ind,1][[1]] <- stats_revised[1,stats_rev_ind,1][[1]] - 1;
#   stats_revised[4,stats_rev_ind,1][[1]] <- stats_revised[4,stats_rev_ind,1][[1]] [indi, ]; #remove that pixel's coordinates
    
    # this case takes care when after the previous deletion we hit one entry alone
#   if(length(stats_revised[4,stats_rev_ind,1][[1]])==2){
#       stats_revised[2,stats_rev_ind,1][[1]] <- stats_revised[4,stats_rev_ind,1][[1]];
#   }else{
#       new_centroids <- colSums(stats_revised[4,stats_rev_ind,1][[1]])/stats_revised[1,stats_rev_ind,1][[1]];
#       stats_revised[2,stats_rev_ind,1][[1]] <- new_centroids;
#   }
    
#
#}

####Update stats_revised#####
#############################
#### Decrease the area by 1, recompute centroid of stats_rev_ind

#recompute centroid

#if(length(stats_revised[4,stats_rev_ind,1][[1]])==2){
#stats_revised[2,stats_rev_ind,1][[1]] <- stats_revised[4,stats_rev_ind,1][[1]]
#stats_revised[1,stats_rev_ind,1][[1]] <- 1

#}else{

#ind <- rowSums(stats_revised[4,stats_rev_ind,1][[1]]  == 0) != ncol(stats_revised[4,stats_rev_ind,1][[1]] )
#stats_revised[4,stats_rev_ind,1][[1]] <- stats_revised[4,stats_rev_ind,1][[1]] [ind, ]; #remove that pixel's coordinates
#stats_revised[1,stats_rev_ind,1][[1]] <- stats_revised[1,stats_rev_ind,1][[1]] - 1
#new_centroids <- colSums(stats_revised[4,stats_rev_ind,1][[1]])/stats_revised[1,stats_rev_ind,1][[1]];
#stats_revised[2,stats_rev_ind,1][[1]] <- new_centroids;
#}

#remove that pixel entry
#stats_revised[3,stats_rev_ind,1][[1]] <- stats_revised[3,stats_rev_ind,1][[1]] [ind]; #remove that pixel id


# update sp to stats_revised cell_label to which it is closest to
#######add this to the new cell assignment, this is not a new cell...it is there already, we are only assigning p to this cell
#print(paste('add superpixel ',stats_WS_superpixel[5,sp,],' to cell ',pixet1_labs[nei_sim_z[min_row]], ' in stats_revised (pixelindlist and pixelidlist)'));
print(paste('add superpixel ',valid_WSS[sp],' to cell ',pixet1_labs[nei_sim_z[min_row]], ' in stats_revised (pixelindlist and pixelidlist)'));



#stats_revised[2,min_cell_row,][[1]] <- colSums(rbind(stats_revised[2,min_cell_row,][[1]],stats_WS_superpixel[2,sp,][[1]]))/2;#pixel cent coord
rr_centroid[min_cell_row,] <- colSums(rbind(as.numeric(rr_centroid[min_cell_row,]),as.numeric(WSSr_centroid[sp,])))/2;#pixel cent coord


#stats_revised[3,min_cell_row,][[1]] <- unique(append(stats_revised[3,min_cell_row,][[1]],stats_WS_superpixel[3,sp,][[1]]));#pixelindlist
#dont have this in Python code anymore

#stats_revised[4,min_cell_row,][[1]] <- unique(rbind(stats_revised[4,min_cell_row,][[1]],stats_WS_superpixel[4,sp,][[1]]));#pixelcoord

#rrp <- rr_pcoord[min_cell_row]
#rrp <- as.numeric(unlist(strsplit(rrp," ")))
#rrp <- rrp[!is.na(rrp)]
#rrp <- t(matrix(rrp,2,length(rrp)/2))
#rrp <- rrp + 1

#WSSp <- WSSr_pcoord[sp]
#WSSp <- as.numeric(unlist(strsplit(WSSp," ")))
#WSSp <- WSSp[!is.na(WSSp)]
#WSSp <- t(matrix(WSSp,2,length(WSSp)/2))
#WSSp <- WSSp + 1

rr_pcoord[[min_cell_row]] <- unique(rbind(rr_pcoord[[min_cell_row]],WSSr_pcoord[[sp]],pixet_1_struct_cp[[min_cell_row]]-1));#pixelcoord

#stats_revised[1,min_cell_row,][[1]] <- stats_revised[1,min_cell_row,][[1]] + stats_WS_superpixel[1,sp,][[1]]; #area
#stats_revised[1,min_cell_row,][[1]] <- length(stats_revised[3,min_cell_row,][[1]]);  #area
rr_area[min_cell_row] <- nrow(rr_pcoord[[min_cell_row]]);  #area

#I AM NOT SURE IF WE NEED THIS. COMMENTING THIS OUT. On 9th Sep 2018
# Since we have a different data structure to deal with Python-based readins
#####
####update pixet_1_struct for the min_cell
####mincell_row_pixet1strct <- which(as.numeric(pixet_1_struct[1,,])==pixet1_labs[nei_sim_z[min_row]]);
#mincell_row_pixet1strct <- pixet1_labs[nei_sim_z[min_row]]; (not the right one)
####if(length(mincell_row_pixet1strct)!=0){
####    pixet_1_struct[2,,mincell_row_pixet1strct][[1]] <- unique(append(pixet_1_struct[2,,mincell_row_pixet1strct][[1]], stats_WS_superpixel[3,sp,][[1]]));
####}

#pixet_1_struct[[min_cell_row]] <- unique(rbind(pixet_1_struct[[min_cell_row]],WSSr_pcoord[[sp]]+1,pixet_1_struct_cp[[min_cell_row]]));#pixelcoord   ..... this is something I added to tide over python swallowing pixel coords :-(

pixet_1_struct[[min_cell_row]] <- unique(rbind(pixet_1_struct[[min_cell_row]],WSSr_pcoord[[sp]]+1));

#remove that pixel entry
#stats_WS_superpixel[1,sp,][[1]] <- 0; #definitely not wise!
WSSr_area[sp] <- 0




#updating to rr_ data structures is just the intermediate..probably not correct to visualise since you are appending the new entries from superpixet2 to original cells. The original cells will have segments in pixet2 waiting to be reassigned. What you ideally need is to add these rr data structures to pixet1 CM and cluster that from iter_pixets >=2










