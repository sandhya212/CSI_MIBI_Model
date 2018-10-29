## 21st Dec 2016
## BISCUIT preparing input data
## (example code for Zeisel et al mouse cortex data
## Downloaded from : https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt
##
## Code author SP

## 28th Dec 2016
## choose genes if needed

## 14th March 2017
## pre-processing data
##

## 18th April 2017
## added logic to split genes based on co-variances rather than variance alone
##
##
##############################################################################################################
#### Modify this portion to accommodate your input data.                                                  ####
##############################################################################################################


print("Loading Data")

########

    
    
# all for pixet_1 based count matrix - coming from Matlab. Do just when iter_pixets==1

if (iter_pixets==1){
    #full.data.0 <- readMat(input_file_name);
    full.data.0 <- read.csv(input_file_name, header=FALSE);
    
    typeof(full.data.0) # to see what type the .mat file is
    str(full.data.0) # to see the $ entries
    
    #full.data <- full.data.0$dataCells;
    #full.data <- full.data.0[[1]];
    full.data <- full.data.0;
    
    dim(full.data)
    #gene_names <- colnames(full.data); #gene_names
    #gene_names <- c('Vimentin','SMA','FoxP3','Lag3','CD4','CD16','CD56','CD31','EGFR','CD209','CD11c','CD138','CD68','CD8','CD3','Keratin17','CD63','CD20','p53','Beta catenin','HLA-DR','CD11b','CD45','Pan-Keratin','MPO','Keratin6');
    
    #gene_names <- c('CD20','Pan-Keratin','CD45','HLA-DR','CD3', 'CD68');
    
    
    #gene_names = c('CD20','Pan-Keratin','CD45','HLA-DR','CD3','CD68','FoxP3','CD4','CD8','CD16','MPO','CD56','CD209','CD11c','CD11b','Keratin6','Keratin17','p53','Beta catenin','EGFR','Vimentin','SMA','CD31')
    
    
    rownames(full.data); # #cell_names
    choose_cells <- nrow(full.data)
    #creating the cellsXgenes data
    full.data.1 <- as.matrix(full.data);

}else{
    #construct a count matrix from the newly updated pixet_1
    
    data_updated <- matrix(0,choose_cells,choose_genes);
    
    #for(si in 1:choose_cells){ #building the revised count matrix from image_tensor using pixet_1_struct with revisions from last run.(those that fell in the same cell or same neighbouring cell. These were removed from pixet_2 and moved to pixet_1)
    #print(paste('s is ', si));
    #sim_pixels_list <- which(stats_revised[3,si,1][[1]] %in% pixet_1[,3]);
    #print(paste('intersection of pixels between stats_revised and pixet_1 are ', stats_revised[3,si,1][[1]][sim_pixels_list]));
        
        #if(length(sim_pixels_list)>0){
        #    currData_updated <- image_tensor[stats_revised[3,si,1][[1]][sim_pixels_list],]
        #    if (length(sim_pixels_list)==1){
        #         data_updated[si,] <- currData_updated;
        #         print(data_updated[si,])
        #     }else{
        #        data_updated[si,] <- colSums(currData_updated);
        #       print(data_updated[si,])
        #   }
        #}#else{###This part occurs when the new pixet_1 has no similar pixels in stats_revised. But stats_revised is the superset. So why is this? - OK, so this is because although stats_revised is all of our valid pixels, when we construct pixet1 from...check this
        #   data_updated[si,] = rep(-1000,numgenes)
        #   print(data_updated[si,])
             #}
             #}
            
    ## to tidy up the last run and make it spacially aware
    #if (iter_pixets==(max_MCMC+1)){ #uncomment this if you want to update the pixet1 at the last stage only
    #    for (ls in 1:choose_cells){
    #        pixet_1_struct[[ls]] <- unique(rbind(pixet_1_struct[[ls]],pixet_1_struct_cp[[ls]]));#pixelcoord
    #    }

    #}

    
    for (si in 1:choose_cells){
        #print(si)
        if(is.matrix(pixet_1_struct[[si]])){
            indi <- logical(length=dim(pixet_1_struct[[si]])[1]);
            #indi <- logical(nrow(pixet_1_struct[[si]]));
        }else{
            indi <- logical(length=length(pixet_1_struct[[si]]));
        }
        indi = !indi
        d1 <- which(pixet_1_struct[[si]][,1] == -9)
        indi[d1] <- FALSE
        
        pixie <- pixet_1_struct[[si]][indi,]
        n_l <- nrow(pixie);
        
        d_temp <- matrix(0,1,choose_genes)
        #if (length(n_l != 0)){
        #if (dim(pixie)[1] != 0){
        if (is.matrix(pixie)){
            if (dim(pixie)[1] != 0){ #we have entries
                for (n_si in 1: n_l){
                    #print(n_si)
                    #print(image_tensor[c(pixet_1_struct[[si]][n_si,1]),c(pixet_1_struct[[si]][n_si,2]),])
                    #d_temp <- d_temp + image_tensor[c(pixet_1_struct[[si]][n_si,1]),c(pixet_1_struct[[si]][n_si,2]),];
                    d_temp <- d_temp + image_tensor[pixie[n_si,1],pixie[n_si,2],];
                }
            }
        }else{ #just one entry with 2 columns, a vector actually
                d_temp <- d_temp + image_tensor[pixie[1],pixie[2],];
        }
        data_updated[si,] <- d_temp;
    }

        full.data.1 <- data_updated;
}





dim(full.data.1);
full.data.1[is.na(full.data.1)] <- 0;




#save(full.data.1,file="full.data.1.RData")

#Idea 1 to get meaningful genes: choose genes based on highest global co-expression
stddev.genes <- apply(full.data.1,2,sd);## find std dev of genes
f <- paste0(getwd(),"/",output_folder_name,"/plots/Stddev_genes.pdf")
pdf(file=f);
plot(sort(stddev.genes,decreasing=T),typ='o')
dev.off();
#full.data.2 <- full.data.1[,order(stddev.genes,decreasing=TRUE)];
#gene_names <- gene_names[order(stddev.genes,decreasing=T)]; 

#Idea 2 to get meaningful genes: perform disparity check between rows and columns
#emp.cov <- cov(full.data.1);
#diag(emp.cov) <- 0;
#rowsums.emp.cov <- rowSums(emp.cov);
#colsums.emp.cov <- colSums(emp.cov);
#gene.disparity <- rowsums.emp.cov + colsums.emp.cov
#f <- paste0(getwd(),"/",output_folder_name,"/plots/Disparity_genes.pdf")
#pdf(file=f);
#plot(sort(gene.disparity,decreasing=T),typ='o')
#dev.off();
#full.data.2 <- full.data.1[,order(gene.disparity,decreasing=TRUE)];
#gene_names <- gene_names[order(gene.disparity,decreasing=T)];

#Idea 3 to get meaningful genes: use Fiedler vector to split genes (graph theoretic partition)
#print("Calculating the Fiedler vector of the data")


#   L.mat <- diag(colsums.emp.cov) - emp.cov; #ensure L.mat is singular p.s.d
#    f.vec <- fiedler.vector(L.mat);
#    save(f.vec,file="fvec.RData");



#load(file="fvec.RData")
#full.data.2 <- full.data.1[,order(f.vec,decreasing=TRUE)]
#gene_names <- gene_names[order(f.vec,decreasing=TRUE)];

#choose cells

lib_size <- rowSums(full.data.1);

f <- paste0(getwd(),"/",output_folder_name,"/plots/lib_size_hist.pdf")
pdf(file=f);
par(mfrow=c(2,1))
plot(sort(lib_size,decreasing=T),typ='l')
hist(lib_size,breaks=500)
dev.off();

#if (input_file_name=="MERGEDtumors_subsetgenes_counts.csv"){
#    full.data.2 <- full.data.2[which(lib_size>400),]
#}

##log transform data

print("Ensuring entire data is numeric and then log transforming it")

#if(iter_pixets==1){
#    X_1 <- log(full.data.2+1);
#}else{
#    X_1 <- full.data.2;
#}
##

X_1 <- log(full.data.1+1);

if(exists("choose_cells")){
    numcells <- choose_cells;
}else{
    numcells <- dim(X_1)[1];
}
print(paste("numcells is", numcells))

if(exists("choose_genes")){
    numgenes <- choose_genes;
}else{
    numgenes <- dim(X_1)[2];
}
print(paste("numgenes is", numgenes))
###

##write the genes used in this run into a file
##$$$$
f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_means/Genes_selected.csv")
write.csv(gene_names[1:numgenes], file=f);
f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_Sigmas/Genes_selected.csv")
write.csv(gene_names[1:numgenes], file=f);


#############################################################################################
#############################################################################################
#############################################################################################

tot_numgenes <- dim(full.data.1)[2];
#rm(full.data.1)
num_gene_batches <- floor(numgenes/gene_batch);
print(paste0('Number of gene batches is ', num_gene_batches));

numgenes <- num_gene_batches * gene_batch;

num_gene_sub_batches <- sub_batch(num_gene_batches);
print(paste0('Number of gene subbatches is ', num_gene_sub_batches));


#############################################################################################
#############################################################################################
#############################################################################################
# preparing the dataset as per user-defined number of cells and genes


    
full.data.3 <- full.data.1[1:numcells,1:numgenes];

## Ensure data is numeric
print("Ensuring user-specified data is numeric")

X_all <- matrix(as.numeric(full.data.3),nrow=numcells,ncol=numgenes);
#if(iter_pixets==1){
    X_all <- log(X_all + 0.5); # log normalisation. + 0.1 to account for zero entries in X that cannot be log transformed.
#}#else{
#   X_all[which(rowSums(X_all)==0),] = rep(min(X_all)-100,numgenes); #keep this for now, but see why this is happening
#}

###
# Visualisation

## centering X and plotting
#X_c_all<- center_colmeans(X_all);
N <- numcells;
D <- numgenes;
n <- rep(0,N)
#for( i in 1:N){
#    n[i] <- norm_vec(X_c_all[i,])
#}
#X_c_norm_all <- X_c_all/max(n)

print('Computing t-sne projection of the data')
X_tsne_all <- Rtsne(X_all,check_duplicates = FALSE);


## plotting standardised X
#X_std_all <- project.data(X_all,D);


##Global normalised data

#log_lib_size <- rowSums(X_all);
#X_all_global_norm <- X_all/(log_lib_size + 0.0001);
#X_tsne_all_global_norm <- Rtsne(X_all_global_norm,check_duplicates = FALSE);

#rm(X_all)

