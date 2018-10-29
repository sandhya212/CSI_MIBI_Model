## 28th Feb 2018
## MIBI model main and helper functions
## main()
## Code author SP
##

###
###


############## helper functions ##############
#gene_names <- c('CD20','Pan-Keratin','CD45','HLA-DR','CD3', 'CD68');
#gene_names_1 <- c('CD20','CD45','HLA-DR','CD3', 'CD68');


#gene_names = c('CD20','Pan-Keratin','CD45','HLA-DR','CD3','CD68','FoxP3','CD4','CD8','CD16','MPO','CD56','CD209','CD11c','CD11b','Keratin6','Keratin17','p53','Beta catenin','EGFR','Vimentin','SMA','CD31')

gene_names_1 = c('CD20','CD45','HLA-DR','CD3','CD68','FoxP3','CD4','CD8','CD16','MPO','CD56','CD209','CD11c','CD11b')

f_rows <- as.integer(length(gene_names_1)/2)

choose_genes <- length(gene_names)

scalar1 <- function(x) {x / sqrt(sum(x^2))}

my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
nbins <- 25

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

#####
##For Leeat data
#####

###Plot heatmap for covariance
ygPal <- colorRampPalette(c('yellow','blue'))
cols <- brewer.pal(7, "Spectral")
cols3 <- colorRampPalette(brewer.pal(10, "Spectral"))(256)

cov_plot <- readMat(paste0(working_path,'/output_6_markers_Leeat_original_py/plots/extras/',"count_matrix_revised.mat"));
cov_plot <- cov_plot[[1]];
scaled_A <- t(apply(cov_plot, MARGIN=1, FUN=scale))
scaled_A[is.nan(scaled_A)] <- 0

colnames(scaled_A) <- gene_names

z_inf <- read.csv(paste0(working_path,'/output_6_markers_Leeat_original_py/plots/extras/',"cluster_probabilities.csv"), header=TRUE)
z_inf$z_inferred



f <- paste0(paste0(working_path,'/output_6_markers_Leeat_original_py/plots/extras/HM_count_matrix_scaled.pdf'));
pdf(file=f)
heatmap(cov(t(scaled_A[order(z_inf$z_inferred),])),col=cols3, main='Heatmap of Cov(CM)')
dev.off()


f <- paste0(paste0(working_path,'/output_6_markers_Leeat_original_py/plots/extras/cov_scaled.pdf'));
pdf(file=f)
image(cov(t(scaled_A[order(z_inf$z_inferred),])))
dev.off()

#f <- paste0(paste0(working_path,'/output_6_markers_Leeat_original_py/plots/extras/cov_unscaled.pdf'));
#pdf(file=f)
#image(cov(t(cov_plot[order(z_inf$z_inferred),])))
#dev.off()

###Plot t-SNE with marker distribution
#X_tsne <- Rtsne(cov_plot,check_duplicates = FALSE,perplexity=30);
#read in Leeat's tsne coordinates
X_tsne <- read.csv(paste0(working_path,'/output_6_markers_Leeat_original_py/plots/extras/',"X_tsne_all.csv"), header=TRUE)
X_tsne_x <- X_tsne[,2]
X_tsne_y <- X_tsne[,3]

#f <- paste0(working_path,'/output_6_markers_Leeat_original_py/plots/extras/patient1_per_marker_before_1.pdf')
#pdf(file=f)

#par(mfrow=c(2,f_rows));
#for (m in 1:choose_genes){
#    d1 <- scalar1(cov_plot[,m]);
#    print(max(d1));
#    print(min(d1));
#    print(m);
#    c1 <- ygPal(10)[as.numeric(cut(d1,breaks = 40))]
    #c1 <- 'jet'
    #c1 =cols[as.numeric(cut(d1,breaks = 40))]
#    plot(X_tsne$Y[,1],X_tsne$Y[,2],col = c1,pch=16,  main=gene_names[m]);
#
#}
#dev.off()

##ld <- cov(t(scaled_A[order(z_inf$z_inferred),]))
##heatmap(cov(ld),col=cols3, main='Heatmap of Cov(CM)')

###plot the marker distribution across all channels for Leeat
#f <- paste0(getwd(),'/output_6_markers_Leeat_original_py/plots/extras/patient1_marker_all_1.pdf')
#pdf(file=f)
#par(mfrow=c(2,f_rows));
#for (m in 1:length(gene_names)){
#    hist(cov_plot[,m],main = gene_names[m],breaks=50,xlab='all cells',col='lightblue');
#}
#dev.copy(p,'myplot.png')
#dev.off()

#f <- paste0(getwd(),'/output_6_markers_Leeat_original_py/plots/extras/patient1_marker_all_scaled.pdf')
#pdf(file=f)
#par(mfrow=c(2,f_rows));
#for (m in 1:length(gene_names)){
#    hist(scaled_A[,m],main = gene_names[m],breaks=50,xlab='all cells',col='lightblue');
#}
#dev.copy(p,'myplot.png')
#dev.off()


for ( j in 1:4){
    f <- paste0(getwd(),'/output_6_markers_Leeat_original_py/plots/extras/patient1_marker_all_scaled_',j,'.pdf')
    pdf(file=f)
    par(mfrow=c(3,2));
    for (m in ((j-1)*6+1):(j*6)){
        if(m!=24){
            hist(scaled_A[,m],main = gene_names[m],breaks=50,xlab='all cells',col='lightblue');
        }
        
    }
    
    dev.off()
}





for ( j in 1:2){
    f <- paste0(getwd(),'/output_6_markers_Leeat_original_py/plots/extras/biaxial_scaled_',j,'.pdf')
    pdf(file=f)
    par(mfrow=c(2,4));
    for (m in ((j-1)*8+1):(j*8)){
        if(m<15){
            lab1 = paste0(gene_names_1[m])
            ind = which(colnames(scaled_A)==gene_names_1[m])
            plot(scaled_A[,2],scaled_A[,ind], xlab = 'Pan_CK', ylab = lab1, pch = 16, frame = FALSE, col = 'blue')
            
        }
        
    }
    
    dev.off()
}


for ( j in 1:2){
    f <- paste0(getwd(),'/output_6_markers_Leeat_original_py/plots/extras/biaxial_',j,'.pdf')
    pdf(file=f)
    par(mfrow=c(2,4));
    for (m in ((j-1)*8+1):(j*8)){
        if(m<15){
            lab1 = paste0(gene_names_1[m])
            ind = which(colnames(scaled_A)==gene_names_1[m])
            plot(cov_plot[,2],cov_plot[,ind], xlab = 'Pan_CK', ylab = lab1, pch = 16, frame = FALSE, col = 'blue')
            
        }
        
    }
    
    dev.off()
}




#f <- paste0(getwd(),'/output_6_markers_Leeat_original_py/plots/extras/biaxial_marginals.pdf')
#pdf(file=f)

#par(mfrow=c(2,f_rows));
#for (m in 1:length(gene_names_1)){
#    lab1 = paste0(gene_names_1[m])
#    ind = which(colnames(scaled_A)==gene_names_1[m])
#    scatterhist(cov_plot[,2],cov_plot[,ind], xlab = "Pan_CK", ylab = lab1, frame = FALSE)
#
#}
#dev.off()

#######
#######



for ( j in 1:2){
    f <- paste0(getwd(),'/output_6_markers_Leeat_original_py/plots/extras/kde_2D_',j,'.pdf')
    pdf(file=f)
    df=cov_plot
    par(mfrow=c(2,4));
    for (m in ((j-1)*8+1):(j*8)){
        if(m<15){
            lab1 = paste0(gene_names_1[m])
            ind = which(colnames(scaled_A)==gene_names_1[m])
    
        x.bin <- seq(floor(min(df[,2])), ceiling(max(df[,2])), length=nbins)
        y.bin <- seq(floor(min(df[,ind])), ceiling(max(df[,ind])), length=nbins)
    
        freq <-  as.data.frame(table(findInterval(df[,2], x.bin),findInterval(df[,ind], y.bin)))
        freq[,1] <- as.numeric(freq[,1])
        freq[,2] <- as.numeric(freq[,2])
    
        freq2D <- diag(nbins)*0
        freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
    
    
        image(x.bin, y.bin, log(freq2D), col=r,xlab = 'Pan_CK', ylab = lab1,pch = 16)
        }
    }
    dev.off()
}



for ( j in 1:2){
    f <- paste0(getwd(),'/output_6_markers_Leeat_original_py/plots/extras/kde_2D_1_',j,'.pdf')
    pdf(file=f)
    df=cov_plot
    par(mfrow=c(2,4));
    for (m in ((j-1)*8+1):(j*8)){
        if(m<15){
            
            lab1 = paste0(gene_names_1[m])
            ind = which(colnames(scaled_A)==gene_names_1[m])
            
            x1 <- cov_plot[,2]
            x2 <- cov_plot[,ind]
            df <- data.frame(x1,x2)
            
            ## Use densCols() output to get density at each point
            x <- densCols(x1,x2, colramp=colorRampPalette(c('black', 'white')))
            df$dens <- col2rgb(x)[1,] + 1L
            
            ## Map densities to colors
            cols <-  colorRampPalette(c('#000099', '#00FEFF', '#45FE4F',
            '#FCFF00', '#FF9400', '#FF3100'))(256)
            df$col <- cols[df$dens]
            
            ## Plot it, reordering rows so that densest points are plotted on top
            plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2,xlab = 'Pan_CK', ylab = lab1)
        }
    }
    dev.off()
}





#####
##For iterations
#####

#for (it in 1:(iter_pixets-1)){ #max_MCMC){
for (it in 1:mcmc_counter){ #max_MCMC){
    cov_plot <- readMat(paste0(working_path,'/output_6_markers_updated_1_',it,'/plots/extras/',"count_matrix_revised.mat"));
    cov_plot <- cov_plot[[1]];
    scaled_A <- t(apply(cov_plot, MARGIN=1, FUN=scale))
    scaled_A[is.nan(scaled_A)] <- 0
    
    colnames(scaled_A) <- gene_names

    z_inf <- read.csv(paste0(working_path,'/output_6_markers_updated_1_',it,'/plots/extras/',"cluster_probabilities.csv"), header=TRUE)
    z_inf$z_inferred
    
    ld <- cov(t(scaled_A[order(z_inf$z_inferred),]))
    
    if (it==mcmc_counter){
        f <- paste0(paste0(working_path,'/output_6_markers_updated_1_',it,'/plots/extras/HM_count_matrix_scaled.pdf'));
        pdf(file=f)
        heatmap(cov(t(scaled_A[order(z_inf$z_inferred),])),col=cols3, main='Heatmap of Cov(CM)')
        #heatmap(cov(ld[order(z_inf$z_inferred),]),col=cols3, main='Heatmap of Cov(CM)')
        dev.off()
    
    }
    
    f <- paste0(paste0(working_path,'/output_6_markers_updated_1_',it,'/plots/extras/cov_scaled.pdf'));
    pdf(file=f)
    image(cov(t(scaled_A[order(z_inf$z_inferred),])))
    dev.off()
    
    #f <- paste0(paste0(working_path,'/output_6_markers_updated_1_',it,'/plots/extras/cov_unscaled.pdf'));
    #pdf(file=f)
    #image(cov(t(cov_plot[order(z_inf$z_inferred),])))
    #dev.off()
    
    
    ###Plot t-SNE with marker distribution
    #X_tsne <- Rtsne(cov_plot,check_duplicates = FALSE,perplexity=30);
    

    
    
    
    
    for ( j in 1:4){
        f <- paste0(working_path,'/output_6_markers_updated_1_',it,'/plots/extras/patient1_per_marker_before_1_Leeat_coord_',j,'.pdf')
        pdf(file=f)
        par(mfrow=c(3,2));
        for (m in ((j-1)*6+1):(j*6)){
            if(m!=24){
                d1 <- scalar1(cov_plot[,m]);
                print(m)
                #print(max(d1));
                #print(min(d1));
                #print(m);
                c1 <- ygPal(10)[as.numeric(cut(d1,breaks = 40))]
                #c1 <- 'jet'
                #c1 =cols[as.numeric(cut(d1,breaks = 40))]
                plot(X_tsne_x,X_tsne_y,col = c1,pch=16,  main=gene_names[m]);
            }
            
        }
        
        dev.off()
    }

    
    
    
    
    
    ###plot the marker distribution across all channels
    #f <- paste0(getwd(),'/output_6_markers_updated_1_',it,"/plots/extras/patient1_marker_all_1.pdf")
    #pdf(file=f)
    #par(mfrow=c(2,f_rows));
    #for (m in 1:length(gene_names)){
    #    hist(cov_plot[,m],main = gene_names[m],breaks=50,xlab='all cells',col='lightblue');
    #}
    
    #dev.off()
    
    for ( j in 1:4){
        f <- paste0(getwd(),'/output_6_markers_updated_1_',it,'/plots/extras/patient1_marker_all_scaled_',j,'.pdf')
        pdf(file=f)
        par(mfrow=c(3,2));
        for (m in ((j-1)*6+1):(j*6)){
            if(m!=24){
                hist(scaled_A[,m],main = gene_names[m],breaks=50,xlab='all cells',col='lightblue');
            }
        
        }
    
        dev.off()
    }


    for ( j in 1:2){
        f <- paste0(getwd(),'/output_6_markers_updated_1_',it,'/plots/extras/biaxial_scaled_',j,'.pdf')
        pdf(file=f)
        par(mfrow=c(2,4));
        for (m in ((j-1)*8+1):(j*8)){
            if(m<15){
                lab1 = paste0(gene_names_1[m])
                ind = which(colnames(scaled_A)==gene_names_1[m])
                plot(scaled_A[,2],scaled_A[,ind], xlab = "Pan_CK", ylab = lab1, pch = 16, frame = FALSE, col = 'blue')
            
            }
        
        }
    
        dev.off()
    }



    for ( j in 1:2){
        f <- paste0(getwd(),'/output_6_markers_updated_1_',it,'/plots/extras/biaxial_',j,'.pdf')
        pdf(file=f)
        par(mfrow=c(2,4));
        for (m in ((j-1)*8+1):(j*8)){
            if(m<15){
                lab1 = paste0(gene_names_1[m])
                ind = which(colnames(scaled_A)==gene_names_1[m])
                plot(cov_plot[,2],cov_plot[,ind], xlab = "Pan_CK", ylab = lab1, pch = 16, frame = FALSE, col = 'blue')
            
            }
        
        }
    
        dev.off()
    }




    #f <- paste0(getwd(),'/output_6_markers_updated_1_',it,"/plots/extras/biaxial_marginals.pdf")
    #pdf(file=f)
    
    #par(mfrow=c(2,f_rows));
    #for (m in 1:length(gene_names_1)){
    #    lab1 = paste0(gene_names_1[m])
    #    ind = which(colnames(scaled_A)==gene_names_1[m])
    #    scatterhist(cov_plot[,2],cov_plot[,ind], xlab = "Pan_CK", ylab = lab1, frame = FALSE)
    #
    #}
    #dev.off()
    
    #######

    
 
    
    
    for ( j in 1:2){
        f <- paste0(getwd(),'/output_6_markers_updated_1_',it,'/plots/extras/kde_2D_',j,'.pdf')
        pdf(file=f)
        df=cov_plot
        par(mfrow=c(2,4));
        for (m in ((j-1)*8+1):(j*8)){
            if(m<15){
                lab1 = paste0(gene_names_1[m])
                ind = which(colnames(scaled_A)==gene_names_1[m])
                
                x.bin <- seq(floor(min(df[,2])), ceiling(max(df[,2])), length=nbins)
                y.bin <- seq(floor(min(df[,ind])), ceiling(max(df[,ind])), length=nbins)
                
                freq <-  as.data.frame(table(findInterval(df[,2], x.bin),findInterval(df[,ind], y.bin)))
                freq[,1] <- as.numeric(freq[,1])
                freq[,2] <- as.numeric(freq[,2])
                
                freq2D <- diag(nbins)*0
                freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
                
                
                image(x.bin, y.bin, log(freq2D), col=r,xlab = 'Pan_CK', ylab = lab1,pch = 16)
            }
        }
        dev.off()
    }
    
    
    
    for ( j in 1:2){
        f <- paste0(getwd(),'/output_6_markers_updated_1_',it,'/plots/extras/kde_2D_1_',j,'.pdf')
        pdf(file=f)
        df=cov_plot
        par(mfrow=c(2,4));
        for (m in ((j-1)*8+1):(j*8)){
            if(m<15){
                
                lab1 = paste0(gene_names_1[m])
                ind = which(colnames(scaled_A)==gene_names_1[m])
                
                x1 <- cov_plot[,2]
                x2 <- cov_plot[,ind]
                df <- data.frame(x1,x2)
                
                ## Use densCols() output to get density at each point
                x <- densCols(x1,x2, colramp=colorRampPalette(c('black', 'white')))
                df$dens <- col2rgb(x)[1,] + 1L
                
                ## Map densities to colors
                cols <-  colorRampPalette(c('#000099', '#00FEFF', '#45FE4F',
                '#FCFF00', '#FF9400', '#FF3100'))(256)
                df$col <- cols[df$dens]
                
                ## Plot it, reordering rows so that densest points are plotted on top
                plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2,xlab = 'Pan_CK', ylab = lab1)
            }
        }
        dev.off()
    }
    

    
    
    
    
    
    
}


