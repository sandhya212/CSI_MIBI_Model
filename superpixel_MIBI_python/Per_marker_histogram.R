## 3rd Jan 2018
## Per marker histogram per cluster
##
## Code author SP
##
## 14th Feb 2018
## commenting out the x11()
###
###

############## helper functions ##############BEFORE####
# plot the per marker histograms per cluster

for (m in 1:choose_genes){
    if (all(full.data.3[,m]==0)==TRUE){
        full.data.3[,m]= 0.001;
    }
}

#f_rows <- 7#as.integer(choose_genes/2)


for (i in 1:final_num_K){
    print(paste("numbers of cells in cluster ",uq_z[i],"are ",length(which(z_inferred_final==uq_z[i]))));
    cells_in_cluster <- which(z_inferred_final==uq_z[i]);
    #x11(); par(mfrow=c(2,3));
    #title(main = paste("Cluster ",uq_z[i]))
    #for (m in 1:length(gene_names)){
    #    hist(full.data.3[cells_in_cluster,m],main = gene_names[m],breaks=50,xlab=paste0("cells in k=",i),col='lightgreen');
        
        # }
    
    #plot as pdf
    
    for ( j in 1:4){
        f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/patient1_M",j,"cluster",i,".pdf")
        pdf(file=f);
        par(mfrow=c(2,3));
        for (m in ((j-1)*6+1):(j*6)){
            if(m!=24){
                hist(full.data.3[cells_in_cluster,m],main = gene_names[m],breaks=50,xlab=paste0("cells in k=",i),col='lightgreen');
            }
        }
        dev.off()
    }

}



# plot the per marker histograms across all cells

#x11(); par(mfrow=c(2,3));
#for (m in 1:length(gene_names)){
#    hist(full.data.3[,m],main = gene_names[m],breaks=50,xlab='all cells',col='lightblue');
#}





for ( j in 1:4){
    f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/patient1_marker_all_",j,".pdf")
    pdf(file=f)
    par(mfrow=c(2,3));
    for (m in ((j-1)*6+1):(j*6)){
        if(m!=24){
            hist(full.data.3[,m],main = gene_names[m],breaks=50,xlab='all cells',col='lightblue');
        }
        
    }
    
    dev.off()
}













# plot per marker distribution on t-sne of X

scalar1 <- function(x) {x / sqrt(sum(x^2))} # to normalise the vector
#ygPal <- colorRampPalette(c('yellow','red'))
ygPal <- colorRampPalette(c('yellow','blue'))
cols <- brewer.pal(7, "Spectral")
cols3 <- colorRampPalette(brewer.pal(10, "Spectral"))(256)
#This adds a column of color values
# based on the y values

#X_tsne <- Rtsne(X_1,check_duplicates = FALSE,perplexity=30);


#x11(); par(mfrow=c(2,3));
#for (m in 1:6){
#    print(m)
#    d1 <- scalar1(full.data.3[,m]);
#    c1 <- ygPal(10)[as.numeric(cut(d1,breaks = 400))]
#    plot(X_tsne$Y[,1],X_tsne$Y[,2],col = c1,  main=gene_names[m]);

#}

#redblue<-colorRampPalette(c("red","orange","blue"))
#filled.contour(m,color.palette=redblue)


for ( j in 1:4){
    f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/patient1_per_marker_before_1_",j,".pdf")
    pdf(file=f)
    par(mfrow=c(3,2));
    for (m in ((j-1)*6+1):(j*6)){
        if(m!=24){
            d1 <- scalar1(full.data.3[,m]);
            print(m)
            #print(max(d1));
            #print(min(d1));
            #print(m);
            c1 <- ygPal(10)[as.numeric(cut(d1,breaks = 40))]
            #c1 <- 'jet'
            #c1 =cols[as.numeric(cut(d1,breaks = 40))]
            plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = c1,pch=16,  main=gene_names[m]);
        }
    
    }

    dev.off()
}

#x11();
#plot(X_tsne$Y[,1],X_tsne$Y[,2],col = col_palette[1*1*(z_inferred_final_plot)],  main="Final_plot");
f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/Final_plot.pdf")
pdf(file=f)
plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],pch=16, col = col_palette[1*1*(z_inferred_final_plot)],  main="Final_plot");
dev.off()

##
#x11(); pie(rep(1,length(uq_z)), col=(col_vector[c(uq_z)]))
f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/final_pie.pdf")
pdf(file=f)
pie(rep(1,length(uq_z)), col=(col_vector[c(uq_z)]))
dev.off()



f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/count_matrix.txt");
write.matrix(full.data.3,file=f,sep="\t")
f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/count_matrix_revised.mat");
writeMat(f,count_matrix_revised = full.data.3);


#heatmap(t(full.data.3[order(z_inferred_final),])%*%full.data.3[order(z_inferred_final),],col=col)
#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/LP_count_matrix.pdf");
#pdf(file=f)
#levelplot(full.data.3[order(z_inferred_final),]%*%t(full.data.3[order(z_inferred_final),]),col=cols3, main='Levelplot of Cov(CM)')
#dev.off()

#f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/HM_count_matrix.pdf");
#pdf(file=f)
#heatmap(cov(t(full.data.3[order(z_inferred_final),])),col=cols3, main='Heatmap of Cov(CM)')
#dev.off()



##commented early on
#par(mfrow=c(1,2))
#col <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
#heatmap(t(full.data.3[order(z_inferred_final),])%*%full.data.3[order(z_inferred_final),],col=col)
#levelplot(full.data.3[order(z_inferred_final),]%*%t(full.data.3[order(z_inferred_final),]),col=col)
#heatmap(full.data.3[order(z_inferred_final),]%*%t(full.data.3[order(z_inferred_final),]),col=col)

#cols2 <- colorRampPalette(c("blue","white","red"))(256)
#levelplot(full.data.3[order(z_inferred_final),]%*%t(full.data.3[order(z_inferred_final),]),col=cols2)
#heatmap(full.data.3[order(z_inferred_final),]%*%t(full.data.3[order(z_inferred_final),]),col=cols2)





