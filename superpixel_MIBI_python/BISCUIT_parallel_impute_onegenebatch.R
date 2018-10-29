## 24th May 2017
## Imputing in the case of a single gene split
##
## Code author SP
##################


X_std_all <- X_all;

#commented it for MIBI model on Feb 8th 2018

z_inferred_final <- factor(results.all.MCMC[[num_gene_batches]][[num_gene_sub_batches]][[1]]);
final_num_K <- length(unique(z_inferred_final));

#mean_alpha_inferred <- mean(alpha_inferred_final)
#mean_beta_inferred <- mean(beta_inferred_final)
uq_z <- unique(factor(results.all.MCMC[[num_gene_batches]][[num_gene_sub_batches]][[1]]));


final_num_K <- length(unique(uq_z));






































