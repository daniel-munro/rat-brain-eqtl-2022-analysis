 
args <- commandArgs(trailingOnly = TRUE)
#print(args)
dir <- args[1]
trait <- args[2]
chr <- args[3]
pos <- args[4]

dir<-paste0(dir,"/")

library(data.table)
#Calculate beta

beta_filename=paste0(dir,"betas/",chr,"_",pos,"_",trait,"_beta.txt")
#Calculate correlation matrix (LD matrix) from dosages

dosage_filename=paste0(dir,"dosages/",chr,"_",pos,"_",trait,".dosage")

  
data=fread(dosage_filename,header=F,stringsAsFactors=F)
LD_matrix=data[,c(-1:-3)]
t_LD_matrix=t(LD_matrix)
transposed_corr=cor(t_LD_matrix)





library(susieR)
z = read.csv(beta_filename, header=T, stringsAsFactors=F)
rs = z$SNP




ld=transposed_corr

ld <- data.matrix(ld)


fitted_rss <- susie_rss(z$z, ld, L = 10,
                        estimate_residual_variance = TRUE, 
                        estimate_prior_variance = TRUE)


outfile=paste0(chr,"_",pos,"_",trait)
save(fitted_rss, file=paste0(dir,"fitted_rss/",outfile,".RData"))



variants = c()
for( set in names(fitted_rss$sets$cs)){
  vars = fitted_rss$sets$cs[set] 
  variants <- c(variants, vars)
}
variants = unlist(variants)



summary(fitted_rss)$vars
CS_table<-summary(fitted_rss)$vars
CS_table[which(CS_table$variable %in% variants),]
rs[variants]
fitted_rss$sets$cs


#extract the variants in the causal sets
CS_table<-summary(fitted_rss)$vars
final=CS_table[which(CS_table$variable %in% variants),]

CS_set=cbind(rs=rs[final$variable],final)


write.csv(CS_set,file=paste0(dir,"causal_sets/",outfile,".csv"),row.names=F,quote=F)

