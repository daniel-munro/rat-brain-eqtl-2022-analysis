library(data.table)

setwd("/projects/ps-palmer/apurva/susieR")

expname="physiological"
filename=paste0("./intervals/",expname,"_susieR_intervals.csv")

qtls=read.csv(filename,header=T,stringsAsFactors = F)

gwas_susie_dir="/projects/ps-palmer/apurva/susieR/causal_sets/"

fastenloc_dir="/projects/ps-palmer/apurva/fastenloc/"


#reformat susieR input file

qtls$susieR_causal_sets<-"present"


for(i in 1:nrow(qtls)){
  chr=qtls$chr[i]
  assoc=qtls$assoc[i]
  start=qtls$susie_interval_start[i]
  stop=qtls$susie_interval_stop[i]
  trait=qtls$trait[i]
  pos=qtls$pos[i]
  
  outfile=paste0(chr,"_",pos,"_",trait)
  raw=read.csv(paste0(gwas_susie_dir,outfile,".csv"),header = T,stringsAsFactors = F)
  if(nrow(raw)>0){
    
    raw$signalID="Loc1"
    data=raw[,c("rs","signalID","variable_prob")]
    
    write.table(data,paste0(fastenloc_dir,"/gwas_input/",outfile,".txt"),row.names=F,quote=F,col.names=F)
    system(paste0("gzip ",fastenloc_dir,"/gwas_input/",outfile,".txt"))
  }else{
    qtls$susieR_causal_sets[i]<-"absent"
  }

}

outdir<-fastenloc_dir
for(i in 1:nrow(qtls)){
  
  if(qtls$susieR_causal_sets[i]=="present"){
    
 
  chr=qtls$chr[i]
  assoc=qtls$assoc[i]
  start=qtls$susie_interval_start[i]
  stop=qtls$susie_interval_stop[i]
  trait=qtls$trait[i]
  pos=qtls$pos[i]
  
  outfile=paste0(chr,"_",pos,"_",trait)
  
  
  
  f<-paste0(outdir,"code/",trait,"_",chr,"_",pos,".sh")
  if (file.exists(f)) file.remove(f)
  d1<-paste0("#!/bin/bash")
  d2<-paste0("#PBS -N ",trait,"_",chr)  ##This needs to be changed for each pheno file name
  d3<-paste0("#PBS -S /bin/bash")
  d4<-paste0("#PBS -l walltime=3:00:00")
  d5<-paste0("#PBS -l nodes=1:ppn=5")
  d6<-paste0("#PBS -j oe")
  d7<-paste0("#PBS -o ",outdir,"pbs_log/",trait,"_",chr,"_",pos,".out")
  d8<-paste0("#PBS -q condo")
  

  bash1<- paste0("for tissue in Acbc LHB VoLo IL PL")
  bash2<-paste0("do")
  part1<-paste0("/home/aschitre/fastenloc/fastenloc/src/fastenloc.static -eqtl /home/dmunro/coloc/tensorqtl/fastenloc.eqtl.anno.$tissue.vcf.gz -gwas ")
  part2 <- paste0(outdir,"/gwas_input/",outfile,".txt.gz -t $tissue -thread 10 -prefix ",outdir,"output/",outfile,".$tissue")
  

  
  command <- paste0(part1,part2)
bash3<-"done"
  line=c(d1,d2,d3,d4,d5,d6,d7,d8,bash1,bash2,command,bash3)
  write(line,file=paste0(outdir,"code/",trait,"_",chr,"_",pos,".sh"),append=TRUE)
  
}
}


  
for tissue in Acbc LHB VoLo IL PL
do
/home/aschitre/fastenloc/fastenloc/src/fastenloc.static -eqtl /home/dmunro/coloc/tensorqtl/fastenloc.eqtl.anno.$tissue.vcf.gz -gwas /projects/ps-palmer/apurva/fastenloc/gwas_input/susie_chr16_mean_iop.txt.gz -t $tissue -thread 10 -prefix /projects/ps-palmer/apurva/fastenloc/output/chr16_mean_iop.$tissue
done
}