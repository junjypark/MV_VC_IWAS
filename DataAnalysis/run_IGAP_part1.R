library(foreach)

chr=19

n.genes=length(list.files(paste0("/gpfs/fs0/project/j/junpark/tianyu47/LD_0.5/chr",chr)))

#Extract g2i
# setwd(paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_UKBGWAS_0.5/chr",chr))
setwd(paste0("/gpfs/fs0/project/j/junpark/tianyu47/UKBGWAS_0.5/chr",chr))
lst=list.files()
lst_split=strsplit(lst, "_")
imgs=do.call("c",lapply(lst_split, function(x)(x[4])))
genes=do.call("c",lapply(lst_split, function(x)(x[5])))

for (gene in 1:n.genes){
  tryCatch({
  print(gene)
  index=which(genes==paste0("gene",gene,".txt"))
  IDP_index=imgs[index]
  IDP_index2=as.numeric(substr(IDP_index,4,7))
  
  result=foreach(j=1:length(IDP_index2),.combine="rbind")%do%{
    file=paste0("UKBGWAS_0.5_chr",chr,"_",IDP_index[j],"_gene", gene,".txt")
    
    result=read.table(file, header=T)
    tstat=result$beta/result$se
    result$gti=sign(tstat)*sqrt((tstat^2/40000)/(tstat^2/(40000)+1))
    result$idpid=IDP_index2[j]
    result  
  }
  
  LD=as.matrix(readRDS(paste0("/gpfs/fs0/project/j/junpark/tianyu47/LD_0.5/chr", chr, "/LD_chr",chr,"_gene",gene,".rds")))
  
  mat=match(result$rsID, rownames(LD))
  combine=list(result=result, LD=LD[mat,mat])
  
  saveRDS(combine, paste0("/gpfs/fs1/home/j/junpark/junpark/imgGenet/chr",chr,"/combine_gene", gene,".rds"))
  }, error=function(e){})
}







##
library(MASS)
library(data.table)
library(Matrix)
for (chr in 16:18){
  n.genes=length(list.files(paste0("/gpfs/fs0/project/j/junpark/tianyu47/LD_0.5/chr",chr)))
  print(paste0("chr=",chr))
  for (gene in 1:n.genes){
    tryCatch({
      print(gene)
      setwd(paste0("/gpfs/fs0/project/j/junpark/tianyu47/intersection_snplist_0.5_bygene/chr",chr))
      # setwd(paste0("/gpfs/fs0/scratch/j/junpark/tianyu47/UKBB/UKBB_IGAP_1000G_snplist_0.5_bygene/chr",chr))
      tbl=fread(paste0("SNPresults_0.5_chr",chr,"_gene",gene,".txt"))
      setwd(paste0("/gpfs/fs1/home/j/junpark/junpark/imgGenet/chr",chr))
      result=readRDS(paste0("combine_gene",gene,".rds"))
      unique.snp=unique(result$result[,1])
      n.unique.snp=length(unique.snp)
      n.idp.s=sum(tbl$IDP_type=="S")
      n.idp.d=sum(tbl$IDP_type=="D")
      n.idp.f=sum(tbl$IDP_type=="F")
      id.idp.s=tbl$IDP_ID[tbl$IDP_type=="S"]
      id.idp.d=tbl$IDP_ID[tbl$IDP_type=="D"]
      id.idp.f=tbl$IDP_ID[tbl$IDP_type=="F"]
      coef.mat=list(S=matrix(0,n.unique.snp,n.idp.s),
                    D=matrix(0,n.unique.snp,n.idp.d),
                    F=matrix(0,n.unique.snp,n.idp.f))
      rownames(coef.mat$S)=rownames(coef.mat$D)=rownames(coef.mat$F)=unique.snp
      colnames(coef.mat$S)=id.idp.s
      colnames(coef.mat$D)=id.idp.d
      colnames(coef.mat$F)=id.idp.f
      
      for (j in 1:n.idp.s){
        ind=which(result$result$idpid==id.idp.s[j] )
        if (length(ind)>0){
          coef=result$result$gti[ind]
          corr=result$LD[ind,ind]
          snpid=result$result$rsID[ind]
          coef.mat$S[snpid,j]=as.numeric(ginv(corr)%*%coef)
        }
      }
      
      for (j in 1:n.idp.d){
        ind=which(result$result$idpid==id.idp.d[j] )
        if (length(ind)>0){
          coef=result$result$gti[ind]
          corr=result$LD[ind,ind]
          snpid=result$result$rsID[ind]
          coef.mat$D[snpid,j]=as.numeric(ginv(corr)%*%coef)
        }
      }
      
      for (j in 1:n.idp.f){
        ind=which(result$result$idpid==id.idp.f[j] )
        if (length(ind)>0){
          coef=result$result$gti[ind]
          corr=result$LD[ind,ind]
          snpid=result$result$rsID[ind]
          coef.mat$F[snpid,j]=as.numeric(ginv(corr)%*%coef)
        }
      }
      coef.mat$S=Matrix(coef.mat$S)
      coef.mat$D=Matrix(coef.mat$D)
      coef.mat$F=Matrix(coef.mat$F)
      
      setwd(paste0("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/imgGenet_coefMat/chr",chr))
      saveRDS(coef.mat, paste0("gene",gene,".rds"))
    }, error=function(e){})
  }
  
}
