args=(commandArgs(TRUE))

job_num=as.numeric(args[1])

chr = job_num

print(chr)

source("/gpfs/fs0/scratch/j/junpark/junpark/MV_VC_IWAS/files/method.R")

chrinfo=read.csv(paste0("/gpfs/fs1/home/j/junpark/junpark/1000G_rsid/rsid_chr",chr,"_converted.csv"))
chrstat=readRDS(paste0("/gpfs/fs0/scratch/j/junpark/junpark/neale/ukb_neale_chr", chr, ".rds"))
chrstat_pos=do.call("c",lapply(strsplit(chrstat$variant, ":"), function(x){x[2]}))
chrstat_pos=as.numeric(chrstat_pos)
hg=readRDS(paste0("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/UKBB_ImgSampleGenotype/hg19info/hg19_chr",chr,".rds"))
pvecs=pvecs2=pvecs3=pvecs4=vector("list", nrow(hg))
for (gene in 1:nrow(hg)){ # nrow(hg)
  tryCatch({
    #Load coefficients
    setwd(paste0("/gpfs/fs0/scratch/j/junpark/junpark/UKBB_reprocess/imgGenet_coefMat/chr",chr))
    coef=readRDS(paste0("gene",gene,".rds"))
    pvecs_gene <- rep(NA, ncol(coef$S))
    pvecs2_gene <- rep(NA, ncol(coef$S))
    pvecs3_gene <- rep(NA, ncol(coef$S))
    pvecs4_gene <- rep(NA, ncol(coef$S))
    
    for (IDPi in 1:ncol(coef$S)){
      tryCatch({
        A1=as.matrix(cbind(coef$S))
        IDPnames <- colnames(A1)
        A1=as.matrix(A1[, IDPi])
        A2=as.matrix(cbind(coef$S[, -IDPi], coef$D, coef$F))
        
        #Load LD
        setwd(paste0("/gpfs/fs1/home/j/junpark/junpark/imgGenet/chr",chr))
        LD=readRDS(paste0("combine_gene",gene,".rds"))$LD
        index=match(rownames(A1), rownames(LD))
        LD=LD[index,index]
        
        #Load GWAS summary statistics
        setwd(paste0("/gpfs/fs0/project/j/junpark/tianyu47/UKB_IGAPset_ADGWAS/chr",chr))
        summ=readRDS(paste0("UKB_IGAPset_ADGWAS_chr",chr,"_gene", gene,".rds"))
        index=match(rownames(A1), summ$rsID)
        summ=summ[index,]
        
        pos=chrinfo$chromEnd[match(summ$rsID, chrinfo$name)]
        summ=chrstat[match(pos, chrstat_pos),c(1,11)]
        corr=as.numeric(as.matrix(summ[,2]))
        pvecs_gene[IDPi]=mv_vc_iwas(corr,LD,A1,A2)
        pvecs2_gene[IDPi]=mv_vc_iwas(corr,LD,A1,A2=NULL)
        pvecs3_gene[IDPi]=mv_vc_iwas(corr,LD,A1,A2=NULL, egger = F)
        pvecs4_gene[IDPi]=mv_vc_iwas(corr,LD,A1,A2, egger=F)
        names(pvecs_gene)[IDPi]=IDPnames[IDPi]
        names(pvecs2_gene)[IDPi]=IDPnames[IDPi]
        names(pvecs3_gene)[IDPi]=IDPnames[IDPi]
        names(pvecs4_gene)[IDPi]=IDPnames[IDPi]
      }, error = function(e) {print("error")})
    }
    pvecs[[gene]]=pvecs_gene
    pvecs2[[gene]]=pvecs2_gene
    pvecs3[[gene]]=pvecs3_gene
    pvecs4[[gene]]=pvecs4_gene
  }, error = function(e) {
    pvecs[[gene]]=NA
    pvecs2[[gene]]=NA
    pvecs3[[gene]]=NA
    pvecs4[[gene]]=NA
  })
}

pvalues <- list(pvecs=pvecs, pvecs2=pvecs2, pvecs3=pvecs3, pvecs4=pvecs4)

setwd("/gpfs/fs0/scratch/j/junpark/tianyu47/MV_VC_IWAS/MVIWAS/results")
saveRDS(pvalues, paste0("chr",chr,"_UKBresults_S.rds"))