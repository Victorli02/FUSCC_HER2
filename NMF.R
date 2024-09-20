






# NMF for HER2 project

library(MyEasyNMF)
library(stringr)
library(parallel)

setwd(NMF_dir)

# Organize data -----------------------------------------------------

EXP_Tumor <- RNA_data[, str_detect(colnames(RNA_data), fixed("_T"))]
colnames(EXP_Tumor) <- substr( colnames(EXP_Tumor), 1, 4)


Tes_Cases <- clinical_data$PatientCode[!is.na(clinical_data$HER2_status) & clinical_data$HER2_status == "HER2_positive"  ]


  Tes_Cases <- intersect( Tes_Cases, colnames(EXP_Tumor) ) 
  Tes_EXP  <-  EXP_Tumor[, Tes_Cases]



# Prepare gct files  ------------------------------------------------------

dir.create("tmp")
setwd("tmp")



## Expression --------------------------------------------------------------------------
exp.fpkm.TT      <- log2(Tes_EXP+1)
exp.fpkm.TT.SD = exp.fpkm.TT[order(apply(exp.fpkm.TT, 1, sd),decreasing = T)[1:mRNA_count_filter],] 

exp.fpkm.TT.SD <- t(apply(t(exp.fpkm.TT.SD), 2, scale)) 
colnames(exp.fpkm.TT.SD) <- colnames(Tes_EXP)

tmp = cbind(NAME = rownames(exp.fpkm.TT.SD), Description = rownames(exp.fpkm.TT.SD), exp.fpkm.TT.SD)
tmp = rbind(c("#1.2",rep( "" , ncol(tmp)-1) ),
            c(nrow(exp.fpkm.TT.SD),ncol(exp.fpkm.TT.SD),rep( "" , ncol(tmp)-2) ),
            colnames(tmp), tmp)

write.table(tmp,sep = "\t",
            file = paste0("./exp.fpkm.TT.SD.gct"),quote = F,col.names = F,row.names = F)

save(exp.fpkm.TT.SD,file = paste("RNA",RNA_SD_percent,'_percent_NMF_',mRNA_count_filter,'RNA.Rdata',sep ="" ))



# Generate gct files ------------------------------------------------------

NMF_taskname <- paste("RNA_",mRNA_count_filter,'_',RNA_SD_percent,'_K8_percent_SD_RNA_NMF',sep ="" )


  tmp = data.frame(V1 = c("RNA"),
                   V2 = c("exp.fpkm.TT.SD.gct")  )


write.table(tmp ,file = './nmf.conf',sep = "\t", row.names = F, col.names =F, quote = F)
tar(tarfile = paste(NMF_dir, NMF_taskname, ".tar", sep = ""), files = "./")


# NMF ---------------------------------------------------------------------

setwd(NMF_dir)
unlink(paste(NMF_dir, "tmp", sep = ""),recursive = T)

time1 = Sys.time()
easyNMF(tar_file = paste(NMF_dir, NMF_taskname, ".tar", sep = ""),
        output_dir = paste(NMF_dir, NMF_taskname, "/",sep = ""),
        kmax = kmax,nrun = nrun )
time2 = Sys.time()
print(time2 - time1)
