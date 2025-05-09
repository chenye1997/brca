setwd("~/brca/")
T59_annotations <- read.delim("~/brca/T59_annotations.txt", stringsAsFactors=TRUE)
features <- read.table(file = 'features.tsv', sep = '\t', header = TRUE)
library(Seurat)

datobj <- ReadMtx( 
  mtx = "matrix.mtx",
  features = "features.tsv",
  cells = "barcodes.tsv.gz")

features$MIR1302.2HG[1:5]

grep("BRCA",features$MIR1302.2HG)

grep("BRCA",c("12BRCA12","12BRCA"))
features$MIR1302.2HG[c(22030, 27662)]

dim(datobj)

length(T59_annotations$barcode)
length(T59_annotations$celltyperBPEType)
length(features$MIR1302.2HG)






figure1 <- datobj[c(which(rownames(datobj) == "BRCA1"),
                    which(rownames(datobj) == "BRCA1")),]
idr <- order(rowMeans(datobj), decreasing = TRUE)
plot(datobj[idr[1],],figure1[2,])
ctypes <- NA
for (i in 1:ncol(datobj)){
  x <- which(T59_annotations$barcode == colnames(datobj)[i])
  if (length(x) == 0) {temp <- NA} else {temp <- T59_annotations$celltyperBPEType[x]}
  ctypes[i] <- paste(temp)
}

alctypes <- paste(unique(T59_annotations$celltyperBPEType))
output1 <- read.csv("heatmap_data.csv")

grpcrostab <- function(dat, grps = 10, decreasing = FALSE){
  for (i in ncol(dat):2){
    dat <- dat[order(dat[,i], decreasing = decreasing),]
  }
  dat <- dat[order(dat[,1]),]
  result <- matrix(NA,ncol(dat)-1,grps)
  rownames(result) <- colnames(dat[,-1])
  grpsid <- round(seq(1,nrow(dat),length = grps+1))
  for (i in 1:grps){
    result[,i] <- colMeans(dat[grpsid[i]:grpsid[i+1],-1])
  }
  result
}

genelist <- unique(c(unique(output1$j),unique(output1$k)))
library(dagbagM)
celllist <- (22)
cellindex <-  NA
for (cell_name_index in 1:length(celllist)){
  for(num in 1:length(ctypes)){
      if (ctypes[num]==alctypes[celllist[cell_name_index]]){
        cellindex <- append(cellindex,num)
      }
  }
}
brca_index <- which(rownames(datobj) == "BRCA1")
cellindex<-cellindex[!is.na(cellindex)]
dag_matrix <- t(as.matrix(datobj[c(genelist,brca_index),cellindex]))
#dag_matrix1 <- datobj[genelist,cellindex]
dag_matrix <-grpcrostab(dag_matrix,grps=5)
dag_matrix <-t( dag_matrix[rowSums(dag_matrix[])>0,])
library(doParallel)
dag_model<- dagbagM::hc(Y=dag_matrix, whiteList=NULL, blackList=NULL, tol = 1e-6, standardize=TRUE, maxStep = 100, restart=10, seed = 1,  verbose = FALSE)
library(igraph)
dag <- graph_from_adjacency_matrix(dag_model$adjacency,mode="directed")
plot(dag)
