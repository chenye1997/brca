setwd("~/Downloads/R21Data")

T59_annotations <- read.delim("~/Downloads/R21Data/T59_annotations.txt", stringsAsFactors=TRUE)
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

pdf(file = "figure1temp.pdf",   
    width = 6, 
    height = 12)
par(mfrow = c(4,2))
for (i in 1:length(alctypes)){ #7
  id <- which(ctypes == alctypes[i])
  for (j in 11:200){ #10
    for (k in 11:200){ #1:8
      if (k!=j){
        temp <- cbind(datobj[idr[k],id],datobj[idr[j],id],figure1[1,id]) 
        tmp <- median(temp[temp[,1]!=0,1])
        dat1 <- as.data.frame((temp[which(temp[,1] < tmp),2:3]))
        dat2 <- as.data.frame((temp[which(temp[,1] >= tmp),2:3]))
        colnames(dat1)[1] <- colnames(dat2)[1] <- "y"
        tryCatch({  
          tmp1 <- summary(lm(y~.,dat1))$coefficients
          tmp2 <- summary(lm(y~.,dat2))$coefficients
          if (nrow(tmp1) + nrow(tmp2) == 4){
            if ((sign(tmp1[2,1])!=sign(tmp2[2,1]))&(tmp1[2,4]<0.05)&(tmp2[2,4]<0.05)) {
              print(paste(i,j,k))      
              plot(dat1, 
                   sub = paste("celltype:",alctypes[i],i,",",j,",",k),
                   main = paste(rownames(datobj)[idr[k]],"low"),
                   ylab = "BRCA1", xlab = rownames(datobj)[idr[j]])
              abline(lm(y~.,dat1), col='red')
              # lines(lowess(dat1[,1], dat1[,2]), col='red')
              plot(dat2, 
                   sub = paste("celltype:",alctypes[i],i,",",j,",",k),
                   main = paste(rownames(datobj)[idr[k]],"high"),
                   ylab = "BRCA1", xlab = rownames(datobj)[idr[j]])   
              abline(lm(y~.,dat2), col='blue')
              # lines(lowess(dat2[,1], dat2[,2]), col='blue')
            }}}, error=function(e){})  
      }
    }
  }
}
dev.off()


## output of above

[1] "1 54 18"
[1] "1 134 173"
[1] "1 159 30"
[1] "1 174 187"
[1] "1 194 36"
[1] "10 89 49"
[1] "10 102 27"
[1] "10 114 22"
[1] "10 114 27"
[1] "10 114 31"
[1] "10 114 41"
[1] "10 114 82"
[1] "10 114 108"
[1] "10 181 28"
[1] "15 134 76"
[1] "15 134 108"
[1] "15 134 130"
[1] "16 75 112"
[1] "17 49 12"
[1] "17 49 56"
[1] "17 156 85"
[1] "17 197 174"
[1] "22 23 167"
[1] "22 41 12"
[1] "22 41 39"
[1] "22 41 49"
[1] "22 41 78"
[1] "22 41 95"
[1] "22 41 96"
[1] "22 41 112"
[1] "22 41 154"
[1] "22 41 160"
