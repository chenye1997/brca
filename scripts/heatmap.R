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

#pdf(file = "figure1temp.pdf",   
#    width = 6, 
#    height = 12)
#par(mfrow = c(4,2))
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
## there are two heatmaps mapping figdat1 and figdat2
sepBRCA <- 10
output <- output1# output <- output1[-which(output1$k==39),]  ##39 has too large value
pdf(file = "new_figure1heat.pdf",   
    width = 22, 
    height = 22)
par(mfrow = c(1,2))
library(ggplot2)
library(reshape2)
library(gplots) 
library(scales)
makeheatmap <- function(){for (jj in unique(output$j)) { #
  fig1 <- fig2 <- NULL#as.data.frame(matrix(NA, 1, sepBRCA))
  for (i in unique(output$i)){ 
    id <- which(ctypes == alctypes[i])
    j <- jj
    temp <- t(as.matrix(datobj[c(idr[j],which(rownames(datobj) == "BRCA1"),idr[unique(output$k)]),id]) )
    medianid <- which(temp[,1]<median(temp[temp[,1]!=0,1]))
    if ((length(medianid) >= sepBRCA+1)&((nrow(temp) - length(medianid) )>= sepBRCA+1)) {
      tmp1 <- temp[medianid,-1] ; tmp2 <- temp[-medianid,-1]
      figdat1 <- grpcrostab(dat = tmp1, grps = sepBRCA, decreasing = FALSE)
      figdat2 <- grpcrostab(dat = tmp2, grps = sepBRCA, decreasing = TRUE)
      mm <- round(nrow(figdat1)/2)
      mn <- nrow(figdat1)
      # I add i for the column name cause ggplot will ignore the same row name
      rownames(figdat1)[1:(mm-1)] <- rownames(figdat2)[1:(mm-1)] <- paste(rownames(figdat1)[1:(mm-1)],"/",i)
      rownames(figdat1)[(mm+1):mn] <- rownames(figdat2)[(mm+1):mn] <- paste(rownames(figdat1)[(mm+1):mn],"\\",i)
      rownames(figdat1)[mm] <- rownames(figdat2)[mm] <- paste(rownames(figdat1)[mm],
                                                              alctypes[i], sep = "  Cell type: ")
      colnames(figdat1) <- colnames(figdat2) <- paste(100*1:sepBRCA/sepBRCA,"%th quantile",sep="")
      fig1 <- rbind(fig1,figdat1)
      fig2 <- rbind(fig2,figdat2)}
  }
  #  fig1_data<-append(fig1_data,fig1)
  #  fig2_data <- append(fig2_data,fig2)
  library(gplots)  
  fig1 <- log(1+cbind(fig1))
  fig1 <- melt(fig1)
  m1 <- median(fig1$value)
  f <- ggplot(fig1, aes( Var2,
                         Var1,fill=value
  ))+geom_tile() + scale_fill_gradient2(
    low = "green", 
    mid = "black", 
    high = "red", 
    midpoint = m1
  )+labs(title = paste("Low level of", rownames(datobj)[idr[j]], collapse = " "),
         x = "BRCA",y = "")
  fig2 <- log(1+cbind(fig2))
  fig2 <- melt(fig2)
  m2 <- median(fig2$value)
  f1 <- ggplot(fig2, aes(x = Var2,
                         y = Var1,
                         fill = value))+geom_tile()+scale_fill_gradient2(
                           low = "green", 
                           mid = "black", 
                           high = "red", 
                           midpoint = m2
                         )+labs(title = paste("High level of", rownames(datobj)[idr[j]], collapse = " "),
                                x = "BRCA",y="")
  combined <- f+f1
  show(combined)
  # heatmap(log(1+cbind(fig1,fig2)), col = redgreen(250),
  #         Rowv = NA, Colv = NA, scale = "row", xlab = "BRCA",
  #         main = paste(c("Low level of", "|       High level of"), rownames(datobj)[idr[j]], collapse = "       "),
  #         margins = c(10,16))
  # 
  # heatmap(log(1+cbind(fig1)), col = redgreen(250),
  #         Rowv = NA, Colv = NA, scale = "row", xlab = "BRCA",
  #         main = paste("low level of", rownames(datobj)[idr[j]], collapse = " "),
  #         margins = c(10,16))
  # 
  # heatmap(log(1+cbind(fig2)), col = redgreen(250),
  #         Rowv = NA, Colv = NA, scale = "row", xlab = "BRCA",
  #         main = paste("High level of", rownames(datobj)[idr[j]], collapse = " "),
  #         margins = c(10,16))
}
}
makeheatmap()
dev.off()


#par(mfrow = c(1,2))
#barplot(fig1[10,])
#barplot(fig2[10,])
