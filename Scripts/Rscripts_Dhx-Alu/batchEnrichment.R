
## Vivek bhardwaj, 28th Aug 2015
## Set of functions to calculate and plot MLE estimates for human mouse and fly
#load stuff
library("ggplot2")
library("reshape2")
library("gridExtra")
library('plyr')
#library("gplots")
library('RColorBrewer')
source("/data/akhtar/bhardwaj/2015_Ibu_dhx9-alu/Scripts/Dey_et.al_pipeline/R/estimate_repenr.R")

## FUNCTIONS

# Calculating MLEs in batch
calcMLE <- function(repeatCountFolder,mappedCountFile,matchFile,outfolder){
  # prepare totalcoutn file
  bowtie.out <- read.table(mappedCountFile,header = T)
  bowtie.out$mappedreads <- bowtie.out[,4] + bowtie.out[,5]
  bowtie.out = bowtie.out[c(1,6)]
  colnames(bowtie.out) = c("sample","mappedreads")
  # calculate MLE one by one
  matchings <- read.table(matchFile)
  for(i in 1:nrow(matchings)){
    tryCatch({
      chip = as.character(matchings[i,1])
      cont = as.character(matchings[i,2])
      chipcounts = list.files(repeatCountFolder,pattern = chip,full.names = TRUE)
      contcounts = list.files(repeatCountFolder,pattern = cont,full.names = TRUE)
    
      head = as.vector(c("rep","c"))
      chipcounts <- read.table(chipcounts,skip=1,as.is = TRUE)
      colnames(chipcounts) = head
      contcounts <- read.table(contcounts,skip=1,as.is = TRUE)
      colnames(contcounts) = head
      
      chipmapped <- bowtie.out[grep(chip,bowtie.out$sample),2]
      contmapped <- bowtie.out[grep(cont,bowtie.out$sample),2]
      outbase <- paste0(outfolder,"/",chip)
      enrichment <- estimate.repenr(output.base = outbase,sc = chipcounts,ic = contcounts,ss = chipmapped,is = contmapped)
    },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  }
}

# merging files
multmerge <- function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE,pattern=".csv")
  datalist = lapply(filenames, function(x){
    file = read.csv(file=x,header=T)
    nam = gsub(".csv","",basename(x))
    colnames(file) <- c("repfam",paste0(nam,".",colnames(file[2:7])))
    return(file)
  })
  out = Reduce(function(x,y) {merge(x,y,by=1)}, datalist)
  return(out)
}

# plotting mappings
plotmappings <- function(calcMLEoutput,mappedCountFile,select=FALSE,selectArg=NULL){
  # prepare totalcoutn file
  counts <- read.table(mappedCountFile,header = T)
  counts$perc <- ((counts[,5] + counts[,4])/counts[,2])*100
  counts = counts[c(1,2,6)]
  calcMLEoutput = calcMLEoutput[,grep("nChIP",colnames(calcMLEoutput))]
  calcMLEoutput = as.data.frame(colSums(calcMLEoutput))
  calcMLEoutput$sample = gsub(".nChIP","",rownames(calcMLEoutput))
  counts <- merge(counts,calcMLEoutput,by.x = 1,by.y = 2)
  colnames(counts) <- c("sample","totreads","genome","repeat.lib")
  counts$repeat.lib <- (counts$repeat.lib/counts$totreads)*100
  counts <- counts[c(1,3:4)]
  if(select == "TRUE"){
    counts <- counts[grep(selectArg,counts[,1]),]
    counts[,1] <- gsub(selectArg,"",counts[,1])
  } 
  
  counts = melt(counts)
  counts$repmap = counts$value
  counts[which(counts$variable != "repeat.lib"),"repmap"] <- 0
  counts = counts[order(counts$repmap,decreasing = TRUE),]
  counts$repmap <- factor(counts$sample, as.character(counts$sample))
  #levels(counts$sample) = counts$sample
  out <- ggplot(counts,aes(repmap,value,fill=variable)) + geom_bar(stat="identity",position = "dodge") + 
    labs(x="Sample",y="% mapping",fill="Reads mapping to",main="Mapping stats: whole genome vs repeat library") + 
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x=element_text(angle = 80,vjust=0.5,size=12))
  return(out)
}

# plotting MLE
plotMLE <- function(calcMLEoutput,subset=FALSE,whatTosubset=NULL,aluOnly=FALSE,
                    RowCorrClustering = FALSE,RowCorrClustering_Alu = FALSE,
                    dendrogram = "row", scale = "row",...){
    # Clustering
    res.clustering = calcMLEoutput
    # remove repeats with 0 and negative counts
    rownames(res.clustering) <- res.clustering$repfam
    res.clustering = res.clustering[,grep("mle",colnames(res.clustering))]
    # take only some proteins for clustering
    if(subset==TRUE){
      res.clustering <- res.clustering[,grep(paste(whatTosubset,collapse="|"),colnames(res.clustering),value = TRUE)]
    }
    res.clustering[res.clustering < 0] <- 0
    rows <- rowSums(res.clustering)
    res.clustering = res.clustering[rows > 0,]
    colnames(res.clustering) = gsub("([A|B]_merged_.*)_HML.mle","\\1",colnames(res.clustering))
    
    samples = data.frame(tot = colnames(res.clustering),
                         samp=gsub('([A|B])_merged_(.*)','\\2',colnames(res.clustering)),
                         rep=gsub('([A|B])_merged_(.*)','\\1',colnames(res.clustering)))
    samples = samples[order(samples$samp),]
    res.clustering = res.clustering[,order(samples$tot)]
    res.clustering = res.clustering[,order(samples$tot)] # I don't know why I have to do it 2 times, but that's how it works
    
    cols <- colorRampPalette(brewer.pal(n = 10,name = "RdBu"))(256)
    
    if(aluOnly==TRUE){
      res.clustering <- res.clustering[grep("Alu",rownames(res.clustering)),]
      color = ifelse(grepl("AluJ",rownames(res.clustering))==T,"yellow",
                     ifelse(grepl("AluS",rownames(res.clustering))==T, "forestgreen" ,
                            ifelse(grepl("AluY",rownames(res.clustering))==T, "steelblue4","grey40")))
      
      notations <- data.frame(legend= c("AluJ","AluS","AluY","Other"),stringsAsFactors = FALSE,
                              col = c("yellow", "forestgreen" ,"steelblue4","grey40"))
      # Change clustering option if asked
      hr <- hclust(as.dist(1-cor(t(as.matrix(res.clustering)), method="pearson")),method="complete")
      if(RowCorrClustering_Alu){
        rowvalue = as.dendrogram(hr)}
      else{rowvalue = TRUE}
      #plot
      heatmap.2(as.matrix(res.clustering),Rowv=rowvalue,Colv = NULL,dendrogram= "row",srtCol = 40,RowSideColors = color,col = cols,
                density.info = "none",trace="none",cexCol=0.7,keysize=1.0,cexRow=0.7,scale = scale,...)
      legend("topright",fill = notations$col, legend = notations$legend,bty = "n",cex = 0.4,title = "Alu Family")  
                                          
    }else{
      hrdist <- as.dist(1-cor(t(as.matrix(res.clustering)), method="pearson"))
      annot <- data.frame(row.names = rownames(res.clustering),
                          Repeat_type = ifelse(grepl("Alu",rownames(res.clustering))==T,"ALU",
                                         ifelse(grepl("LTR",rownames(res.clustering))==T, "LTR" ,
                                                ifelse(grepl("MER",rownames(res.clustering))==T, "MER",
                                                       ifelse(grepl("Charlie",rownames(res.clustering))==T, "HAT/Charlie",
                                                              ifelse(grepl("LINE",rownames(res.clustering))==T, "LINE",
                                                                     ifelse(grepl('[A-Z]*)n',rownames(res.clustering))==T, 
                                                                            "Simple Repeat","Other")
                                                              ))))))
      #change annotation color
      Var1        <- c("yellow", "forestgreen" ,"steelblue4","darkmagenta","orange3","purple","grey40")
      names(Var1) <- c("ALU","LTR","MER","HAT/Charlie","LINE","Simple Repeat","Other")
      anno_colors <- list(Repeat_type = Var1)
      
      if(RowCorrClustering){
        rowvalue = as.dendrogram(hr)}
      else{rowvalue = TRUE}
      #plot

      pheatmap(res.clustering, clustering_distance_rows = hrdist, fontsize_row = 2, 
               annotation_row = annot, cluster_cols = FALSE, annotation_colors = anno_colors)
    }
    
}

## Plotting MLE for MOUSE samples
plotMLE_mice <- function(calcMLEoutput,subset=FALSE,whatTosubset=NULL,B_only=FALSE,RowCorrClustering= FALSE,RowCorrClustering_B1= FALSE,scale= "row",...){
  # Clustering
  res.clustering = calcMLEoutput
  # remove repeats with 0 and negative counts
  rownames(res.clustering) <- res.clustering$repfam
  res.clustering = res.clustering[,grep("mle",colnames(res.clustering))]
  # take only some proteins for clustering
  if(subset==TRUE){
    res.clustering <- res.clustering[,grep(paste(whatTosubset,collapse="|"),colnames(res.clustering),value = TRUE)]
  }
  res.clustering[res.clustering < 0] <- 0
  rows <- rowSums(res.clustering)
  res.clustering = res.clustering[rows > 0,]
  colnames(res.clustering) = gsub("([A|B]_.*).mle","\\1",colnames(res.clustering))
  
  samples = data.frame(tot = colnames(res.clustering),
                       samp=gsub('([A|B])_(.*)','\\2',colnames(res.clustering)),
                       rep=gsub('([A|B])_(.*)','\\1',colnames(res.clustering)))
  samples = samples[order(samples$samp),]
  res.clustering = res.clustering[,order(samples$tot)]
  res.clustering = res.clustering[,order(samples$tot)] # I don't know why I have to do it 2 times, but that's how it works
  
  cols <- colorRampPalette(brewer.pal(n = 10,name = "RdBu"))(256)
  
  if(B_only==TRUE){
    rowN <- substr(rownames(res.clustering),1,2)
    res.clustering <- res.clustering[grep("B[1|2|3|4]",rowN),]
    color = ifelse(grepl("B1",rownames(res.clustering))==T,"yellow",
                   ifelse(grepl("B2",rownames(res.clustering))==T, "forestgreen" ,
                          ifelse(grepl("B3",rownames(res.clustering))==T, "steelblue",
                                 ifelse(grepl("B4",rownames(res.clustering))==T, "purple","grey40"))))
    
    notations <- data.frame(legend= c("B1","B2","B3","B4","Other"),stringsAsFactors = FALSE,
                            col = c("yellow", "forestgreen" ,"steelblue","purple","grey40"))
    # Change clustering option if asked
    hr <- hclust(as.dist(1-cor(t(as.matrix(res.clustering)), method="pearson")),method="complete")
    if(RowCorrClustering_B1){
      rowvalue = as.dendrogram(hr)}
    else{rowvalue = TRUE}
    #plot
    heatmap.2(as.matrix(res.clustering),Rowv=rowvalue,dendrogram= "row",Colv=NULL,srtCol = 40,RowSideColors = color,col = cols,
              density.info = "none",trace="none",cexCol=0.7,keysize=1.0,cexRow=0.7,scale = scale,...)
    legend("topright",fill = notations$col, legend = notations$legend,bty = "n",cex = 0.4,title = "SINE-B Family")  
    
  }else{
    rowN <- substr(rownames(res.clustering),1,2)
    color = ifelse(grepl("B1",rowN)==T,"yellow",
                   ifelse(grepl("LT",rowN) == T, "forestgreen" ,
                          ifelse(grepl("ME",rowN)==T, "steelblue4",
                                 ifelse(grepl("ML",rowN)==T, "darkmagenta",
                                        ifelse(grepl("HA|X",rowN)==T, "black",
                                               ifelse(grepl("RL",rowN)==T, "lightpink3" ,
                                                      ifelse(grepl("\\([A|T|G|C]*",rowN)==T, "purple","grey80")
                                                      
                                               ))))))
    
    notations <- data.frame(legend= c("B1-SINEs","LTRs","MERs","MLTs","LINEs","RLTRs","Simple Repeats","Other"),stringsAsFactors = FALSE,
                            col = c("yellow", "forestgreen" ,"steelblue4","darkmagenta","black","lightpink3" ,"purple","grey80"))
    # Change clustering option if asked
    hr <- hclust(as.dist(1-cor(t(as.matrix(res.clustering)), method="pearson")),method="complete")
    if(RowCorrClustering){
      rowvalue = as.dendrogram(hr)}
    else{rowvalue = TRUE}
    #plot
    heatmap.2(as.matrix(res.clustering),Rowv=rowvalue,dendrogram= "row",srtCol = 40,RowSideColors = color,col = cols,Colv=NULL,
              density.info = "none",trace="none",cexCol=0.5,keysize=1.0,cexRow=0.001,labRow=NULL,scale = scale,...)
    legend("topright",fill = notations$col, legend = notations$legend,bty = "n",cex = 0.3,title = "Repeat Family")  
  }
  
}

## Plotting MLE for drosophila samples
plotMLE_droso <- function(calcMLEoutput,subset=FALSE,whatTosubset=NULL,aluOnly=FALSE,RowCorrClustering= FALSE,RowCorrClustering_Alu= FALSE,dendrogram= "row", scale="row",...){
  # Clustering
  res.clustering = calcMLEoutput
  # remove repeats with 0 and negative counts
  rownames(res.clustering) <- res.clustering$repfam
  res.clustering = res.clustering[,grep("mle",colnames(res.clustering))]
  # take only some proteins for clustering
  if(subset==TRUE){
    res.clustering <- res.clustering[,grep(paste(whatTosubset,collapse="|"),colnames(res.clustering),value = TRUE)]
  }
  res.clustering[res.clustering < 0] <- 0
  rows <- rowSums(res.clustering)
  res.clustering = res.clustering[rows > 0,]
  colnames(res.clustering) = gsub("(.*).mle","\\1",colnames(res.clustering))
  
  samples = data.frame(tot = colnames(res.clustering),
                       samp=gsub('([A|B])_(.*)','\\2',colnames(res.clustering)),
                       rep=gsub('([A|B])_(.*)','\\1',colnames(res.clustering)))
  samples = samples[order(samples$samp),]
  res.clustering = res.clustering[,order(samples$tot)]
  res.clustering = res.clustering[,order(samples$tot)] # I don't know why I have to do it 2 times, but that's how it works
  
  cols <- colorRampPalette(brewer.pal(n = 10,name = "RdBu"))(256)
  
  if(aluOnly==TRUE){
    res.clustering <- res.clustering[grep("Alu",rownames(res.clustering)),]
    color = ifelse(grepl("AluJ",rownames(res.clustering))==T,"yellow",
                   ifelse(grepl("AluS",rownames(res.clustering))==T, "forestgreen" ,
                          ifelse(grepl("AluY",rownames(res.clustering))==T, "steelblue4","grey40")))
    
    notations <- data.frame(legend= c("AluJ","AluS","AluY","Other"),stringsAsFactors = FALSE,
                            col = c("yellow", "forestgreen" ,"steelblue4","grey40"))
    # Change clustering option if asked
    hr <- hclust(as.dist(1-cor(t(as.matrix(res.clustering)), method="pearson")),method="complete")
    if(RowCorrClustering_Alu){
      rowvalue = as.dendrogram(hr)}
    else{rowvalue = TRUE}
    #plot
    heatmap.2(as.matrix(res.clustering),Rowv=rowvalue,dendrogram= "row",srtCol = 40,RowSideColors = color,col = cols,
              density.info = "none",trace="none",cexCol=0.7,keysize=1.0,cexRow=0.7,scale = scale,...)
    legend("topright",fill = notations$col, legend = notations$legend,bty = "n",cex = 0.4,title = "Alu Family")  
    
  }else{
    color = ifelse(grepl("_DM",rownames(res.clustering))==T,"yellow",
                   ifelse(grepl("LTR",rownames(res.clustering))==T, "forestgreen" ,
                          ifelse(grepl("LINE",rownames(res.clustering))==T, "orange3",
                                    ifelse(grepl('[A-Z]*)n',rownames(res.clustering))==T, "purple","grey40")
                                        )))
    
    notations <- data.frame(legend= c("DMs","LTRs","LINEs","Simple Repeats","Other"),stringsAsFactors = FALSE,
                            col = c("yellow", "forestgreen", "orange3","purple","grey40"))
    # Change clustering option if asked
    hr <- hclust(as.dist(1-cor(t(as.matrix(res.clustering)), method="pearson")),method="complete")
    if(RowCorrClustering_Alu){
      rowvalue = as.dendrogram(hr)}
    else{rowvalue = TRUE}
    #plot
    heatmap.2(as.matrix(res.clustering),Rowv=rowvalue,dendrogram= "row",srtCol = 40,RowSideColors = color,col = cols,
              density.info = "none",trace="none",cexCol=0.5,keysize=1.0,cexRow=0.001,labRow=NULL,scale = scale,...)
    legend("topright",fill = notations$col, legend = notations$legend,bty = "n",cex = 0.3,title = "Repeat Family")  
  }
  
}


