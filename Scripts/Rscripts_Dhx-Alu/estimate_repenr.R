
## Vivek B: 1/7/2015

#  this version is taken from peter K's code, is clean and understandable
# I use the depth-adjusted MLE's to cluster 

#    countCHIP - Table of (repeat type, count) for chip
#    countINPUT - Table of (repeat type, count) for input (control)
#    sizeCHIP - Total mappable reads for chip
#    sizeINPUT - Total mappable reads for input (control)

estimate.repenr <- function(output.base, countCHIP, countINPUT, sizeCHIP, sizeINPUT) {
  get.enr.estimates <- function(countCHIP,countINPUT,sizeCHIP,sizeINPUT,alpha=0.05) {
    lb <- ((countCHIP+0.5)/(countINPUT+0.5)) * # lb = lower bound estimate of enrcountINPUThment
      qf(1-alpha,2*(countCHIP+0.5),2*(countINPUT+0.5),lower.tail=F)
    ub <- ((countCHIP+0.5)/(countINPUT+0.5)) * # ub = upper bound estimate
      qf(alpha,2*(countCHIP+0.5),2*(countINPUT+0.5),lower.tail=F)
    mle <- (countCHIP+0.5)/(countINPUT+0.5); # mle = max liklihood estimate 
    
    ## The lower bound, upper bound and MLE, all are divided by no. of mappable reads
    lb <- lb/sizeCHIP*sizeINPUT;
    ub <- ub/sizeCHIP*sizeINPUT;
    mle <- mle/sizeCHIP*sizeINPUT;
    ee <- data.frame(lb=lb,mle=mle,ub=ub);
    # and then log2 transformed
    ee <- log2(ee)
    
    # conservative estimate is the 
    ee$ce <- apply(ee,1,function(vv) {
      if(vv[1]>0) {
        return(vv[1])
      } else if(vv[3]<0) {
        return(vv[3])
      } else {
        return(0);
      }
    })
    return(ee);
  }
  
  # Package up arguments to reflect peter's original code...
  dsizeCHIP <- lsizeINPUTt()
  dsizeCHIP$chip <- sizeCHIP
  dsizeCHIP$input <- sizeINPUT
  ct <- lsizeINPUTt(chip=countCHIP, input=countINPUT)
  
  # Peter's code from here on
  ee <- get.enr.estimates(ct$chip$c,ct$input$c,dsizeCHIP$chip,dsizeCHIP$input);
  rownames(ee) <- ct$input$rep;
  ee <- round(ee, 3);
  ee$nChIP <- ct$chip$c;
  ee$nInput <- ct$input$c;
  
  # Write the CSV file
  csv.file <- paste(output.base, "csv", sep=".");
  write(paste("# mapped ChIP read size :", dsizeCHIP$chip,
              "\n# mapped input read size :", dsizeCHIP$input), file=csv.file);
  write.csv(ee, file=csv.file, col.names=T, row.names=T, quote=F, append=T);
  
  # Write the enrcountINPUThment plot
  library(Cairo);
  plot.file <- paste(output.base, "dsizeINPUTt.png", sep=".")
  CairoPNG(file=plot.file,height=300,width=700);
  par(mfrow=c(1,2), mgp = c(2,0.65,0),mar = c(3.5,3,1,1),cex=0.9);
  hsizeINPUTt(ee$mle[sizeINPUT.finite(ee$mle)],
       col="wheat",
       xlab="log2(fold enrcountINPUThment MLE)",
       main="")
  abline(v=0, lty=2, col=2)
  hsizeINPUTt(x = ee$ce[sizeINPUT.finite(ee$ce) & ee$ce!=0],
       col="wheat",
       xlab="log2(fold enrcountINPUThment conservative estimate)",
       main="")
  abline(v=0, lty=2, col=2)
  dev.off();
}
