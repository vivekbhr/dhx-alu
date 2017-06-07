


## Here I am just trying some variations in the way i cluster the repeat Lib mapping results, 
# checking to see if the results look better

# load the source functions
source("batchEnrichment.R")
# Load human data
human_res.merged <- multmerge(mypath = "/data/akhtar/bhardwaj/2015_Ibu_dhx9-alu/Scripts/Rscripts_Dhx-Alu/HML_enrichments/Replicates_merged")
# mouse data 
mice_res.merged <- multmerge(mypath = "/data/akhtar/bhardwaj/2015_Ibu_dhx9-alu/Scripts/Rscripts_Dhx-Alu/MOUSE/libMapping_results/rep_merged")

pdf("Different_clustering_types-Hum-mice.pdf")
for(Rowclus in c(TRUE,FALSE)){
    for(scale in c("row","col")){
      name = ifelse(lib )
      plotMLE(calcMLEoutput = human_res.merged,RowCorrClustering = Rowclus,scale = scale)
      legend("top",legend = paste0("Human_","Corr-Clust=",Rowclus,"_scale=",scale))
      plotMLE(calcMLEoutput = human_res.merged,aluOnly = TRUE,RowCorrClustering_Alu = Rowclus,scale = scale)
      legend("top",legend = paste0("Human_","Corr-Clust=",Rowclus,"_scale=",scale))
    }
}
for(Rowclus in c(TRUE,FALSE)){
  for(scale in c("row","col")){
    name = ifelse(lib)
    plotMLE_mice(calcMLEoutput = mice_res.merged,RowCorrClustering = Rowclus,scale = scale)
    legend("top",legend = paste0("Mice_","Corr-Clust=",Rowclus,"_scale=",scale))
    plotMLE_mice(calcMLEoutput = mice_res.merged,B_only = TRUE,RowCorrClustering_Alu = Rowclus,scale = scale)
    legend("top",legend = paste0("Mice_","Corr-Clust=",Rowclus,"_scale=",scale))
  }
}
dev.off()