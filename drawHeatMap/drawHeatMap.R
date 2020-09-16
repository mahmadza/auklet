
library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )

#import output from flattenHeatMap.sh
d<-read.table(input_file)          #define your input file here
experiment="STARR-seq"             #define the experiment here, ex. ChIP-seq or STARR-seq or DHS-seq etc
ref_coordinates="XXX ChIP-seq"     #define where you get the reference coordinates from, ex. ChIP-seq or STARR-seq or DHS-seq peaks etc

#create files
png(filename=output_file,width=8.5,height=11,pointsize=10,units="in",res=150)     #define your output file here, as PNG
#if do PDF, the thousands of points drawn can cause rendering in ex. Adobe Illustrator REALLY slow

title=paste0("STARR-seq centered around\n",ref_coordinates," coordinates")      #put your title here

#draw heatmap
heatmap.2(as.matrix(d),
main=title,
xlab=paste0(experiment," coverage"),
ylab=paste0(ref_coordinates," peaks"),
Rowv=F,Colv=F,
dendrogram="none",trace="none",density.info="none",
col=colorpanel(20,rgb(49,105,154,25,maxColorValue=255),rgb(255,255,255,25,maxColorValue=255),rgb(164,25,32,25,maxColorValue=255)))

dev.off()



