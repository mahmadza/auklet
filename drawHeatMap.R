
library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )

#impurt output from flattenHeatMap.sh
d<-read.table(input_file)     #define your input file here

#create files
png(filename=output_file,width=8.5,height=11,pointsize=10,units="in",res=150)     #define your output file here, as PNG

title=paste0("STARR-seq centered around\nXXX coordinates")      #put your title here

heatmap.2(as.matrix(d),
main=title,
xlab="STARR-seq coverage",
ylab=paste0("XXX peaks"),
Rowv=F,Colv=F,
dendrogram="none",trace="none",density.info="none",
col=colorpanel(20,rgb(49,105,154,25,maxColorValue=255),rgb(255,255,255,25,maxColorValue=255),rgb(164,25,32,25,maxColorValue=255)))

dev.off()



