library('topOnto')
topOnto::initONT('HDO')

a<-system.file("extdata/annotation","human_gene2HDO", package ="topOnto")
g<-system.file("extdata/genelist","age", package ="topOnto")

g2t <- readMappings(file = a)


##filter the annotation with your list of terms
t2g=revmap(g2t)
ts=scan(file="/home/xin/Desktop/colin/150terms",what='character')
g2t=revmap(t2g[ts])

##write the annotation to file
file="/home/xin/Desktop/colin/neuro_all"
lapply(names(g2t), function(x){
  id=x
  value=paste(g2t[[x]],collapse=',')
  write(paste(id,value,sep='\t'),file, append=TRUE, ncolumns=1000)
})

##



geneNames=names(g2t)
myInterestingGenes=(read.csv(header = FALSE, file = g))$V1
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames


GOdata <- new("topGOdata", ontology = "HDO", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultElimFis<- runTest(GOdata, algorithm = "elim", statistic = "fisher")
allRes <- GenTable(GOdata,  elim = resultElimFis,elim = resultElimFis,topNodes = 30,useLevels=TRUE)
print(allRes)