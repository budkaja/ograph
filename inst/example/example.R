library('topOnto')
topOnto::initONT('HDO')

a<-system.file("extdata/annotation","human_gene2HDO", package ="topOnto")
g<-system.file("extdata/genelist","age", package ="topOnto")

geneID2GO <- readMappings(file = a)
geneNames=names(geneID2GO)
myInterestingGenes=(read.csv(header = FALSE, file = g))$V1
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames


GOdata <- new("topGOdata", ontology = "HDO", allGenes = geneList,annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultElimFis<- runTest(GOdata, algorithm = "elim", statistic = "fisher")
allRes <- GenTable(GOdata,  elim = resultElimFis,elim = resultElimFis,topNodes = 30,useLevels=TRUE)
print(allRes)


ts=names(sort(pValue.classic)[1:10])



library('ograph')
ograph::initOGraph('HDO')
g<-new("ontGraph",ontology='HDO')
mg<-ograph::subGraphByNodes(graph=g@graph,ts)
treeplot(mg,label=1)

pValue.classic <- score(resultFis)
pValue.elim <- score(resultElimFis)[names(pValue.classic)]

gstat <- termStat(GOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)


showSigOfNodes(GOdata, score(resultElimFis), firstSigNodes = 5, useInfo = 'all')


graph=subGraphByLevel(g@graph,l=5)
graph=mapGene2Graph(graph)
calculate.upper.adj.matrix.print(graph)
system.time(calculate.upper.adj.matrix(graph))
system.time(calculate.upper.adj.matrix3(graph))

ds=c("DOID:10652","DOID:3312","DOID:12849","DOID:5419","DOID:14330","DOID:12858")
x=read.table(file = "/home/xin/Desktop/SYNSYS_17_02_2015.csv", sep = "\n", header=FALSE)
colin=as.vector(x$V1)
geneID2GO=geneID2GO[colin]
GO2geneID=revmap(geneID2GO)
gs=GO2geneID[ds] 
def=Term(ONTTERM)
def=def[names(gs)]
sink("/home/xin/Desktop/colin.txt")
for(i in names(gs)){
  id=i
  name=def[i]
  genes=gs[[i]]
  for(j in genes){
    cat(paste(id,name,j,sep="\t"))
    cat("\n")
  }
}
sink()



library('ograph')
ograph::initOGraph('HDO')
g<-new("ontGraph",ontology='HDO')
mg<-subGraphByNodes(g@graph,nodes=c('DOID:10652'))
mg<-mapGene2Graph(mg,file=system.file("extdata/annotation","human_gene2HDO", package ="topOnto"))
GO2geneID=.get.genesInNodes(mg,V(mg)$name)
e2s<-entrez2symbol()
l<-lapply(GO2geneID,function(x){
  unname(e2s[x][!is.na(e2s[x])])
})
