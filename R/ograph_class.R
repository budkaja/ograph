##
## The class represent an ontology in igraph
##
setClass("ontGraph",
         slots = c(
           ## some description of the data/experiment
           #description = "character",
           ## which of the ontology to be used: HDO...
           ontology = "character",
           ##igraph pbject
           graph = "ANY",
           ##node level list
           levels = "ANY",
           ##color2level
           color2level = "ANY",
           ##term id2def
           term2def = "ANY"
         )
)

setMethod("initialize", "ontGraph",function(.Object,ontgraph=NULL,ontology='HDO',annotation_file=''){
  
  .Object@ontology<-ontology
  #.Object@description<-description
  
  initOGraph(ontology = ontology)
  .Object@term2def=Term(ONTTERM)
  
  ##build the ontology
  if(!is.null(ontgraph)){
    .Object@graph<-ontgraph
  }else{
    ##init the topOnto.xx.db for ontology data that can be used to 
    ##build ograpg object
    initOGraph(ontology=ontology)
    .Object@graph<-buildOntGraph(toTable(ONTTERM),toTable(ONTPARENTS))
  }
  
  ##build nodes levels
  .Object@levels <- buildLevels(.Object@graph)
  
  ##add graph attributes
  .Object@graph<-node.addDefAttribute(.Object@graph)
  .Object@graph<-node.addLevelAttribute(.Object@graph,.Object@levels)
  .Object@graph<-node.addLeafAttribute(.Object@graph)
  .Object@graph<-node.addColorAttributeByLevel(.Object@graph)
  
  
  
  ##.Object@graph
  .Object@color2level=rainbow(.Object@levels$noOfLevels)
  names(.Object@color2level)<-1:.Object@levels$noOfLevels
  
  ##add annotation to node?
  if(file.exists(annotation_file))
    .Object@graph<-mapGene2Graph(.Object@graph,annotation_file)
  
  .Object
})




setMethod("print", "ontGraph", function(x, ...){
  cat("\n------------------------- ontGraph object -------------------------\n")

  cat("\n Ontology:\n")
  cat("   - ", x@ontology, "\n")
  
  cat("\n Graph:\n")
  summary(x@graph)
  
  cat("\n levels:\n")
  summary(x@levels)
})
