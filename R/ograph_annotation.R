
##################################
## read annotation file and add a 'genes' attribute to graph node
##################################
readMappings <- function(file, sep = "\t", IDsep = ",") {
  a <- read.delim(file = file, header = FALSE,
                  quote = "", sep = sep, colClasses = "character")
  
  ## a bit of preprocesing to get it in a nicer form
  map <- a[, 2]
  names(map) <- gsub(" ", "", a[, 1]) ## trim the spaces
  
  ## split the IDs
  return(lapply(map, function(x) gsub(" ", "", strsplit(x, split = IDsep)[[1]])))
}



#' @title mapGene2Graph
#' @description read annotation file and add a 'genes' attribute to graph node
#' @aliases readMappings mapGene2Graph
#' @usage 
#' ...
#' @param graph igraph object represent the ontology graph
#' @param rollup whether the gene will be rolled up
#' @param file annotation file
#' @author Xin He
mapGene2Graph<-function(graph,file=system.file("extdata/annotation","human_gene2HDO_o", package ="topOnto"),rollup=TRUE){
  geneID2Node <- readMappings(file)
  rId <- unlist(geneID2Node, use.names = FALSE)
  lId <- rep(names(geneID2Node), sapply(geneID2Node, length))
  GO2geneID<-split(lId, rId)
  
  geneTerms <- new.env(hash = TRUE, parent = emptyenv())
  allNodes<-V(graph)$name
  lapply(allNodes,
         function(x) {
           e <- new.env(hash = TRUE, parent = emptyenv())
           if(length(GO2geneID[[x]])>0)
             multiassign(GO2geneID[[x]], rep(TRUE, length(GO2geneID[[x]])), envir = e)           
           assign(x, e, envir = geneTerms)
         }) 
  
  
  if(rollup==TRUE){
    ## get the levels list
    nodeLevel=buildLevels(graph)
    levelsLookUp <- nodeLevel$level2nodes
    noOfLevels <- nodeLevel$noOfLevels
    
    for(i in noOfLevels:1) {
      cat(paste("rolling up annotation of level",i,"/",noOfLevels,"\n",sep=' '))
      currentNodes <- levelsLookUp[[as.character(i)]]
      
      ## get all the adjacent nodes (teoreticaly nodes from level i - 1)  
      adjList<-get.adjlist(graph,mode='out')
      .p<-adjList[currentNodes]
      
      
      ## push the genes from level i to level i - 1
      lapply(currentNodes,
             function(node) {
               ## get the genes from this node
               genesID <- ls(get(node, envir = geneTerms, mode = 'environment'))
               
               ## debug option, just in case something goes wrong
               #if(length(genesID) == 0)
               #   print(i)
               
               ## for each adiacent node mapp the genesID to them
               if(length(genesID) > 0){
                 lapply(V(graph)[.p[[node]]]$name,
                        function(destNode) {
                          destEnv <- get(destNode, envir = geneTerms, mode = 'environment')
                          multiassign(genesID, rep(FALSE, length(genesID)), envir = destEnv)
                          return(NULL)
                        })
                 return(NULL)
               }})
    }
  }

  
  V(graph)$genes<-as.list(geneTerms)[allNodes]
  
  graph
}

mapGene2Graph2<-function(graph,file){
  #geneID2Node <- readMappings(file = system.file("examples/annotation","human_d2g", package ="topOnto"))
  geneID2Node <- readMappings(file)
  rId <- unlist(geneID2Node, use.names = FALSE)
  lId <- rep(names(geneID2Node), sapply(geneID2Node, length))
  GO2geneID<-split(lId, rId)
  
  allNodes <- V(graph)$name
  V(graph)$genes<-GO2geneID[allNodes]
  graph
}

##################################
# Roll up annotation in igraph
#
# mode=c(1,2,3)
# 1 abs value : will roll up regardless of the leaf deep
# 2 even : will roll up the deepest leavs only
#
# limit : stop rolling when reach a certain level

rollUpAnnotation<-function(graph,mode=1,times=1,limit=3){

  
  ##build levels
  levels<-buildLevels(graph=graph)
  
  for(x in 1:times){
    ##all leaves
    leaves<-findLeafNode(graph)
    ##only roll up leaves that are deep enough
    leaves<-names(levels$nodes2level[leaves][levels$nodes2level[leaves]>limit])
    ##only roll up the deepest leaves
    if(mode==2){
      #find the deepest leaf    
      deep<-sort(unlist(levels$nodes2level[leaves]),decreasing=TRUE)
      leaves<-names(deep[deep==deep[[1]]])
      cat(paste("Rolling up from level ",deep[[1]]," to ",deep[[1]]-1,"...\n",sep=''))
    }
    
    ##if reach top level, stop
    if(length(leaves)==0) stop('can not roll up any more. Top researched!')
    
    cat(paste(length(leaves)," nodes will be rolled up.\n",sep=''))
    
    adj<-get.adjlist(graph,mode='out')
    
    geneTerms <- new.env(hash = TRUE, parent = emptyenv())
    allNodes<-V(graph)$name
    multiassign(allNodes,value=V(graph)$genes,envir=geneTerms)
    
    #     ## get the levels list
    #     levelsLookUp <- levels$level2nodes
    #     noOfLevels <- levels$noOfLevels
    
    lapply(leaves,
           function(node) {
             ## get the genes from this node
             genesID <- ls(get(node, envir = geneTerms, mode = 'environment'))
             
             ## for each adjacent node mapp the genesID to them
             if(length(genesID)>0)
               lapply(adj[[node]],
                      function(destNode) {
                        dest<-V(graph)[destNode]$name
                        destEnv <- get(dest, envir = geneTerms, mode = 'environment')
                        multiassign(genesID, rep(FALSE, length(genesID)), envir = destEnv)
                        V(graph)[dest]$genes<-list(destEnv)
                        return(NULL)
                      })
             return(NULL)
           })
    ##delete leaves
    graph<-delete.vertices(graph,which(V(graph)$name %in% leaves))
  }
  
  
  graph
}

rollUpAnnotation2<-function(graph,mode=1,times=1,limit=3){
  levels<-buildLevels(graph=graph)
  for(x in 1:times){
    ##all leaves
    leaves<-findLeafNode(graph)
    ##only keep leaves that are deep enough
    leaves<-names(levels$nodes2level[leaves][levels$nodes2level[leaves]>limit])
    
    ##keep leave base on deep %
    if(mode==2){
      #find the deepest leaf    
      deep<-sort(unlist(levels$nodes2level[leaves]),decreasing=TRUE)
      leaves<-names(deep[deep==deep[[1]]])
      cat(paste("Rolling up from level ",deep[[1]]," to ",deep[[1]]-1,"...\n",sep=''))
    }
    
    
    ##if reach top level, stop
    if(length(leaves)==0)
      stop('can not roll up any more. Top researched!')
    
    
    cat(paste(length(leaves)," nodes will be rolled up.\n",sep=''))
    
    adj<-get.adjlist(graph,mode='out')
    
    geneTerms <- new.env(hash = TRUE, parent = emptyenv())
    allNodes<-V(graph)$name
    lapply(allNodes,
           function(x) {
             e <- new.env(hash = TRUE, parent = emptyenv())
             if(length(V(graph)[x]$genes[[1]])>0)
               multiassign(V(graph)[x]$genes[[1]], rep(TRUE, length(V(graph)[x]$genes[[1]])), envir = e)           
             assign(x, e, envir = geneTerms)
           })
    ls(get('DOID:1307',geneTerms))
    
    for(node in leaves){
      genesID <- ls(get(node, envir = geneTerms))
      if(length(genesID)>0)
        lapply(adj[[node]],function(destNode) {
          dest<-V(graph)[destNode]$name
          destEnv <- get(dest, envir = geneTerms, mode = 'environment')
          multiassign(genesID, rep(FALSE, length(genesID)), envir = destEnv)
          return(NULL)
        })
    }
    
    ##delete leaves
    graph<-delete.vertices(graph,which(V(graph)$name %in% leaves))
    
    ##copy back the annotation
    V(graph)$genes<-as.list(geneTerms)[V(graph)$name]
  }
  
  graph
}

.get.genesInNodes<-function(graph,nodes){
  x <- nodes %in% V(graph)$name
  if(!all(x)) {
    warning("Nodes not present in the graph:", nodes[!x])
    nodes <- nodes[x]
  }
  retValue<-lapply(V(graph)[nodes]$genes,ls)
  names(retValue)<-nodes
  retValue
}

##################################
# Roll up to a certain level
rollUpToLevel<-function(graph,mode=2,level=5){
  levels<-buildLevels(graph)
  #V(graph)$name %in% unlist(levels$level2nodes[names(levels$level2nodes)>level],use.names=FALSE)
  while(length(intersect(V(graph)$name,unlist(levels$level2nodes[names(levels$level2nodes)>level],use.names=FALSE)))>0)
    graph<-rollUpAnnotation(graph,mode=mode,times=1)
  graph
}


##################################
# output gene annotation from graph
saveAnnotationFromGraph<-function(graph,file){
  sink(file)
  graph=subg
  allnodes=V(graph)$name
  GO2geneID=.get.genesInNodes(graph,allnodes)
  for(i in allnodes){
    gs=GO2geneID[[i]]
    if(length(gs)>0){
      ps_str=paste(gs,collapse=',')
      cat(paste(i,ps_str,sep='\t'),"\n")
    }
  }
  sink()
}

##################################
#map entrez id two gene symbol
entrez2symbol<-function(file=''){  
  if(file=='')
    file=system.file("extdata","entrez2symbol.txt", package ="ograph")
  
  lines<-read.delim(file,header=TRUE,sep=",")
  index<-which(!as.vector(lines$HGNC.symbol) %in% '')
  entrez_ids<-as.vector(lines$EntrezGene.ID)[index]
  symbols<-as.vector(lines$HGNC.symbol)[index]
  e2s<-symbols
  names(e2s)<-entrez_ids
  e2s
}


# l<-lapply(GO2geneID,function(x){
#   unname(e2s[x][!is.na(e2s[x])])
# })
# terms<-Term(ONTTERM)
# for(i in names(l)){
#   if(length(l[[i]])>0){
#   genes<-paste(l[[i]],collapse='\t')
#   cat(paste(i,terms[i],genes,sep='\t'))
#   cat("\n")
#   }
# }



.getTermsDefinition <- function(igraph,whichTerms, numChar = 20, multipLines = FALSE) {
  
  termsNames=ograph@term2def[whichTerms]
  
  if(!multipLines) 
    shortNames <- paste(substr(termsNames, 1, numChar),
                        ifelse(nchar(termsNames) > numChar, '...', ''), sep = '')
  else
    shortNames <- sapply(termsNames,
                         function(x) {
                           a <- strwrap(x, numChar)
                           return(paste(a, sep = "", collapse = "\\\n"))
                         })
  
  names(shortNames) <- names(termsNames)
  return(shortNames[whichTerms])
}
