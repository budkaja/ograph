##################################
#' @name buildOntGraph
#' @title buildOntGraph
#' @aliases buildOntGraph
#' @description build the ontology with igraph. 
#' Need to use topOnto.xxx.db
#' @usage 
#' buildOntGraph(toTable(ONTTERM),toTable(ONTPARENTS))
#' @details build the igraph object of the ontology
#' @param nodes a data frame contain the nodes and their definition. 
#' This is usually generate by using the ONTTERM object from the 
#' topOnto.xxx.db package.
#' @param edges a data frame contain the relationship between nodes. 
#' This is usually generate by using the ONTPARENTS object from the 
#' topOnto.xxx.db package.
#' @return an igrph object represent the ontology
#' @author Xin he
#' @keywords buildOntGraph
buildOntGraph<-function(terms,edges){
  
  ##first column needs to be parent terms to use layout.reingold.tilfor
  edges<-edges[(c(1,2,3))]
  names(edges)<-c('to','from','type')
  
  ##For the terms, the table is a one-to-many mapping, because it contains the (potentially many) synonyms and secondary IDS of the category. 
  ##But we throw away these columns anyway, so we throw away the rows that appear many time in the table as well.
  terms <- terms[ !duplicated(terms[,1]), ][2:3]
  
  ##remove nodes with no edges
  goodterms <- terms[terms$id %in%  unique(c(edges$to,edges$from)),]
  ##remove edge with undef nodes
  edges <- edges[edges$from %in% terms$id & edges$to %in% terms$id,]
  
  if(length(terms$id)==0 | length(edges$from)==0)
    stop("no edge or term found!")
  
  ##node name and def
  ontnodes <- terms$id
  ontnames <- terms$Term
  nodes=data.frame(name=ontnodes,def=ontnames)
  
  ##create igraph object
  g=graph.data.frame(edges, directed=TRUE, vertices=nodes)
  
  g
}


##################################
## build term levels
##################################
#' @name buildLevels
#' @title buildLevels
#' @aliases buildLevels
#' @description build a set of object to represent the igraph structure
#' @usage 
#' levels<-buildLevels(graph)
#' @details build a set of object to represent the igraph structure
#' @param graph an igraph object.
#' @return a list object contains:
#'     \item{\code{nodes2level}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{level2nodes}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{noOfLevels}:}{Object of class \code{"character"} ~~ }
#'     \item{\code{noOfNodes}:}{Object of class \code{"ANY"} ~~ }
#' @author Xin he
#' @keywords buildLevels
buildLevels<- function(graph){
  nodes2level <- new.env(hash = T, parent = emptyenv())
  queue <- findRoot(graph)
  level <- 1
  adjList <- get.adjlist(graph, mode = "in")
  
  while(length(queue)) {
    ## assign the curent level
    multiassign(as.character(queue), rep(level, length(queue)), envir = nodes2level)
    
    ## move to the next level
    level <- level + 1
    queue <- names(adjList)[unique(unlist(adjList[queue], use.names = TRUE))]
  }
  
  ## revert the node2level
  nl <- unlist(as.list(nodes2level))
  f.index <- rep(sort(unique(nl)), table(nl))
  level2nodes <- split(names(sort(nl)), f.index)
  
  list(nodes2level = as.list.environment(nodes2level),
              level2nodes = level2nodes,
              noOfLevels = level - 1,
              noOfNodes = length(nodes2level))
}


##################################
## give the full ontoGraph and a target level,
## this method return a sub ontoGraph with the nodes above the level while maintain the overall tree structure
## levels is generate by the method buildLevels(). If levels is not provide, it will calaulate levels for the graph.
##################################
#' @name subGraphByLevel
#' @title subGraphByLevel
#' @aliases subGraphByLevel
#' @description build a subGraph with the specific level. Nodes deeper then level l will be removed
#' @usage subGraphByLevel(graph,level)
#' @param graph an igraph object.
#' @param l an integer specify the level you want to keep from root to bottom
#' @return an pruned igraph object
#' @author Xin he
#' @export
#' @keywords subGraphByLevel
subGraphByLevel<-function(graph,l=3){
  levels=buildLevels(graph)
  myterms<-unique(unlist(levels$level2nodes[1:l],use.names=FALSE))
  gsub<-subGraphByNodes(graph,myterms)
  gsub
}


##################################
## give the full ontoGraph and a set of nodes,
## this method return a sub ontoGraph with the nodes of interest while maintain the overall tree structure
##################################
subGraphByNodes<-function(graph,nodes,include_path_to_root=TRUE){
  ifelse(include_path_to_root,mynodes <- findInducedSubGraphNodes(graph,nodes),mynodes<-nodes)
  subvids<-V(graph)[V(graph)$name %in% mynodes]
  gsub<-induced.subgraph(graph,vids=subvids,impl='auto')  
  gsub
}

##################################
## add node attribute color
## color the node base on node level
##################################
node.addColorAttributeByLevel<-function(graph,levels=NULL){
  if(is.null(levels)){
    V(graph)$color='grey'
  }else{
    V(graph)$color=rainbow(levels$noOfLevels)[unlist(levels$nodes2level[V(graph)$name],use.names=FALSE)]
  }
  graph
}

##################################
## add node attribute level
##################################
node.addLevelAttribute<-function(graph,levels){
  V(graph)$level=unlist(levels$nodes2level[V(graph)$name],use.names=FALSE)
  graph
}

##################################
## add node attribute levf
##################################
node.addLeafAttribute<-function(graph){
  leaves<-findLeafNode(graph)
  V(graph)$is_leaf=FALSE
  V(graph)[V(graph)$name %in% leaves]$is_leaf=TRUE
  graph
}


##################################
## add node attribute Def
##################################
node.addDefAttribute<-function(graph){
  def=Definition(ONTTERM)
  key=match(x = names(def),V(graph)$name)
  V(graph)[key]$description=def
  graph
}

##################################
## method to find all leaves node in a graph
## need to be a directed igraph
##################################
findLeafNode<-function(graph){
  if(!is.directed(graph))
    stop("not a directed graph!")
  adj<-get.adjlist(graph,mode='in')
  leaf.list<-lapply(X=adj,function(x) ifelse(length(adj[x])==0,TRUE,FALSE))
  names(leaf.list[leaf.list==TRUE])
}

##################################
## method root node
##################################
findRoot<-function(graph){
  aux <- sapply(get.adjlist(graph,mode='out'), length)
  names(which(aux == 0))
}

##################################
## method to find a node's direct parents
## need to be a directed igraph and direction should be from parents to children
##################################
findParentNodes<-function(graph,node){
  adj<-get.adjlist(graph,mode='out')
  .p<-adj[[node]]
  parents<-V(graph)[.p]$name
  unique(parents)
}

##################################
## find parent nodes recursively till root reached
##################################
findAllParentNodes<-function(graph,node){
  root<-findRoot(graph)
  parents=findParentNodes(graph,node)
  for(p in parents){
    if(!p %in% root){
        ps = findAllParentNodes(graph,p)
        parents=c(parents,ps)
    }
  }
  unique(parents)
}

##################################
## method to find a node's direct children
## need to be a directed igraph and direction should be from parents to children
##################################
findChildrenNodes<-function(graph,node){
  rg=reverseArch(graph)
  findParentNodes(rg,node)
}

##################################
## find children nodes recursively till leaf reached
##################################
findAllChildrenNodes<-function(graph,node){
  rg=reverseArch(graph)
  findAllParentNodes(rg,node)
}


##################################
## find the shorest path to root
##################################
shortest_path_to_root<-function(graph,node){
  .r<-get.shortest.paths(graph,from=V(graph)[node],to=V(graph)['all'],mode='out',output='vpath')
  nodes<-.r[['vpath']][[1]]
  V(graph)[nodes]$name
}


##################################
## find node level
##################################
find.node.level<-function(graph,nodes=c()){
  if(!is.nodes.in.graph(graph,nodes))
    warning("some of the nodes are not in the graph!")
  
  l<-buildLevels(graph)
  l$nodes2level[nodes]
}

##################################
## check if the nodes are in the graph
##################################
is.nodes.in.graph<-function(graph,nodes=c()){
  if(sum(!nodes %in% get.node.attribute(graph,'name'))>0)
    FALSE
  else
    TRUE
}
##################################
## check which the nodes are in the graph
##################################
which.node.in.graph<-function(graph,nodes=c()){
  nodes %in% get.node.attribute(graph,'name')
}

##################################
## check if the nodes are leaf nodes
##################################
is.leaf<-function(graph,nodes){
  ln<-findLeafNode(graph)
  r<-nodes %in% ln
  names(r)<-nodes
  r
}

##################################
## set node color
##################################
set.node.color<-function(graph,nodes=c(),color){
  set.node.attribute(graph,'color',color,nodes=nodes)
}


##################################
## Get nodes attrs
##################################
get.node.attribute<-function(graph,attr,nodes=c()){
  if(length(nodes)==0) 
    nodes=V(graph)$name
  
  x <- nodes %in% V(graph)$name
  if(!all(x)) {
    warning("Nodes not present in the graph:", nNames[!x])
    nodes <- nodes[x]
  }
  
  retValue <- get.vertex.attribute(graph, attr,nodes )
  
  if(length(retValue) == 1)
    return(retValue[[1]])
  
  return(retValue)
}

##################################
## Set nodes attrs
##################################
set.node.attribute<-function(graph,attr_name,attr_value,nodes=c()){
  if(length(nodes)==0) nodes=V(graph)$name
  x <- nodes %in% V(graph)$name
  if(!all(x)) {
    warning("Nodes not present in the graph:", nNames[!x])
    nodes <- nodes[x]
  }
  graph<-set.vertex.attribute(graph, attr_name, nodes, attr_value)
  graph
}

##################################
## reverse the direction of the graph
## tree plot layout 'layout.reingold.tilfor' require direction to be from children to parents
##################################
reverseArch <- function(graph) {
  if(!is.directed(graph))
    stop("not a directed graph!")
  
  nodes<-get.data.frame(graph,what='vertices')
  edges<-get.data.frame(graph,what='edges')
  
  edges.rev<-data.frame(edges[,2],edges[,1],edges[,3:ncol(edges)])
  colnames(edges.rev)=colnames(edges)
  
  g=graph.data.frame(edges.rev, directed=TRUE, vertices=nodes)
  g
}

############################   findInducedSubGraphNodes   ############################
## Given a term (or a list of terms) this function is returning
## the nodes in the subgraph induced by startNodes
findInducedSubGraphNodes <- function(graph,startNodes) {
  
  ## build a lookUp table with the nodes in the graph
  nodeLookUp <- new.env(hash = T, parent = emptyenv())
  
  nodesGraph <- V(graph)$name
  
  adjList <- get.adjlist(graph, mode = "out")
  
  ## recursivly build the list of induced nodes
  buildInducedGraph <- function(node) {
    ## if we have visited the node, there is nothing to do
    if(exists(node, envir = nodeLookUp, mode = 'logical', inherits = FALSE))
      return(1)
    
    ## we put the node in the graph and we get his parents
    assign(node, TRUE, envir = nodeLookUp)
    
    adjNodes <- names(adjList)[unlist(adjList[node])]
    
    if(length(adjNodes) == 0)
      return(2)
    
    for(i in 1:length(adjNodes))
      buildInducedGraph(adjNodes[i])
    
    return(0)
  }
  
  ## we start from the specified nodes
  lapply(startNodes, buildInducedGraph)
  return(ls(nodeLookUp))
}



############################   heat color node base on some value   ############################
## 
colorMapNode<-function(graph,nodes,values){
  values<-sub(pattern='< 1e-30',replacement='1e-30',x=values)
  values<-as.numeric(values)
  names(values)=nodes
  
  x=sort(values+10^-20)
  
  log.x=log10(x)
  color <- round(log.x - range(log.x)[1] + 1)
  colorMap <- heat.colors(max(color))
  
  color<-sapply(names(color),function(x){
    colorMap[color[x]]
  })
  
  graph.result=set.node.attribute(graph,attr_name='color',attr_value=color,nodes=names(color))
  graph.result
}



############################   Debug function   ############################
peekNode<-function(graph,node){
  node<-'DOID:10652'
  attrs<-list.vertex.attributes(graph)
  for(attr in attrs){
    cat(paste(attr,get.node.attribute(graph,attr,node),sep=":"))
    cat("\n\n")
  }
}


############################   Debug function   ############################
to.latex<-function(topOntoresult,table.name='',table.des='',number_of_row=5){
  #head(T1$GOBP_P)
  #t<-T1$GOBP_P[1:5,]
  
  
  
  
  if(number_of_row>nrow(topOntoresult))
    number_of_row=nrow(topOntoresult)
  
  t<-topOntoresult[1:number_of_row,c(1,2,4,5,6,8,9)]
  colnames(t)<-c("TERM.ID","TERM.NAME","Annotated","Significant","Expected","elim-p","p")
  cat(s)
  cat("\\begin{table}[ht]\n")
  cat("\\rowcolors{2}{gray!25}{white}\n")
  cat("{\\footnotesize")
  cat("\n")
  cat("\\begin{tabular}{*{2}{l}*{5}{c}}\n")
  cat(paste(colnames(t),collapse="&"))
  cat("\\\\")
  cat("\n")
  cat("\\hline")
  cat("\n")
  for(i in 1:nrow(t)){
    cat(paste(t[i,],collapse="&"))
    cat("\\\\")
    cat("\n")
  }
  cat("\\end{tabular}\n")
  cat("}")
  cat(paste("\\caption{",table.des,"}\n",sep=""))
  cat(paste("\\label{tab:",table.name,"}\n",sep=""))
  cat("\\end{table}\n")
}

to.latex.content<-function(topOntoresult,number_of_row=5){
 
  if(number_of_row>nrow(topOntoresult))
    number_of_row=nrow(topOntoresult)
  
  t<-topOntoresult[1:number_of_row,c(1,2,4,5,6,8,9)]
 
  for(i in 1:nrow(t)){
    cat(paste(i,paste(t[i,],collapse="&"),sep="AAA"))
    cat("\\\\")
    cat("\n")
  }

}




# to.latex(T1$GOBP_P,5)
# topOntoresult=T1$GOBP_P
# setwd('/home/xin/Desktop/usercase/fly/analysis/')
# 
# setwd('atowtt1pVSn')
# T1<-new.env(hash = TRUE)
# p.n<-readRDS('foldchange.rds')
# p.n.enriched<-p.n[p.n$fc>1.5 & p.n$fdr<0.01,]
# ##332(paper 330)
# length(unique(p.n.enriched$entrez_id))
# assign('p.n',p.n,T1)
# assign('p.n.enriched',p.n.enriched,T1)
# 
# GOBP_P<-readRDS('topOnto/GOBP_T1P.rds')
# assign('GOBP_P',GOBP_P[as.numeric(GOBP_P$elimfisher)<0.05,],T1)
# GOBP_N<-readRDS('topOnto/GOBP_T1N.rds')
# assign('GOBP_N',GOBP_N[as.numeric(GOBP_N$elimfisher)<0.05,],T1)
# HDO_P<-readRDS('topOnto/HDO_T1P.rds')
# assign('HDO_P',HDO_P[as.numeric(HDO_P$elimfisher)<0.05,],T1)
# HDO_N<-readRDS('topOnto/HDO_T1N.rds')
# assign('HDO_N',HDO_N[as.numeric(HDO_N$elimfisher)<0.05,],T1)
# HPO_P<-readRDS('topOnto/HPO_T1P.rds')
# assign('HPO_P',HPO_P[as.numeric(HPO_P$elimfisher)<0.05,],T1)
# HPO_N<-readRDS('topOnto/HPO_T1N.rds')
# assign('HPO_N',HPO_N[as.numeric(HPO_N$elimfisher)<0.05,],T1)
# setwd('/home/xin/Desktop/usercase/fly/analysis/')
# 

