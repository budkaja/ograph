##################################
##simple method to plot the tree with pre-defined parms
treeplot<-function(graph,label=TRUE,label.nodes=c(),
                   vertex.size=3,
                   vertex.label.cex=1,
                   edge.arrow.size=0.2,
                   edge.width=0.5,
                   vertex.label.dist=0,
                   vertex.label.degree=-pi/4,
                   show.genes=FALSE,
                   only.gene=FALSE,
                   root='all'){
  graph<-reverseArch(graph)
  
  #text=paste(V(graph)$name,V(graph)$def,sep="\n")
  label.text=vector(mode='character',length(V(graph)))

  if(label){
    if(length(label.nodes)>0){
      #always plot the first three level
      default.nodes<-V(subGraphByLevel(reverseArch(graph),3))$name
      label.nodes=unique(c(default.nodes,label.nodes))
      index<-which(V(graph)$name %in% label.nodes)
      label.text[index]=paste(V(graph)$name[index],V(graph)$def[index],sep="\n")
    }else{
      label.text=paste(V(graph)$name,V(graph)$def,sep="\n")
    }
    
    
    if(show.genes){
      if(only.gene)
        label.text=sapply(V(graph)$genes,length)
      else
        label.text=paste(label.text,sapply(V(graph)$genes,length),sep="\n")
    }
    
  }

  
  if(label){
    #plot(graph,vertex.size=vertex.size,vertex.label.cex=vertex.label.cex,vertex.label=paste(V(graph)$name,V(graph)$name,sapply(V(graph)$genes,length),sep="\n"),edge.arrow.size=edge.arrow.size,edge.width=edge.width,edge.color='black',layout=layout.reingold.tilford(graph,flip.y=TRUE,root=which(V(graph)$name=='all')))
    plot(graph,vertex.size=vertex.size,vertex.label.dist=vertex.label.dist,vertex.label.degree=vertex.label.degree,
         vertex.label.cex=vertex.label.cex,
         vertex.label=label.text,
         edge.arrow.size=edge.arrow.size,
         edge.width=edge.width,
         edge.color='black',
         layout=layout.reingold.tilford(graph,flip.y=TRUE,root=which(V(graph)$name==root)))
  }else{
    plot(graph,vertex.size=vertex.size,vertex.label=NA,edge.arrow.size=edge.arrow.size,edge.width=edge.width,edge.color='black',layout=layout.reingold.tilford(graph,flip.y=TRUE,root=which(V(graph)$name=='all')))
  }
}

##similar to tree plot, but use interactive plotting. work only with small graph.
##
tktreeplot<-function(graph,label=0,vertex.size=3,vertex.label.cex=1,edge.arrow.size=0.2,edge.width=0.5,vertex.label.dist=0,vertex.label.degree=-pi/4,show.genes=FALSE,only.gene=FALSE,root='all'){
  if(length(V(graph))>500){
    die('too many nodes!')
  }
  graph<-reverseArch(graph)
  
  text=paste(V(graph)$name,V(graph)$def,sep="\n")
  if(show.genes){
    if(only.gene)
      text=sapply(V(graph)$genes,length)
    else
      text=paste(text,sapply(V(graph)$genes,length),sep="\n")
  }
  
  if(label==1){
    #plot(graph,vertex.size=vertex.size,vertex.label.cex=vertex.label.cex,vertex.label=paste(V(graph)$name,V(graph)$name,sapply(V(graph)$genes,length),sep="\n"),edge.arrow.size=edge.arrow.size,edge.width=edge.width,edge.color='black',layout=layout.reingold.tilford(graph,flip.y=TRUE,root=which(V(graph)$name=='all')))
    tkplot(graph,vertex.size=vertex.size,vertex.label.dist=vertex.label.dist,vertex.label.degree=vertex.label.degree,
         vertex.label.cex=vertex.label.cex,
         vertex.label=text,
         edge.arrow.size=edge.arrow.size,
         edge.width=edge.width,
         edge.color='black',
         layout=layout.reingold.tilford(graph,flip.y=TRUE,root=which(V(graph)$name==root)))
  }else{
    tkplot(graph,vertex.size=vertex.size,vertex.label=NA,edge.arrow.size=edge.arrow.size,edge.width=edge.width,edge.color='black',layout=layout.reingold.tilford(graph,flip.y=TRUE,root=which(V(graph)$name=='all')))
  }
}


##################################
##simple method to plot the tree with pre-defined parms
nomalplot<-function(graph,label=0){
  if(label==1){
    plot(graph,vertex.size=3,vertex.label.cex=1,vertex.label=paste(V(graph)$name,V(graph)$name,sapply(V(graph)$genes,length),sep="\n"),edge.arrow.size=0.2,edge.width=0.5,edge.color='black',layout=layout.fruchterman.reingold)
  }else{
    plot(graph,vertex.size=3,vertex.label=NA,edge.arrow.size=0.2,edge.width=0.5,edge.color='black',layout=layout.fruchterman.reingold)
  }
}


##################################
##value is a named vector with nodes and their values.
##plot the tree and color the significant nodes
##plot2file(filename,width=50,heigth=20)
##plotSig(graph=g@graph,testresult=resultElimFis,number_of_node=50,label=1)
##dev.off()
plotSig<-function(graph,value,number_of_node=0,only.plot.sig=T,...){
  #turn to numeric
  tmp<-names(value)
  value<-sub(pattern='< 1e-30',replacement='1e-30',x=value)
  value<-as.numeric(value)
  names(value)<-tmp
  
  x=sort(value+10^-20)
  #x=sort(score(testresult))
  
  if(length(x)>number_of_node & number_of_node>0 )
    x=x[1:number_of_node]
  
  log.x=log10(x)
  
  color <- round(log.x - range(log.x)[1] + 1,3)
  index=unique(color)
  colorMap <- heat.colors(length(index))
  
  color<-sapply(names(color),function(x){
    colorMap[which(index==color[x])]
  })
  
  if(!exists('label.nodes'))
    label.nodes=names(x)
  
  g=subGraphByNodes(graph,nodes=c(names(x),label.nodes))
  g=ograph::set.node.attribute(g,attr_name='color',attr_value=color,nodes=names(color))
  
  if(only.plot.sig)
    treeplot(g,label.nodes=names(x),...)
  else
    treeplot(g,...) 
  
}


##################################
##function to plot the following graph into a file
##need to call dev.off when finish plotting
plot2file<-function(filename,width=12,heigth=8,units='in',res=300){
  png(filename, width,heigth, units=units, res=res)
}

##################################
##plot the wordcloud base on p-value
##value is a named vector with nodes and their values.
plotWordcloud<-function(value,number_of_node=Inf,scale=c(3,0.1),filename='',width=12,heigth=8,units='in',res=300){
  require(wordcloud)
  def=Term(ONTTERM)
  
  ns<-names(value)
  value<-sub(pattern='< 1e-30',replacement='1e-30',x=value)
  value<-as.numeric(value)
  names(value)<-ns
  x=sort(value+10^-20)
  
  if(!is.infinite(number_of_node)){
    x=x[1:number_of_node]
  }
  
  y=-log(x)
  freq=y/sum(y)
  min.freq=sort(freq[freq>0])[1]
  if(filename!=''){
    png(filename, width,heigth, units=units, res=res)
    wordcloud(words=def[names(y)],freq=freq,scale=scale,min.freq=min.freq,random.order=FALSE, max.words=Inf,rot.per=0, use.r.layout=FALSE, colors=brewer.pal(8, 'Dark2'))
    dev.off() 
  }else{
    wordcloud(words=def[names(y)],freq=freq,scale=scale,min.freq=min.freq,random.order=FALSE, max.words=Inf,rot.per=0, use.r.layout=FALSE, colors=brewer.pal(8, 'Dark2'))
  }
}

###############################################
## save igraph object to format
saveGraph<-function(graph,name,loc,graphml=TRUE,e=TRUE,v=TRUE){
  
  if(graphml)
    write.graph(graph,paste(loc,name,".graphml",sep=''), format="graphml")
  
  tmp=get.data.frame(graph,what='both')
  if(e)
    write.table(tmp$edges, file = paste(loc,name,"_edges.txt",sep=''), sep = "\t",row.names = FALSE,quote=FALSE)
  if(v){
    ##parse the genes coloum
    tmp$vertices$genes<-.listOfEnv2ListOfList(tmp$vertices$genes)
    if(length(tmp$vertices$genes)>0){
    tmp$vertices$genes<-sapply(tmp$vertices$genes,FUN=function(x){
     paste(x,collapse=',')
    })
    }
    write.table(tmp$vertices, file = paste(loc,name,"_nodes.txt",sep=''), sep = "\t",row.names = FALSE,quote=FALSE)
  }
  
}

loadGraph<-function(edges_file,nodes_file){
  edges=read.table(edges_file,header=TRUE,sep='\t')
  nodes=read.table(nodes_file,header=TRUE,sep='\t')
  
  unique(union(edges$from,edges$to))
  nodes=nodes[,1:2]
  levels(nodes$name)
  
  setdiff(unique(union(edges$from,edges$to)),levels(nodes$name))
  
  g=graph.data.frame(edges, directed=FALSE, vertices=nodes)
}

.listOfEnv2ListOfList<-function(listOfEnv){
  lapply(listOfEnv,ls)
}


##################################
##plot the graph in text format 
##giving an ider how the structure looks like
plotGraphStructure<-function(graph,indent='--',text=c('name')){
  root<-findRoot(graph)
  levels<-buildLevels(graph)
  
  f<-function(graph,node){
    level<-levels$nodes2level[[node]]
    string=paste(rep(indent,level),collapse='')
    t=''
    for(i in text){
      t = c(t,get.node.attribute(graph,i,c(node)))
    }
    t=paste(t,collapse=' ')
    cat(paste(string,t,"\n",sep=' '))
    cs<-findChildrenNodes(graph,node)
    if(length(cs)>0){
      for(i in cs){
        f(graph,i)
      }
    }
  }
  
  f(graph,root)
}



################
#method to turn an igraph to graphNEL
.to.GraphNEL<-function(igraph){
  
  nel<-igraph.to.graphNEL(igraph)
  nAttrs<-list()
  
  v.n <- list.vertex.attributes(igraph)
  v.n <- v.n[v.n != "name"]
  index<-get.vertex.attribute(igraph, 'name')
  
  for (n in v.n) {
    nAttrs[[n]]<-unlist(nodeData(nel, attr = n))
  }
  
  
  dic<-c('color'='fillcolor')
  names(nAttrs)<-unname(sapply(names(nAttrs),function(x){
    if(is.na(dic[x]))
      x
    else
      dic[x]
  }))
  
  list(graph=nel,nodeAttrs=nAttrs)
}


################
#method to turn an igraph to graphNEL and plot
plot.graphNEL<-function(igraph,label=FALSE,showEdges = TRUE,node.shape='circle',node.fontsize = 9,edge.fontsize = 9,node.height = 0.45,label.only.def=T){
  require(Rgraphviz)
  r<-.to.GraphNEL(igraph)
  nel<-r$graph
  nodeAttrs<-r$nodeAttrs
   
  
  node.names <- nodes(nel)
  if(label.only.def==1)
    nodeAttrs$label <- .getTermsDefinition(igraph,node.names,multipLines=T)
  else
    nodeAttrs$label <- paste(node.names,nodeAttrs$def, sep = "\\\n")
  names(nodeAttrs$label) <- node.names
  
#   nodeAttrs$shape<-'circle'
#   nodeAttrs$height<-0.45
  
  ## we set the global Graphviz attributes
  graphAttrs <- getDefaultAttrs(layoutType = 'dot')
  graphAttrs$cluster <- NULL
  
  #graphAttrs$graph$splines <- FALSE
  graphAttrs$graph$size <- "6.99,3.99"
  
  ## set the node shape
#   graphAttrs$node$shape <- 'ellipse'
  graphAttrs$node$shape <- node.shape

  
  ## set the fontsize for the nodes labels
  graphAttrs$node$fontsize <- node.fontsize
  graphAttrs$edge$fontsize <- edge.fontsize
  graphAttrs$node$style <- 'filled'
  graphAttrs$node$height <- node.height
#   graphAttrs$node$width <- '1.5'
  
  if(!showEdges)
    graphAttrs$edge$color <- 'white'
  else
    ## if we want to differentiate between 'part-of' and 'is-a' edges
    ##    0 for a is_a relation,  1 for a part_of relation
    ##  edgeAttrs$color <- ifelse(.getEdgeWeights(dag) == 0, 'black', 'red')
    graphAttrs$edge$color <- 'black'
  
  
  plot(reverseEdgeDirections(nel),nodeAttrs = nodeAttrs,attrs=graphAttrs)
}



