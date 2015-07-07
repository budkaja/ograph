##################################
##simple method to plot the tree with pre-defined parms
treeplot<-function(graph,label=0,vertex.size=3,vertex.label.cex=1,edge.arrow.size=0.2,edge.width=0.5,vertex.label.dist=0,vertex.label.degree=-pi/4,show.genes=FALSE,only.gene=FALSE,root='all'){
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
    plot(graph,vertex.size=vertex.size,vertex.label.dist=vertex.label.dist,vertex.label.degree=vertex.label.degree,
         vertex.label.cex=vertex.label.cex,
         vertex.label=text,
         edge.arrow.size=edge.arrow.size,
         edge.width=edge.width,
         edge.color='black',
         layout=layout.reingold.tilford(graph,flip.y=TRUE,root=which(V(graph)$name==root)))
  }else{
    plot(graph,vertex.size=vertex.size,vertex.label=NA,edge.arrow.size=edge.arrow.size,edge.width=edge.width,edge.color='black',layout=layout.reingold.tilford(graph,flip.y=TRUE,root=which(V(graph)$name=='all')))
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
plotSig<-function(graph,value,label=1,number_of_node=0,...){
  x=sort(testresult+0.00001)
  #x=sort(score(testresult))
  
  if(length(x)>number_of_node & number_of_node>0 )
    x=x[1:number_of_node]
  
  log.x=log10(x)
  
  color <- round(log.x - range(log.x)[1] + 1)
  colorMap <- heat.colors(max(color))
  
  color<-sapply(names(color),function(x){
    colorMap[color[x]]
  })
  
  g=subGraphByNodes(graph,nodes=names(x))
  g=ograph::set.node.attribute(g,attr_name='color',attr_value=color,nodes=names(color))
  treeplot(g,label=label,vertex.size=5,vertex.label.cex=1,...)
}

##################################
##function to plot the following graph into a file
##need to call dev.off when finish plotting
plot2file<-function(filename,width=12,heigth=8,units='in',res=300){
  png(filename, width,heigth, units=units, res=res)
}

##################################
##give in the result generate by topOnto statics test
##plot the wordcloud base on p-value
plotWordcloud<-function(testresult,number_of_node=Inf,scale=c(3,0.1),filename,width=12,heigth=8,units='in',res=300){
  require(wordcloud)
  def=Term(ONTTERM)
  x=sort(score(testresult))
  y=-log(x)
  freq=y/sum(y)
  min.freq=sort(freq[freq>0])[1]
  png(filename, width,heigth, units=units, res=res)
  wordcloud(words=def[names(y)],freq=freq,scale=c(3,0.1),min.freq=min.freq,random.order=FALSE, max.words=Inf,rot.per=0, use.r.layout=FALSE, colors=brewer.pal(8, 'Dark2'))
  dev.off() 
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