
##################################
## brirf information of the graph
analysis<-function(g){
  levels<-buildLevels(g)
  cat("# No of nodes\t")
  print(length(V(g)))
  cat("# No of edges\t")
  print(length(E(g)))
  cat("# max deep\t")
  print(levels$noOfLevels)
  
  cat("# Density (No of edges / possible edges)\t")
  print(graph.density(g))
  cat("# Number of islands\t")
  print(clusters(g)$no)
  cat("# Diameter of the graph (the length of the longest geodesic)\t")
  print(diameter(g))
  
  par(mfrow=c(2,2))
  #plot
  #degree.distribution(g)
  plot(degree.distribution(g), xlab="node degree")
  lines(degree.distribution(g))
  #if(length(V(g)<3000)) nomalplot(g,label=0)
  
  leveldis<-unlist(levels$nodes2level,use.names=FALSE)
  plot(table(leveldis))
  plot(prop.table(table(leveldis)), xlab="deep",ylab='%')
  
  #degree(graph)
  #plot(degree.distribution(graph, cumulative = FALSE) )
  #shortest.paths(g)
  par(mfrow=c(1,1))
}

###############################################
# Define link between nodes base on there annotation
# lookup should should be like this:
# List of 231
# $ all         : NULL
# $ DOID:0050432: chr [1:111] "18" "367" "412" "414" ...
# $ DOID:0050587: chr [1:3] "3218" "58512" "114798"
# ...
findlink<-function(id1=1,id2=2,lookup){
  out=list('weight'=0,'genes'=c())
  #print(paste(id1,id2,sep=" vs "))
  
  if(id1!=id2){
    t1=lookup[[id1]]
    t2=lookup[[id2]]
    if(length(t1)>0 & length(t2)>0){
      share=intersect(ls(t1),ls(t2))
      p1<-length(share)/length(t1)
      p2<-length(share)/length(t2)
      #w <- mean(c(p1,p2))
      out[['p1']]= length(t1)
      out[['p2']]= length(t2)
      out[['share']]= length(share)
      out[['weight']]= round(x=p1*p2,digits=5)
      out[['genes']] = share
    }    
  }
  out
}



###############################################
## Ignore weight, if two nodes share genes, return TRUE
findBiLink<-function(id1=1,id2=2,lookup){
  out=findlink(id1,id2,lookup)
  out=ifelse(length(out)>0,TRUE,FALSE)
  out
}

###############################################
## calculate adj matrix base on annotation data
## if two node has shared annotation, a link is put between these two nodes
##the weight is compute by the method findlink()
calculate.upper.adj.matrix<-function(graph,leaf_only=FALSE){
  ##generate lookup table
  ifelse(leaf_only, 
         names<-findLeafNode(graph), 
         names<-get.node.attribute(graph,attr='name'))
  lookup=get.node.attribute(graph,attr='genes')
  names(lookup)<-names
  
  ##generate the adjacency matrix
  m<-matrix(c(0),nrow=length(names),ncol=length(names),dimnames=list(names,names))
  for(i in 1:length(names)){
    for(j in i:length(names)){
      r=findlink(names[i],names[j],lookup)
      if(r[['weight']]>0)
        m[i,j]<-r[['weight']]
        #print(paste(names[i],names[j],r[['weight']],paste(r[['genes']],collapse=','),sep="\t"))
    }
  }
  m
}

calculate.upper.adj.matrix.print<-function(graph,leaf_only=FALSE){
  ##generate lookup table
  ifelse(leaf_only, 
         names<-findLeafNode(graph), 
         names<-get.node.attribute(graph,attr='name'))
  lookup=get.node.attribute(graph,attr='genes')
  names(lookup)<-names
  
  ##generate the adjacency matrix
  m<-matrix(c(0),nrow=length(names),ncol=length(names),dimnames=list(names,names))
  for(i in 1:length(names)){
    message(paste(i,length(names),sep='/'))
    for(j in i:length(names)){
      r=findlink(names[i],names[j],lookup)
      if(r[['weight']]>0)
        #m[i,j]<-r[['weight']]
        cat(paste(names[i],names[j],r[['p1']],r[['p2']],r[['share']],r[['weight']],paste(r[['genes']],collapse=','),sep="\t"),"\n")
    }
  }
  #m
}
###############################################
## calculate adj matrix base on annotation data
## if two node has shared annotation, a link is put between these two nodes
##the weight is compute by the method findlink()
calculate.upper.adj.matrix3<-function(graph,leaf_only=FALSE){
  ##generate lookup table
  ifelse(leaf_only, 
         names<-findLeafNode(graph), 
         names<-get.node.attribute(graph,attr='name'))
  lookup=get.node.attribute(graph,attr='genes')
  names(lookup)<-names
  
  x=new.env(hash = TRUE, parent = emptyenv())    
  for(i in 1:length(names)){
    id1=new.env(hash = TRUE, parent = emptyenv())
    assign(names[i], id1, envir = x)
    for(j in i:length(names)){
      id2=new.env(hash = TRUE, parent = emptyenv())
      
      r=findlink(names[i],names[j],lookup)
      assign('weight', r[['weight']], envir = id2)
      e_gene <- new.env(hash = TRUE, parent = emptyenv())
      if(length(r[['genes']])>0)
        multiassign(r[['genes']], rep(TRUE, length(r[['genes']])), envir = e_gene)    
      assign('genes', e_gene, envir = id2)
      
      assign(names[j], id2, envir = id1)
    }
  }
  x
}

findResult<-function(r,id1,id2){
  k1=as.list.environment(r)[[id1]]
  k2=as.list.environment(r)[[id2]]
  if(exists(x=id2,envir=k1)){
    retult=as.list.environment(get(x=id2,envir=k1))
  }else{
    retult=as.list.environment(get(x=id1,envir=k2))
  }
  retult
}

calculate.upper.adj.matrix2<-function(graph,leaf_only=FALSE){
  ##generate lookup table
  ifelse(leaf_only, 
         names<-findLeafNode(graph), 
         names<-V(graph)$name)
  lookup <- new.env(hash = T, parent = emptyenv())
  multiassign(x=names,value=V(graph)[names]$genes,envir=lookup)
  get(x='DOID:0060046',envir=lookup)
  
  ##generate the adjacency matrix
  m<-matrix(c(0),nrow=length(names),ncol=length(names),dimnames=list(names,names))
  for(i in 1:length(names)){
    for(j in i:length(names)){
      r=findlink2(names[i],names[j],lookup)
      if(r[['weight']]>0)
        m[i,j]<-r[['weight']]
    }
  }
  
  m
}