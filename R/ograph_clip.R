##################################
##
##
searchDescription4Keyword<-function(graph,keys){
  .des=get.node.attribute(graph,'description')
  pattern=paste(keys,collapse='|')
  hits=grep(pattern, .des, perl=TRUE, value=FALSE)
  terms=c()
  if(length(hits))
    terms=V(graph)[hits]$name
  
  terms
}

searchName4Keyword<-function(graph,keys){
  .des=get.node.attribute(graph,'def')
  pattern=paste(keys,collapse='|')
  hits=grep(pattern, .des, perl=TRUE, value=FALSE)
  terms=c()
  if(length(hits))
    terms=V(graph)[hits]$name
  
  terms
}





#hits=searchDescription4Keyword(g@graph,c('neuro','brain'))
#treeplot(subGraphByNodes(g@graph,nodes=hits))
