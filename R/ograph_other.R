addFCurve<-function(){
for(f in seq(0.1,1,0.1)) {
  curve(f*x/(2*x-f), f/2,1.1, add=TRUE, n=1000, col="lightgrey",)
  text(1,f/(2-f), label=sprintf("f=%.1f",f), pos=1, col="grey")
}
}

calculatePRF<-function(true,pred,filter=FALSE,g=NULL,plot=FALSE){
  if(filter){
    
    if(is.null(g))
      stop("a valide ontgraph object is needed.")
    
    if(length(pred)>0 & pred[1] != ""){
      #filterout obsoleted terms
      true=true[true %in% V(g@graph)$name]
      pred=pred[pred %in% V(g@graph)$name]
      share=intersect(true,pred)
      mg=subGraphByNodes(graph=g@graph,nodes=unique(c(true,pred)))
      
      ##mainly for debugging
      if(plot){
        if(length(true)>2 && length(pred)>2){
          V(mg)$color='grey'
          V(mg)$color[V(mg)$name %in% true]='yellow'
          V(mg)$color[V(mg)$name %in% pred]='red'
          V(mg)$color[V(mg)$name %in% share ]='green'
          treeplot(mg,vertex.size=20,label=1,vertex.label.cex=0.9)
          legend("bottomright",c('true','pred','hit'),fill=c("yellow","red","green"),bty='n',cex = 0.75)
        }
      }
      
      
      pred=sapply(pred,function(x){
        ps=findAllParentNodes(mg,x)
        t=true[true %in% ps]
        ifelse(length(t)>0,t[1],x)   
      })
      
      pred=unique(pred)
    }else{
      ##pred is empty so no filter is needed.
    }
  }
  
  tp <- pred[pred %in% true]
  fp <- pred[!pred %in% true]
  fn <- true[!true %in% pred]
  
  precision <- length(tp) / length(pred)
  recall <- length(tp)/length(true)
  Fmeasure <- 2 * precision * recall / (precision + recall)
  
  if(is.na(Fmeasure))
    Fmeasure = 0
  
  if(is.na(precision))
    precision = 0
  
  if(is.na(recall))
    recall = 0 
  
  r=c(round(precision,3),round(recall,3),round(Fmeasure,3))
  names(r)<-c("p","r","f")
  r
}