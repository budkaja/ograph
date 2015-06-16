init <- function(){
  cat("ograph loaded.\n")
  #cat("Please run initOGraph(ontology_name) to build the ograph object .\n")
  initWHAT()
}

initWHAT <- function(){
  cat("These ontology package(s) are currently available in your libPath:\n")
  x=list.files(path=.libPaths(),pattern='topOnto.*.db')
  cat(sapply(x,function(y){ y=sub('topOnto.','',y,perl=TRUE,);y=sub('.db','',y,perl=TRUE,)},USE.NAMES = FALSE),"\n") 
  #x=list.files(system.file("extdata", package ="topOnto"),pattern='*.sqlite')
}

initOGraph <- function(ontology='HDO'){
  require("methods", quietly=TRUE)
  pkg <- paste('topOnto',ontology,'db',sep='.')
  
  
  ##detach other topOnto.xx.db package to avoid name conflict
  if( length(grep(pkg , search(), perl=TRUE, value=FALSE)) == 0 ){
    while(length(grep("topOnto.\\w+.db", search(), perl=TRUE, value=FALSE)) >0 ){
      detach(pos = grep("topOnto.\\w+.db", search(), perl=TRUE, value=FALSE)[1], unload=TRUE,force=TRUE)
    }
  }
  
  require(pkg, character.only = TRUE) || stop(paste('package ',pkg,' is required',sep=''))
  require("igraph")
  #detach(pos = match(paste("package", pkg, sep = ":"), search()))
}