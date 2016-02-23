
.onAttach <- function(lib, pkg) {
  # where <- match(paste("package:", pkg, sep=""), search())
  # initVar(where) 
  ## some preprocessing
  init()
  if("igraph" %in% rownames(installed.packages())){
    v<-packageVersion("igraph")
    if(!v=='0.7.1'){
      warning("igraph 0.7.1 is required!, Please install it by typing the following into your R prot!
              install.packages('https://cran.r-project.org/src/contrib/Archive/igraph/igraph_0.7.1.tar.gz', repos=NULL, type='source')")
    }
  }else{
    warning("igraph 0.7.1 is required!, Please install it by typing the following into your R prot!
              install.packages('https://cran.r-project.org/src/contrib/Archive/igraph/igraph_0.7.1.tar.gz', repos=NULL, type='source')")
  }
  
  
  
}
packageVersion("topOnto.HPO.db")
