##' A Simple Wrapper to Run Model
##'
##' content for details here 
##'
##' 
##' @title runmodel
##' @param times vector of times to run the model over
##' @param p parameter list
##' @return a matrix of timeseries
##' @author Pete Dodd
##' @export
runmodel <- function(times,p){
  mdl <- tbmod0$new(user=p)
  mdl$run(times)
}
