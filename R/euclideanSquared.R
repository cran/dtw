## Euclidean squared distance.
## should work on vectors (in R^n) too

## `euclideanSquared` <-
## function(a,b) {
##   z<-a-b;
##   z<-drop(z %*% z);                     #inner dot product
##   return (z);
## }

`euclideanSquared` <-
  function(a,b) {
    return((a-b)^2);
  }
