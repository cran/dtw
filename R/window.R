
## no warping window: no restrictions

`noWindow` <-
function(iii,jjj) {
  return(T);
}



# A band around the diagonal. The band includes the diagonal +-
# window.size, measured on the second index

`sakoeChibaWindow` <-
function(iii,jjj) {
  diagj<-(iii*template.size/query.size);
  return(abs(jjj-diagj)<=window.size);
}


## DUMMY!! TODO
## itakura parallelogram: to be done

`itakuraWindow` <-
function(iii,jjj) {
  return(T);
}
