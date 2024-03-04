td = function(mat){return(as.data.frame(t(mat)))}

chunks = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

nbr = function(x){return(length(x))}

nbrs = function(x){return(lengths(x))}

nbrwhich = function(x){return(length(which(x)))}

nbrintersect = function(v1,v2){return(length(intersect(v1,v2)))}

nbrintersect_all <- function(a,b,...){return(length(Reduce(intersect, list(a,b,...))))}

nbrsetdiff = function(v1,v2){return(length(setdiff(v1,v2)))}

nbrunique = function(vector){return(length(unique(vector)))}

nbrunlist = function(vector){return(length(unlist(vector)))}

add = function(x, i){
  if(is.atomic(x)){eval.parent(substitute(x <- c(x,i)))}
  if(is.list(x))  {eval.parent(substitute(x <- append(x,list(i))))}
}

intersect_list = function(a_list, freq_min = NULL){
  if(is.null(freq_min)){
    freq_min = length(a_list)
  }
  a_list = unlist(a_list)
  list_factor = as.factor(a_list)
  result = unique(a_list[which(list_factor %in% names(which(table(list_factor) >= freq_min)))])
  if(length(result) == 0){
    return(NULL)
  }
  return(result)
}

intersect_multi = function(..., y){
  alist = find_object(list(...), test.function = is.atomic)
  sapply(alist, function(x) intersect(x, y))
}

vnames = function(..., type = NULL, negate = FALSE){
  results = lapply(find_igraphs(list(...)), function(x) igraph::V(x)$name)
  if(is.null(type)){
    if(length(results) == 1){
      return(results[[1]])
    } else {
      return(results)
    }
  } else if(is.character(type)){
    if(length(results) == 1){
      return(stringr::str_subset(string = results[[1]], pattern = type, negate = negate))
    } else {
      return(lapply(results, function(x) stringr::str_subset(string = x, pattern = type, negate = negate)))
    }
  }
}
