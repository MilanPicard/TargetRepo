#############################################
#  extract_by_rwr
extract_by_rwr = function(Graph, start_nodes, directionel = TRUE, weight = FALSE, nCores = 1, restart = .5, adj_method = get.adjacency_matrix, processed = TRUE){
  # Indirect neighboordhood from starting nodes about all the graph
  # create output dataframe
  if(nbrunique(start_nodes) != nbr(start_nodes)){
    stop("start_nodes must have only unique values")
  }
  
  # Prepare adjacency matrix for later computation
  if(isFALSE(weight)){
    ADJ = adj_method(Graph, attr = NULL)
  } else {
    ADJ = adj_method(Graph, attr = weight)
  }
  if(nCores == 1){
    result_list = list()
    for(i in 1:nbr(start_nodes)){
      if(i%%1 == 0){
        cat('\r', paste0(i,'/',nbr(start_nodes), " seeds"))
      }
      res = RWR(Graph = Graph, Seeds = start_nodes[i], adjacency = ADJ)
      add(result_list, res$Score[match(vnames(Graph), res$Node)])
    }
    df = do.call(rbind, result_list)
    df = as.data.frame(df)
  } else {
    df = as.data.frame(matrix(ncol = vcount(Graph), nrow = nbr(start_nodes)))
    results = parallel::mclapply(start_nodes, function(seed) RWR(Graph = Graph, Seeds = seed, adjacency = ADJ, r = restart), mc.cores = nCores)
    for(i in 1:nbr(results)){
      df[i, ] = results[[i]]$Score[match(vnames(Graph), results[[i]]$Node)]
    }
  }
  rownames(df) = start_nodes
  colnames(df) = vnames(Graph)
  # Remove seeds from the columns features
  df[is.na(df)] = restart
  message("start post-processing")
  if(processed){
    df = log10(df)
    df = as.data.frame(scales::rescale(as.matrix(df), c(0,1)))
    df[is.infinite(as.matrix(df))] = 0
    df = signif(df, 2)
  } else {
    df = signif(df, 2)
  }
  df = cbind(Target = rownames(df), df)
  rownames(df) = NULL
  return(df)
}

get.adjacency_matrix = function (graph, attr = NULL){
  adjacency = igraph::get.adjacency(graph, type = "both", attr = attr, sparse = getIgraphOpt("sparsematrices"))
  rown = colnames(adjacency)
  coln = rownames(adjacency)
  cla = class(adjacency)
  coeff = Matrix::rowSums(adjacency)
  coeff[coeff == 0] = 1
  adjacency = adjacency/coeff
  colnames(adjacency) = coln
  rownames(adjacency) = rown
  return(as(SparseM::t(adjacency), cla))
}

RWR = function(Graph, Seeds, r = .5, attr = NULL, nCores = 1, adjacency = NULL,  verbose = FALSE, stop_delta = 1e-10){
  if(!is.igraph(Graph)){
    stop("Not igraph object")
  }
  if(!is.character(Seeds)){
    stop("Seeds not characters")
  }
  Seeds = intersect(vnames(Graph), Seeds)
  if(length(Seeds) == 0){
    stop("Seeds not present in graph")
  }
  if(r < 0 | r > 1){
    stop("Restart probability not recognized")
  }
  if(is.null(adjacency)){
    get.adjacency_matrix = function(graph, attr = NULL){
      adjacency = igraph::get.adjacency(graph, type = "both", attr = attr, sparse = getIgraphOpt("sparsematrices"))
      rown = colnames(adjacency)
      coln = rownames(adjacency)
      cla = class(adjacency)
      coeff = Matrix::rowSums(adjacency)
      coeff[coeff == 0] = 1
      adjacency = adjacency/coeff
      colnames(adjacency) = coln
      rownames(adjacency) = rown
      return(as(SparseM::t(adjacency), cla))
    }
    adjacency = get.adjacency_matrix(graph = Graph, attr = attr)
  }
  prox_vector = matrix(0, nrow = ncol(adjacency), ncol = 1)
  prox_vector[which(colnames(adjacency) %in% Seeds)] = 1/length(Seeds)
  delta=1
  restart_vector = prox_vector
  while(delta >= stop_delta){
    old_prox_vector = prox_vector
    prox_vector = (1 - r) * (adjacency %*% prox_vector) + r * restart_vector
    delta = sqrt(sum((prox_vector - old_prox_vector)^2))
  }
  if(isTRUE(verbose)){message("Random Walk with restart finished")}
  Results = data.frame(Node = rownames(prox_vector), Score = prox_vector[, 1])
  Results = Results[order(Results$Score, decreasing = T), ]
  `rownames<-`(Results[-which(Results$Node %in% Seeds), ], NULL)
}

extract_by_rwr_inv = function(Graph, start_nodes, directionel = TRUE, weight = FALSE, nCores = 1, restart = .5, adj_method = get.adjacency_matrix, processed = TRUE){
  Graph_inv = igraph::reverse_edges(Graph)
  extract_by_rwr(Graph = Graph_inv, start_nodes = start_nodes, directionel = directionel, weight = weight, nCores = nCores, restart = restart, adj_method = adj_method, processed = processed)
}

extract_by_rwr_ind = function(Graph, start_nodes, directionel = TRUE, weight = FALSE, nCores = 1, restart = .5, adj_method = get.adjacency_matrix, processed = TRUE){
  Graph_ind = igraph::as.undirected(Graph, mode = "collapse")
  extract_by_rwr(Graph = Graph_ind, start_nodes = start_nodes, directionel = directionel, weight = weight, nCores = nCores, restart = restart, adj_method = adj_method, processed = processed)
}


###########################################
#  extract_by_shortest_paths
extract_by_shp = function(Graph, start_nodes, to = vnames(Graph), weight = FALSE, mode = "out"){
  if (nbrunique(start_nodes) != nbr(start_nodes)) {
    stop("start_nodes must have only unique values")
  }
  # if weight is false
  if(isFALSE(weight)) {
    dist_matrix = distances(graph = Graph, v = start_nodes, to = to, mode = mode, weights = NA)
    # when no path exist, sometimes create Inf values
    dist_matrix[!is.finite(dist_matrix)] = NA
  } else {
    # if weight is character, this assume that all weights are between 0 and 1, and that it represents a similarity, not a distance.
    # That is in the network, two nodes separated by a edge of weight 1, have a distance of 0
    # If  weight  =~ 0, then no information can pass along the edge, the distance should be infinite.
    edge_weight = get.edge.attribute(Graph, name = weight)
    log_weight = -log10(edge_weight)
    dist_matrix = 10^-distances(graph = Graph, v = intersect(start_nodes, vnames(Graph)), to = to, mode = mode, weights = log_weight)
    # when no path exist, sometimes create Inf values
    dist_matrix[dist_matrix == 0] = NA
    dist_matrix = abs(1-dist_matrix)
  }
  # data.frames are better
  dist_matrix = as.data.frame(dist_matrix)
  dist_matrix = cbind(Target = rownames(dist_matrix), dist_matrix)
  rownames(dist_matrix) = NULL
  return(dist_matrix)
}

extract_by_shp_inv = function(Graph, start_nodes, to = vnames(Graph), weight = FALSE, mode = "out"){
  Graph_inv = igraph::reverse_edges(Graph)
  return(extract_by_shp(Graph = Graph_inv, start_nodes = start_nodes, to = to, mode = mode, weight = weight))
}

extract_by_shp_ind = function(Graph, start_nodes, to = vnames(Graph), weight = FALSE, mode = "out"){
  Graph_ind = igraph::as.undirected(Graph, mode = "collapse")
  return(extract_by_shp(Graph = Graph_ind, start_nodes = start_nodes, to = to, mode = mode, weight = weight))
}

#############################################
#  extract_basic_stats
extract_topolo = function(Graph, start_nodes = vnames(Graph), weight = FALSE){
  if(isFALSE(weight)){
    weight = NA
  } else {
    weight = get.edge.attribute(Graph, weight)
  }
  results = list(
    Degree_in = degree(Graph, v = start_nodes, mode = "in"),
    Degree_out = degree(Graph, v = start_nodes, mode = "out"),
    Degree_all = degree(Graph, v = start_nodes, mode = "all"),
    Closeness_out = signif(closeness(Graph, vids = start_nodes, mode = "out"), 3),
    Closeness_all = signif(closeness(Graph, vids = start_nodes, mode = "all"), 3),
    Eccentricity_out = eccentricity(Graph, vids = start_nodes, mode = "out"),
    Eccentricity_all = eccentricity(Graph, vids = start_nodes, mode = "all"),
    Eigen_centrality = signif(eigen_centrality(Graph, scale = 1, weights = weight)$vector[start_nodes], 3),
    Betweenness = signif(betweenness(Graph, v = start_nodes, weights =  weight), 3),
    Harmonic_centrality = signif(harmonic_centrality(graph = Graph, vids = start_nodes, weights = weight, mode = "out"), 4))
  
  # similarity measures
  Similarity.invlogweighted = similarity.invlogweighted(graph = Graph, mode = "out", vids = start_nodes)
  rownames(Similarity.invlogweighted) = start_nodes
  colnames(Similarity.invlogweighted) = V(Graph)$name
  diag(Similarity.invlogweighted[start_nodes, start_nodes]) = 1
  colnames(Similarity.invlogweighted) = paste0("Similarity_", colnames(Similarity.invlogweighted))
  Similarity.invlogweighted = signif(Similarity.invlogweighted, 3)
  
  names(results) = paste0("Topology_", names(results))
  results = c(results, list(Similarity.invlogweighted))
  
  if(!is.na(weight)){
    Diversity = diversity(graph = as.undirected(Graph), weights = weight)
    add(results, Diversity = Diversity)
  }
  results = as.data.frame(do.call(cbind, results))[start_nodes, ]
  results = cbind(Target = rownames(results), results)
  rownames(results) = NULL
  return(results)
}

#####################################
#  extract_cliques
extract_clique = function(Graph, start_nodes = vnames(Graph), keep.cliques = FALSE, min = 15, max = 50){
  if(any(str_detect(vnames(Graph), "Path_"))){
    Pathway_prots = lapply(ego(Graph, mindist = 1, nodes = vnames(Graph, type = "Path_")), names)
    Pathway_prots = Pathway_prots[nbrs(Pathway_prots) >=10 & nbrs(Pathway_prots) <= 50]
    Sub_graph = induced.subgraph(Graph, str_subset(unlist(Pathway_prots), "Prot_"))
  } else {
    stop("No pathways found for the determination of cliques")
  }
  
  # Only proteins are present in the network
  Sub_graph = as.undirected(induced.subgraph(Sub_graph, vnames(Sub_graph, type = "Prot_")))
  G_clique_max = lapply(igraph::max_cliques(graph = Sub_graph, min = min, max = max), names)
  names(G_clique_max) = paste0("Clique_", 1:nbr(G_clique_max))
  return(list(Graph = Sub_graph, Cliques = G_clique_max))
}

extract_clidis = function(Graph, start_nodes = vnames(Graph), dist_matrix, min = 15, max = 50, weight = FALSE, mode = "out", Ncores = 1){
  find_big_clique = function(cliques, node) {
    max_element = NULL
    max_length = 0
    for (element in cliques) {
      if (node %in% element && length(element) > max_length) {
        max_element = element
        max_length = length(element)
      }
    }
    
    return(max_element)
  }
  
  message("Calculating cliques")
  Clique_list = extract_clique(Graph = Graph, start_nodes = start_nodes, keep.cliques = TRUE, min = min, max = max)
  Graph = Clique_list$Graph
  Cliques = Clique_list$Cliques
  
  message("Select bigest cliques")
  nodes = vnames(Graph, type = "Prot_")
  clique_to_keep = parallel::mclapply(seq_along(nodes), function(i) {print(i) ; find_big_clique(cliques = Cliques, node = nodes[i])}, mc.cores = Ncores)
  clique_to_keep = clique_to_keep[nbrs(clique_to_keep) > 0]
  clique_to_keep = unique(clique_to_keep)
  names(clique_to_keep) = paste0("Clique_", 1:nbr(clique_to_keep))
  
  message("Calculate distances to cliques")
  Cliques_Means = lapply(clique_to_keep, function(clic) rowMeans(dist_matrix[, clic]))
  Cliques_Means = do.call(cbind, Cliques_Means) %>% as.data.frame
  colnames(Cliques_Means) = paste0("Clique_", 1:ncol(Cliques_Means))
  Cliques_Means = cbind(Target = dist_matrix$Target, Cliques_Means)
  return(list(Cliques = clique_to_keep, Clique_dist_matrix = Cliques_Means))
}

extract_cludis = function(Graph, start_nodes = vnames(Graph), Clusters, dist_matrix, weight = FALSE, mode = "out"){
  # Calculate distances of each strart nodes to each clusters given by extract_clusters using the median distance of all the shortest path from the start node to each member of the cluster. 
  message("Calculate distances to clusters")
  Clusters_Mean = lapply(Clusters, function(cl) rowMeans(dist_matrix[, intersect(cl, colnames(dist_matrix))]))
  Clusters_Mean = as.data.frame(do.call(cbind, Clusters_Mean))
  Clusters_Mean = cbind(Target = dist_matrix$Target, Clusters_Mean)
  return(Clusters_Mean)
}

######################################
# extract_distance_signature
extract_siprot = function(Graph, start_nodes = vnames(Graph), Clusters, dist_matrix, Signature_prot, Signature_name, Ncores = 1){
  Signature_prot = intersect(Signature_prot, vnames(Graph))
  Clusters_sign = lapply(Clusters, function(clust) intersect(clust, Signature_prot))
  Clusters_sign = Clusters_sign[nbrs(Clusters_sign) > 1]
  
  Clusters_sign_Mean = lapply(Clusters_sign, function(cl) rowMeans(dist_matrix[, intersect(cl, colnames(dist_matrix))]))
  names(Clusters_sign_Mean) = paste0(Signature_name, "_protein_", names(Clusters_sign_Mean))
  dist_full_signature = setNames(list(rowMeans(dist_matrix[, intersect(Signature_prot, colnames(dist_matrix))])), paste0(Signature_name, "_protein_Full"))
  return(data.frame(Target = dist_matrix$Target, dist_full_signature, Clusters_sign_Mean))
}

extract_siprot2 = function(Graph, start_nodes = vnames(Graph), Clusters, dist_matrix, Signatures_prot, Ncores = 1){
  results_per_signature = list()
  for(i in seq_along(Signatures_prot)){
    Signature_prot = Signatures_prot[[i]]
    prot_signature_name = names(Signatures_prot)[i]
    df = extract_siprot(Graph = Graph, start_nodes = start_nodes, dist_matrix = dist_matrix, Clusters = Clusters, 
                        Signature_prot = Signature_prot, Signature_name = prot_signature_name, Ncores = Ncores)
    add(results_per_signature, df[, -1])
  }
  return(data.frame(Target = dist_matrix$Target, do.call(cbind, results_per_signature)))
}

extract_sigene = function(Graph, start_nodes = vnames(Graph), Clusters, dist_matrix, dist_matrix_inv, Signature_gene, Signature_name, Ncores = 1){
  # Prepare data 
  Signature_gene = intersect(Signature_gene, vnames(Graph))
  Clusters_sign = lapply(Clusters, function(clust) intersect(clust, Signature_gene))
  Clusters_sign = Clusters_sign[nbrs(Clusters_sign) > 1]
  
  # Distances from targets to gene signature
  Clusters_sign_Mean_to = lapply(Clusters_sign, function(cl) rowMeans(dist_matrix[, intersect(cl, colnames(dist_matrix))]))
  names(Clusters_sign_Mean_to) = paste0("To_", Signature_name, "_gene_", names(Clusters_sign_Mean_to))
  dist_full_signature_to = setNames(list(rowMeans(dist_matrix[, intersect(Signature_gene, colnames(dist_matrix))])), paste0("To_", Signature_name, "_protein_Full"))
  
  # Distances from gene signature to targets.
  Clusters_sign_Mean_from = lapply(Clusters_sign, function(cl) rowMeans(dist_matrix_inv[, intersect(cl, colnames(dist_matrix_inv))]))
  names(Clusters_sign_Mean_from) = paste0("From_", Signature_name, "_gene_", names(Clusters_sign_Mean_from))
  dist_full_signature_from = setNames(list(rowMeans(dist_matrix_inv[, intersect(Signature_gene, colnames(dist_matrix_inv))])), paste0("From_", Signature_name, "_protein_Full"))
  
  return(data.frame(Target = dist_matrix$Target, dist_full_signature_to, Clusters_sign_Mean_to, dist_full_signature_from, Clusters_sign_Mean_from))
}

extract_sigene2 = function(Graph, start_nodes = vnames(Graph), Clusters, dist_matrix, dist_matrix_inv, Signatures_gene, Ncores = 1){
  results_per_signature = list()
  for(i in seq_along(Signatures_gene)){
    Signature_gene = Signatures_gene[[i]]
    gene_signature_name = names(Signatures_gene)[i]
    df = extract_sigene(Graph = Graph, start_nodes = start_nodes, dist_matrix = dist_matrix, dist_matrix_inv = dist_matrix_inv, Clusters = Clusters, 
                        Signature_gene = Signature_gene, Signature_name = gene_signature_name, Ncores = Ncores)
    add(results_per_signature, df[, -1])
  }
  return(data.frame(Target = dist_matrix$Target, do.call(cbind, results_per_signature)))
}

######################################
# Create fake signature in the protein and gene layers
extract_random = function(Graph, start_nodes = vnames(Graph), dist_matrix, Clusters_prot, Signatures_prot, niter = 1e8, TT){
  available_proteins = setdiff(vnames(Graph, type = "Prot_"), TT)
  available_genes = setdiff(vnames(Graph, type = "Gene_"), paste0("Gene_", mapIds(org.Hs.eg.db, str_remove(TT, "Prot_"), column = "SYMBOL", keytype = "UNIPROT")))
  
  Signatures_prot = lapply(Signatures_prot, function(x) intersect(x, available_proteins))
  
  # Create random clusters
  Fake_cluster = lapply(ego(graph = Graph, nodes = sample(available_proteins, nbr(Clusters_prot))), function(set) str_subset(names(set), "Prot_"))
  Fake_cluster = lapply(Fake_cluster, function(clust) sample(clust, size = round(nbr(clust)*0.8)))
  Fake_cluster = Fake_cluster[nbrs(Fake_cluster) > 2]
  names(Fake_cluster) = paste0("Fake_Clust_", 1:nbr(Fake_cluster))
  
  # Create random signatures
  Fake_signature = lapply(nbrs(Signatures_prot), function(N) sample(available_proteins, size = N))
  names(Fake_signature) = paste0("Fake_signature_", 1:nbr(Fake_signature))
  # Distances to fake clusters or fake signatures
  Fake_cluster_dist = extract_cludis(Graph = Graph, start_nodes = start_nodes, Clusters = Fake_cluster, dist_matrix = dist_matrix, weight = FALSE, mode = "out")[, -1]
  Fake_signature_dist = extract_siprot2(Grap = Graph, start_nodes = start_nodes, Clusters = Clusters_prot, dist_matrix = dist_matrix, Fake_signature, Ncores = 1)[, -1]
  
  # Shuffle network while keeping edge distribution
  message("Random graph 1")
  set.seed(123)
  Random_Graph = rewire(graph = Graph, with = keeping_degseq(niter = niter))
  
  SHP_mat_fake1 = extract_by_shp(Graph = Random_Graph, start_nodes = start_nodes, to = vnames(Random_Graph))[, -1]
  colnames(SHP_mat_fake1) = paste0("Fake_interaction_SHP_", colnames(SHP_mat_fake1))
  
  RWR_mat_fake1 = extract_by_rwr(Graph = Random_Graph, start_nodes = start_nodes)[, -1]
  colnames(RWR_mat_fake1) = paste0("Fake_interaction_RWR_", colnames(RWR_mat_fake1))
  
  
  # Really shuffle network
  message("Random graph 2")
  all_node = vnames(Graph)
  Graph_df = get.data.frame(Graph)
  set.seed(123)
  Graph_df$from = sapply(seq_along(Graph_df$from), function(i) sample(all_node, 1))
  Graph_df$to   = sapply(seq_along(Graph_df$to)  , function(i) sample(all_node, 1))
  Random_Graph2 = graph_from_data_frame(unique(Graph_df), directed = T)
  
  SHP_mat_fake2 = extract_by_shp(Graph = Random_Graph2, start_nodes = start_nodes, to = vnames(Random_Graph2))[, -1]
  colnames(SHP_mat_fake2) = paste0("Fake_random_SHP_", colnames(SHP_mat_fake2))
  
  RWR_mat_fake2 = extract_by_rwr(Graph = Random_Graph2, start_nodes = start_nodes)[, -1]
  colnames(RWR_mat_fake2) = paste0("Fake_random_RWR_", colnames(RWR_mat_fake2))
  
  return(list(Results = data.frame(Target=dist_matrix$Target, 
                                   SHP_mat_fake1, RWR_mat_fake1, 
                                   SHP_mat_fake2, RWR_mat_fake2, 
                                   Fake_cluster_dist, Fake_signature_dist), 
              Fake_Network1 = Random_Graph, Fake_Network2 = Random_Graph2))
}

extract_random_top = function(Graph, start_nodes = vnames(Graph), niter = 1e8){
  # Shuffle network while keeping edge distribution
  message("Random graph 1")
  set.seed(123)
  Random_Graph = rewire(graph = Graph, with = keeping_degseq(niter = niter))
  TOP_mat_1 = extract_topolo(Graph = Random_Graph, start_nodes = start_nodes)
  Target = TOP_mat_1[[1]]
  TOP_mat_1 = TOP_mat_1[, -1]
  names(TOP_mat_1) = paste0("Fake_interaction_TOP_", names(TOP_mat_1))
  
  # Really shuffle network
  message("Random graph 2")
  all_node = vnames(Graph)
  Graph_df = get.data.frame(Graph)
  set.seed(123)
  Graph_df$from = sapply(seq_along(Graph_df$from), function(i) sample(all_node, 1))
  Graph_df$to   = sapply(seq_along(Graph_df$to)  , function(i) sample(all_node, 1))
  Random_Graph2 = graph_from_data_frame(unique(Graph_df), directed = T)
  TOP_mat_2 = extract_topolo(Graph = Random_Graph2, start_nodes = start_nodes)[, -1]
  names(TOP_mat_2) = paste0("Fake_random_TOP_", names(TOP_mat_2))
  
  return(list(Results = data.frame(Target=Target, TOP_mat_1, TOP_mat_2), Fake_Network1 = Random_Graph, Fake_Network2 = Random_Graph2))
}
