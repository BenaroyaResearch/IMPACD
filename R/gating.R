#' Apply a gate to a gating set
#'
#' This function applies a gate/node to a gating set. It can also handle logical
#' gates. Will call `add_new_node` to actually add the node after figuring out 
#' what gates need to be added.
#'
#' @param gs gating set
#' @param gateName name of the gate to be applied, comes from `processGates()`
#' @param gates list with gate objects, `gateName` must be an entry for it
#' @param gate_metadata data frame with information on gates, comes from 
#' `processGates()`
#' @param node parent node to which the gate is applied to
#' @return A gating set (`gs`) with additional node
#' @examples
#' # Apply the CD49b+ gate to the "Teff+/CD45RO+/TEM+" population
#' apply_gate(gs, "CD49b+", gates, gate_metadata, "Teff+/CD45RO+/TEM+")
apply_gate = function(gs, gateName, gates, gate_metadata, node) {
  
  if (gate_metadata$logical_gate[gate_metadata$name == gateName] == "NOT") {
    # if it is a NOT gate, add the gate that the boolean not gate is based upon
    complimentaryGateName = flipGateName(gateName)
    
    gs = add_new_node(gs, complimentaryGateName, gates, node)
    
  } else if (gate_metadata$logical_gate[gate_metadata$name == gateName] == "OR") {
    
    # remove +/- and then add "_", to identify related subgates
    subgates = names(gates)[str_detect(names(gates), paste0(str_remove(gateName, c("\\+|-")), "_"))]
    
    for (subgate_ in subgates) {
      gs = add_new_node(gs, subgate_, gates, node)
    }  
    
  } else if (str_detect(gateName, "\\-")) {
    posGateName = flipGateName(gateName)
    gs = add_new_node(gs, posGateName, gates, node)
  }
  
  gs = add_new_node(gs, gateName, gates, node)
  
  return(gs)
}

#' Adds node to a gating set
#'
#' This function applies a gate/node to a gating set. Node gets only added if it
#' does not already exist. 
#'
#' @param gs gating set
#' @param gateName name of the gate to be applied, comes from `processGates()`
#' @param gates list with gate objects, `gateName` must be an entry for it
#' @param parent_node parent node to which the gate is applied to
#' @return A gating set (`gs`) with additional node
#' @examples
#' # Apply the CD49b+ gate to the "Teff+/CD45RO+/TEM+" population
#' add_new_node(gs, "CD49b+", gates, "Teff+/CD45RO+/TEM+")
add_new_node = function(gs, gateName, gates, parent_node) {
  existingNodes = gs_get_pop_paths(gs, path = "full")
  nodeToAdd = paste(parent_node, gateName, sep = "/")
  
  if (str_remove(nodeToAdd, "root") %in% existingNodes) {return (gs)}
  
  gs_pop_add(gs, gates[[gateName]], parent = parent_node)
  
  return(gs)
  
}

#' Applies root level gating
#'
#' This function applies the root level gating to a gating set. 
#'
#' @param gs gating set
#' @param gating_strategy a data frame that contains columns with gate `name`, 
#' `type` (root or leaf) and `root_oder`, see examples
#' @param gate_metadata data frame with information on gates, comes from 
#' `processGates()`
#' @param gates list with gate objects, `gateName` must be an entry for it
#' @param starting_node optional, the starting node for gating, usually root or "", `default = ""`
#' @return A list, with first element being the gating set (`gs`) and second 
#' element being the root node path
#' @examples
#' # Generate the gating strategy data frame
#' root_path = c("Teff+", "CD45RO+", "TEM+")
#' leafs = c("CD49b+", "TIGIT+")
#' 
#' gating_strategy = data.frame(name = c(root_path, leafs))
#' gating_strategy = mutate(gating_strategy, type = case_when(name %in% root_path ~ 'root',
#                                                            name %in% leafs ~ 'leaf'))
#' gating_strategy$root_order = ""
#' gating_strategy$root_order[1:length(root_path)] = 1:length(root_path)
#' 
#' # Apply root level gating
#' gating_results = rootLevelGating(gs, gating_strategy, gate_metadata, gates)
#' 
#' gs = gating_results[[1]]
#' root_node = gating_results[[2]]
#' print(root_node)
#' [1] "root/Teff+/CD45RO+/TEM+"
#' @export
rootLevelGating = function(gs, gating_strategy, gate_metadata, gates, starting_node = "") {
  
  node = starting_node
  for (step in 1:max(gating_strategy$root_order)) {
    gateName = dplyr::filter(gating_strategy, root_order == step)$name
    
    if (node == "") { # must come after next_path above
      node = "root"
    }
    next_node = paste(node, gateName, sep = "/")
    # print(next_node)
    
    gs = apply_gate(gs, gateName, gates, gate_metadata, node)
    
    node = next_node
    
  }
  
  # recompute(gs)
  
  return(list(gs, node))
}

#' Applies level 1 gating
#'
#' This function applies the level 1 gating to a gating set after root pathing.
#' 
#' @param gs gating set
#' @param leafs an array with all the leafs to be gated
#' @param gate_metadata data frame with information on gates, comes from 
#' `processGates()`
#' @param gates list with gate objects, `gateName` must be an entry for it
#' @param node the starting node for gating, most like the root node returned 
#' from `rootLevelGating`
#' @return A list, with first element being the gating set (`gs`) and second 
#' element being a data frame with all gated paths
#' @examples
#' level1_results = level1Gating(gating_sets_by_batch[[batch_]], level1_nodes, gate_metadata_by_batch[[batch_]], gates_by_batch[[batch_]], root_node) 
#'
#' gs = level1_results[[1]]
#' level1_paths = level1_results[[2]]
#' @export
level1Gating = function(gs, leafs, gate_metadata, gates, node) {
  message( paste("Stage 1 gating for root path '", node,"'") )
  message( paste("Gating", length(leafs), "nodes") )
  
  for (leaf_name in leafs) {
    gs = apply_gate(gs, leaf_name, gates, gate_metadata, node) 
  }
  
  recompute(gs)
  
  paths = data.frame(node = leafs, parent = node, population = paste0(node, "/", leafs))
  paths$gs_path = str_remove(paths$population, "root")
  
  return(list(gs, paths))
  
}

#' Applies level 2 gating
#'
#' This function applies the level 2 gating to a gating set after root pathing.
#' Will take care of the permutations and eliminate duplicate and super 
#' low-yield paths
#' 
#' @param gs gating set
#' @param leafs an array with all the leafs to be gated
#' @param gate_metadata data frame with information on gates, comes from 
#' `processGates()`
#' @param gates list with gate objects, `gateName` must be an entry for it
#' @param starting_node the starting node for gating, most like the root node returned 
#' from `rootLevelGating`
#' @param max_depth maximum depth for gating, `default = 3`. Warning, increasing
#' this might take very, very long.
#' @param level1_nodes can be provided to restrict the space for the first level
#' gating, else all `leafs` will be gated, `default = NULL`
#' @param cellcountCutoff minimum number of cells for a population to continue going down the tree
#' @return A data frame with all gated paths.
#' @examples
#' leafs = c("CD49b", "TIGIT", "CD161")
#' 
#' root_node = "root/Teff+/CD45RO+/TEM+"
#' 
#' # Gate up to three levels deep, doing all permutations at every level
#' paths = level2Gating(gs, leafs, gate_metadata, gates, starting_node = root_node, max_depth = 3)
#' 
#' # Gate up to fours levels deep, but use only "TIGIT+" and "TIGIT-" in first gating depth
#' paths = level2Gating(gs, leafs, gate_metadata, gates, starting_node = root_node, max_depth = 4, level1_nodes = c("TIGIT+", "TIGIT-"))
#' @export
level2Gating = function(gs, leafs, gate_metadata, gates, starting_node = "root", max_depth = 3, level1_nodes = NULL, cellcountCutoff = 30) {
  
  if (starting_node == "") {starting_node = "root"}
  if (starting_node != "root" && str_detect(starting_node, "root")) {
    starting_node_ = str_remove(starting_node, "root")
  } else {
    starting_node_ = starting_node
  }
  
  # check if starting_node exists in gs
  if ( !(starting_node_ %in% gs_get_pop_paths(gs, path = "full")) ) {stop("Starting node not detected in gs, you need to use root_level_gating first if you don't start at 'root'.")}
  
  if (is.null(max_depth)) {max_depth = length(unique(str_remove(leafs, "\\+|\\-"))) } # check for number unique nodes
  
  keep_going = TRUE
  
  depth_ = 1
  paths = c()
  
  for (depth_ in 1: max_depth) {
    
    message( paste("Gating depth:", depth_) )
    
    if (depth_ == 1) {
      parent_nodes = starting_node
      
    } else {
      parent_nodes = dplyr::filter(paths, depth == depth_-1)
      parent_nodes = filter_duplicate_paths(parent_nodes)
      if (nrow(parent_nodes) == 0) {message("Exhausted leaf paths"); break}
    }
    
    if ((depth_ == 1) & !(is.null(level1_nodes))) { 
      paths_ = new_paths(level1_nodes, parent_nodes) 
    } else {
      paths_ = new_paths(leafs, parent_nodes)
    }
    
    paths_$depth = depth_
    
    message( paste("Gating", nrow(paths_), "nodes") )
    
    paths = rbind(paths, paths_)
    
    for (iRow in 1:nrow(paths_)) { # turn into apply?
      gs = apply_gate(gs, paths_$node[iRow], gates, gate_metadata, paths_$parent[iRow])
    }
    
    recompute(gs)
    
    message("Pruning low yield paths")
    counts_ = gs_pop_get_count_fast(gs, subpopulations = paths_$gs_path)
    
    paths_before_pruning = nrow(paths)
    
    paths = filter_low_count_paths(counts_, paths, cellcountCutoff = cellcountCutoff)
    rownames(paths) = NULL 
    message( paste("Eliminated", paths_before_pruning-nrow(paths), "paths with <", cellcountCutoff, "cells" )) 
    message( paste(nrow(paths), "total populations \n" )) 
  }
  
  message("Recalculating population counts")
  recompute(gs)
  
  return(list(gs, paths))
}

filter_low_count_paths = function(counts, gating_paths, cellcountCutoff = 30) {
  summary = counts %>%
    group_by(Population) %>%
    summarise(nCells = sum(Count))
  
  summary$pathsToClear = summary$nCells < cellcountCutoff
  # skip, if there aren't any paths that need to be filtered out
  if (!any(summary$pathsToClear)) {return(gating_paths)}
  
  summary$Population_regexp = str_replace_all(summary$Population, "\\+", "\\\\\\+")
  summary$Population_regexp = str_replace_all(summary$Population_regexp, "\\-", "\\\\\\-")
  
  pathsToRemove = paste0(summary$Population_regexp[summary$pathsToClear], collapse = "|") 
  
  gating_paths = gating_paths[!str_detect(gating_paths$gs_path, pathsToRemove), ] # remove paths that we do not want
  
  return(gating_paths)
}



filter_double_paths = function(df) {
  nodes = unique(str_remove(df$node, "\\+|\\-"))
  
  df$filter = FALSE
  
  # step through all nodes and make sure that each node is only represented 
  # once in each path, e.g. TIGIT+/LAG3+/LAG3- needs to be removed since 
  # LAG3 is represented as pos and negative gate
  # Also remove TIGIT+/LAG3+/LAG3+
  for (i in 1:length(nodes)) {
    df$filter = ifelse(str_detect(df$gs_path, paste0(nodes[i], "\\+")) & str_detect(df$gs_path, paste0(nodes[i], "\\-")), TRUE, df$filter)
    df$filter = ifelse( ( str_count(df$gs_path, paste0(nodes[i], "\\+")) > 1 ) | ( str_count(df$gs_path, paste0(nodes[i], "\\-")) > 1), TRUE, df$filter)
  }
  # remove all rows where filter is TRUE
  df = df[!df$filter,]
  
  # remove the filter column
  df = df %>% select(-filter)
  rownames(df) = NULL 
  
  return(df)
}


filter_duplicate_paths = function(df) {
  nodes = unique(str_remove(df$node, "\\+|\\-"))
  
  nodes = c(paste0(nodes, "\\+"), paste0(nodes, "\\-"))
  
  masks = lapply(df$population, function(x) { ifelse(str_detect(x, nodes), "1", "0") })
  
  for (i in 1:length(masks)) {
    df$mask[i] = unlist(paste0(masks[[i]], collapse = ""))
  }
  
  df = df[!duplicated(df$mask),]
  df = df %>% select(-mask)
  
  rownames(df) = NULL
  return(df)
}

new_paths = function(leafs, parent_nodes) {
  paths = c()
  
  if (is.character(parent_nodes)) {
    paths = data.frame(node = leafs, parent = parent_nodes)
    paths$population = paste(paths$parent, paths$node, sep = "/")
    paths$gs_path = str_remove(paths$population, "root")
    
  } else {
    for (iParent in 1:nrow(parent_nodes)) {
      paths_ = data.frame(node = leafs, parent = parent_nodes$population[iParent])
      paths_$population = paste(paths_$parent, paths_$node, sep = "/")
      paths_$gs_path = str_remove(paths_$population, "root")
      
      paths = rbind(paths, paths_)
    }
  }
  
  paths = filter_double_paths(paths)
  
  return(paths)
}