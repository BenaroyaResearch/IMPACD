#' Transforms values to same space as flowJo workspaces
#'
#' This function takes values and markers and transforms the values into the 
#' same as in the flowJo workspace
#'
#' @param values value or values to transform
#' @param marker string that identifies the transformation for that marker
#' @param transList_wsp a transformer list object
#' @return transformed value(s)
#' @examples
#' # Example on how to obtain the transformation list
#' ws = open_flowjo_xml("flowJo.wsp"))
#' gs_flowJo = flowjo_to_gatingset(ws, name = "batch#1")
#' trans = gh_get_transformations(gs_flowJo[[1]])
#' transList_wsp <- transformList(from = names(trans), trans)
#' 
#' Apply to a value
#' transformed_value = transform_to_wsp_space(140, "Sm152Di", transList_wsp)
transform_to_wsp_space = function(values, marker, transList_wsp) {
  return(transList_wsp@transforms[[marker]]@f(values))
}

#' Flip the gate
#'
#' Flips CD3+ to CD3- and vice versa. If string does not contain `+` or `-`, 
#' same string is returned.
#'
#' @param gateName string with gate name
#' @return string with flipped gate
#' @examples
#' gateName("CD3+")
flipGateName = function(gateName) {
  if (str_detect(gateName, "\\+")) {
    return(str_replace(gateName, pattern = "\\+", replacement = "\\-"))
  } else if (str_detect(gateName, "\\-")) {
    return(str_replace(gateName, pattern = "\\-", replacement = "\\+"))
  } else {
    print("No '+/-' to switch in gate name.")
    return(gateName)
  }
}

#' Convert characters to numerics
#'
#' A function that can take a character and convert it to a numeric value.
#' Needed for processing the gating information in some cases, most usually 
#' logical gates. If the symbol provided is `-`, Inf, or -Inf is returned.
#'
#' @param coordinates array with possible characters or `-` symbols
#' @return array after conversion to numeric values
#' @examples
#' x = convert_to_numeric(c("-", "3"))
#' print(x)
#' [1] -Inf 3
#' 
#' #' x = convert_to_numeric(c("3", "-"))
#' print(x)
#' [1] 3 Inf
convert_to_numeric = function(coordinates) {
  converted = c()
  
  if (coordinates[1] == "-") {
    converted[1] = -Inf
  }
  else {
    converted[1] = as.numeric(coordinates[1])
  }
  
  if (coordinates[2] == "-") {
    converted[2] = Inf
  }
  else {
    converted[2] = as.numeric(coordinates[2])
  }
  
  return(converted)
}

#' Preprocess the gate definitions
#'
#' A function that takes a data frame with gate definition (presumably an excel
#' exported from flowJo), and returns a list. First element in the list is clean
#' data frame that contains all information about the gates and second element 
#' is a list with all gate objects.
#'
#' @param gate_definitions data frame with gate definitions, must contain the 
#' following columns: `Annotation` (name of gate, eg. CD3), `Marker Type`, 
#' `Name`, `Gate Position`
#' @param transList transformer list object, extracted from flowJo workspace
#' @return list with first element being a clean
#' data frame that contains all information about the gates and second element 
#' a list with all gate objects.
#' @examples
#' # Extract the transformation list from flowJo workspace
#' ws = open_flowjo_xml("flowJo.wsp"))
#' gs_flowJo = flowjo_to_gatingset(ws, name = "batch#1")
#' trans = gh_get_transformations(gs_flowJo[[1]])
#' transList_wsp <- transformList(from = names(trans), trans)
#' 
#' # Gate defintions (will usually be loaded from an excel), should look similar to:
#' gate_definitions = data.frame(Annotation = "CD49b", Name = "Dy161Di", "Marker Type" = "1D", "Gate Position" = 76.4)
#' 
#' processed_gates = processGates(gate_definitions_, transList_wsp) 
#` gate_metadata = processed_gates[[1]]
#` gates = processed_gates[[2]]
#' @export
processGates = function(gate_definitions, transList) {
  
  gate_metadata = c()
  gates = list()
  
  for (iRow in 1:nrow(gate_definitions)) {
    gate_definition = gate_definitions[iRow, ]
    # print(gate_definition$Annotation)
    
    metadata = data.frame(name = c(paste0(gate_definition$Annotation, "+"), paste0(gate_definition$Annotation, "-")))
    
    gate_type = gate_definition$`Marker Type`
    metadata$gate_type = gate_type
    
    # Get the markers/channels
    if (gate_type == "1D") {
      marker1 = gate_definition$Name
      marker2 = ""
      
    } else {
      marker1 = str_split(gate_definition$Name, " vs. ")[[1]][1]
      marker2 = str_split(gate_definition$Name, " vs. ")[[1]][2]
    }
    metadata$marker1 = marker1
    metadata$marker2 = marker2
    
    # Extract original gate definition for the positive gate
    original_gate_position = gate_definition$`Gate Position`
    metadata$original_gate_position = original_gate_position
    
    # Handle logic gates (not, and, or)
    if (str_detect(original_gate_position, "\\(NOT\\)")) {
      logical_gate = "NOT"
      coordinates = str_remove(original_gate_position, "\\(NOT\\) ")
      
    } else if (str_detect(original_gate_position, "\\(OR\\)")) { # Still need to work out the logic for this gate
      logical_gate = "OR"
      subgates_ = str_split(str_remove(original_gate_position, "\\(OR\\) "), "/", simplify = TRUE)
      nsubgates = ncol(subgates_)
      subgates_ = subgates_[, 2:nsubgates]
      
      coordinates = list()
      
      # Extract the coordinates from the string in the excel
      for (isubgate in 1: (nsubgates - 1)) {
        coordinates_ = str_extract_all(subgates_[isubgate],"\\(?[-0-9.]+\\)?")[[1]]
        coordinates_ = coordinates_[2:length(coordinates_)]
        # Turn into a matrix,
        coordinates_ = matrix(c(coordinates_[seq(1, length(coordinates_), 2)], 
                                coordinates_[seq(2, length(coordinates_), 2)]), 
                              ncol = 2, 
                              dimnames = list(c("min", "max"), c(marker1, marker2)) )
        
        coordinates[[isubgate]] = coordinates_
      }
      
    } else if (str_detect(original_gate_position, "\\(AND\\)")) { # Still need to work out the logic for this gate
      logical_gate = "AND"
      coordinates = str_remove(original_gate_position, "\\(AND\\) ")
      
    } else {
      logical_gate = FALSE
      coordinates = original_gate_position
    }
    metadata$logical_gate = logical_gate
    
    if (logical_gate %in% c("AND")) {
      print(paste0("skipping ", gate_definition$Annotation, " since AND logic is not yet supported"))
      next
    }
    
    # Process gate coordinates and transform them into the same space as the fcs files
    if (logical_gate == "OR") {
      transformed_coordinates = list()
      for (isubgate in 1: (nsubgates - 1)) {
        coordinates_ = coordinates[[isubgate]]
        transformed_coordinates_ = matrix(nrow = 2, ncol = 2, dimnames = list(c("min", "max"), c(marker1, marker2)) )
        
        transformed_coordinates_[, 1] = transform_to_wsp_space(convert_to_numeric(coordinates_[, 1]), marker1, transList_wsp)
        transformed_coordinates_[, 2] = transform_to_wsp_space(convert_to_numeric(coordinates_[, 2]), marker2, transList_wsp)
        
        transformed_coordinates[[isubgate]] = transformed_coordinates_
        
      }
      
      coordinates = transformed_coordinates
      
    } else {
      if (gate_type == "1D") {
        coordinates = transform_to_wsp_space(as.numeric(coordinates), marker1[1], transList_wsp)
        
        coordinates = list(c(coordinates, Inf))
        names(coordinates) = marker1[1]
        
      } else {
        # Extract the coordinates from the string in the excel
        coordinates = as.numeric(str_extract_all(coordinates,"\\(?[-0-9.]+\\)?")[[1]])
        
        # Turn into a matrix,
        coordinates = matrix(c(coordinates[seq(1, length(coordinates), 2)], 
                               coordinates[seq(2, length(coordinates), 2)]), 
                             ncol = 2)
        
        # convert
        coordinates[ , 1] = transform_to_wsp_space(coordinates[ , 1], marker1, transList_wsp)
        coordinates[ , 2] = transform_to_wsp_space(coordinates[ , 2], marker2, transList_wsp)
        
        colnames(coordinates) = c(marker1, marker2)
        
      }
    }
    
    # Write transformed values to dataframe
    metadata$transformed_gate_position = list(coordinates)
    
    
    if (logical_gate == FALSE) {
      # Constructs gates with flowCore
      # Builds the "positive" gate, assuming this is not a logical gate
      pos_gate = if (gate_type == "1D") rectangleGate(filterId = metadata$name[1], coordinates) else pos_gate = polygonGate(filterId = metadata$name[1], coordinates)
      
      exprs = paste0("!", metadata$name[1])
      call = substitute(booleanFilter(v, filterId = metadata$name[2]), list(v = as.symbol(exprs)))
      neg_gate = eval(call)
    } else if (logical_gate == "NOT") {
      neg_gate = if (gate_type == "1D") rectangleGate(filterId = metadata$name[2], coordinates) else pos_gate = polygonGate(filterId = metadata$name[2], coordinates)
      
      exprs = paste0("!", metadata$name[2])
      call = substitute(booleanFilter(v, filterId = metadata$name[1]), list(v = as.symbol(exprs)))
      pos_gate = eval(call)
      
    } else if (logical_gate == "OR") {
      gatename_ = str_remove(metadata$name[1], "\\+|-")
      
      subgates = c()
      for (isubgate in 1: (nsubgates - 1)) {
        coordinates_ = coordinates[[isubgate]] 
        subgate = paste(gatename_, isubgate, sep = "_")
        subgates = append(subgates, subgate)
        gates[[subgate]] = rectangleGate(filterId = subgate, .gate = coordinates_)
      }
      
      exprs_pos = paste0(subgates, collapse = "|")
      call = substitute(booleanFilter(v, filterId = metadata$name[1]), list(v = as.symbol(exprs_pos)))
      pos_gate = eval(call)
      
      exprs_neg = paste0("!", subgates, collapse = "&")
      call = substitute(booleanFilter(v, filterId = metadata$name[2]), list(v = as.symbol(exprs_neg)))
      neg_gate = eval(call)
      
    } else {
      stop("AND not yet supported.")
    }
    
    gates[[pos_gate@filterId]] = pos_gate
    gates[[neg_gate@filterId]] = neg_gate
    
    gate_metadata = rbind(gate_metadata, metadata)
    
  }
  
  return(list(gate_metadata, gates))
}