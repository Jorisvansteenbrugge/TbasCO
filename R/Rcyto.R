library('RCy3')
library(RJSONIO)

Parse_ruleframe <- function(rules.dataframe){
  .Parse_ruleframe_line <- function(line){
    rule <- unlist(strsplit(line[1], " => "))
    rule[1] <- gsub("[{]", "",  rule[1])
    rule[1] <- gsub("[}]", "",  rule[1])
    rule[2] <- gsub("[{]", "",  rule[2])
    rule[2] <- gsub("[}]", "",  rule[2])


    lhs <- as.character(unlist(strsplit(rule[1], ",")))
    rhs <- as.character(unlist(strsplit(rule[2], ",")))

    if (rule[1] == ""){
      rule[1] = "_"
    }
    if (rule[2] == ""){
      rule[2] = "_"
    }

    return(list(
      "lhs"        = lhs,
      "rhs"        = rhs,
      "support"    = as.character(line[2]),
      "confidence" = as.character(line[3]),
      "lift"       = as.character(line[4])
    ))
  }

  .Get_unique_nodes <- function(list_entry){
    return(c(list_entry$lhs, list_entry$rhs))
  }
  rule_list <- apply(rules.dataframe, 1, .Parse_ruleframe_line)
  node_list <- unique(unlist(lapply(rule_list, .Get_unique_nodes)))
  names(node_list) <- NULL
  names(rule_list$support) <- NULL
  return(list("node_list"        = node_list,
              "interaction_list" = rule_list))
}

#' Draw a network of association rules in Cytoscape
#' @name Network_Association_Rules
#' @description Association rules between trait-attributes can be visualized in
#' a directed graph network.
#' @param rules.dataframe A dataframe containing association rules,
#' as calculated by \code{\link{Association_Rules}}.
#' @param annotation.db List containing a dictionary like structure with traits as
#' names and annotation as values.
#' @param N The number of association rules to visualise
#' @export
#' @author JJM van Steenbrugge
Network_Association_Rules <- function(rules.dataframe, annotation.db, N = 200){
  if(!"RCy3" %in% installed.packages()){
    source("https://bioconductor.org/biocLite.R")
    biocLite("RCy3")
  }
  library(RCy3)
  g <- new('graphNEL', edgemode='directed')
  g <- initNodeAttribute(graph = g,
                         attribute.name  = 'attribute',
                         attribute.type  = 'char',
                         default.value   = 'undefined')
  g <- initNodeAttribute(graph = g,
                         attribute.name = 'pathway_module',
                         attribute.type = 'char',
                         default.value  =  'undefined')
  g <- initNodeAttribute(graph = g,
                         attribute.name = 'label',
                         attribute.type = 'char',
                         default.value = 'default_vals')

  g <- initEdgeAttribute(graph = g,
                         attribute.name = 'edgeType',
                         attribute.type = 'char',
                         default.value  = 'Arrow')
  g <- initEdgeAttribute(graph = g,
                         attribute.name = 'support',
                         attribute.type = 'char',
                         default.value  = 'undefined')
  g <- initEdgeAttribute(graph = g,
                         attribute.name = 'confidence',
                         attribute.type = 'char',
                         default.value  = 'undefined')
  g <- initEdgeAttribute(graph = g,
                         attribute.name = 'lift',
                         attribute.type = 'char',
                         default.value  = 'undefined')


  rules_list <- Parse_ruleframe(rules.dataframe[1:N,])
 # pathway_modules <- getPathwayName(annotation.db)


  ## Adding Nodes
  added_nodes <- c()
  for(nodes in rules_list$node_list){
    if(! is.na(nodes) ){

      nodes <- unlist(strsplit(nodes, ","))
      for(node in nodes){

        if(! node %in% added_nodes){
          g <- graph::addNode(node, g)

          added_nodes <- c(node, added_nodes)

          trait <- unlist(strsplit(node, "[.]"))[1]
      #    nodeData(g, node, 'pathway_module') <- pathway_modules[[trait]]
          nodeData(g, node, 'attribute') <- 'module'
        }
      }
    }
  }

  print("NODES ADDED")
  added_interactions <- c()
  inter_nodes <- c()
  ## Adding Edges
  for(interaction in rules_list$interaction_list){
    inter <- paste(paste(interaction$lhs, collapse = "&"),
                   paste(interaction$rhs, collapse = "&"),
                   sep = "+")
    if(! inter == "NA+NA"){
      # Make sure we don't overwrite duplicates just to be sure
      if(! inter %in% added_interactions){

        if (length(interaction$lhs) == 1 && length(interaction$rhs) == 1) {
          g <- graph::addEdge(interaction$lhs, interaction$rhs, g)
          edgeData(g, interaction$lhs, interaction$rhs, 'support') <- as.character(interaction$support)
          edgeData(g, interaction$lhs, interaction$rhs, 'confidence') <- as.character(interaction$confidence)
          edgeData(g, interaction$lhs, interaction$rhs, 'lift') <- as.character(interaction$lift)
        } else{
          # Add interaction node
          g <- graph::addNode(inter, g)
          inter_nodes <- c(inter, inter_nodes)
          nodeData(g, inter, 'attribute') <- 'interaction'
          nodeData(g, inter, 'label') <-  ""

          # Attach LHS to the interaction node
          lhs <- unlist(strsplit(interaction$lhs, ","))
          for(node in lhs){
            node <- gsub("[{]", "", node)
            node <- gsub("[}]", "", node)
            g <- graph::addEdge(node, inter, g)
            edgeData(g, node, inter, 'support')    <- as.character(interaction$support)
            edgeData(g, node, inter, 'confidence') <- as.character(interaction$confidence)
            edgeData(g, node, inter, 'lift')       <- as.character(interaction$lift)
          }
        # Attach RHS from the interaction node to the others
          rhs <- unlist(strsplit(interaction$rhs, ","))
          for(node in rhs){
            node <- gsub("[{]", "", node)
            node <- gsub("[}]", "", node)
            g <- graph::addEdge(inter, node, g)
            edgeData(g, inter, node, 'support')    <- as.character(interaction$support)
            edgeData(g, inter, node, 'confidence') <- as.character(interaction$confidence)
            edgeData(g, inter, node, 'lift')       <- as.character(interaction$lift)
          }
        }
        added_interactions <- c(inter, added_interactions)
      }
    }
  }

  cw <- CytoscapeWindow('vignette',
                        graph = g,
                        overwrite = TRUE)
  setEdgeTargetArrowRule(cw, 'edgeType', 'Associates_with', 'Associates_with', default = 'Arrow')


  top_10_modules <- c("Biosynthesis_of_secondary metabolites", "Cell_signaling",
                      "Two-component_regulatory_system","Cofactor_and_vitamin_biosynthesis",
                      "Central_carbohydrate_metabolism", "ATP_synthesis", "Ribosome",
                      "Carbon_fixation", "Other_carbohydrate_metabolism", "Drug_resistance",
                      "Drug efflux transporter/pump", "Aromatic amino acid metabolism",
                      "RNA processing", "Bacterial secretion system", "Methane metabolism",
                      "Spliceosome", "Aromatics degradation", "Saccharide, polyol, and lipid transport system",
                      "Other amino acid metabolism", "Serine and threonine metabolism")

  colors <- c("#070D00", "#00FF00", "#341D83", "#483D84", "#5E568D", "#746E9A",
              "#8A86A7", "#A09EB5", "#B6B5C2", "#CBCBCF", "#CBCBCF", "#B4B6C2",
              "#9D9FB5", "#8488A8", "#6C719A", "#525A8E", "#374284", "#062983",
              "#0000FF", "#0F0A00")

  setNodeColorRule (
    cw,
    node.attribute.name = 'pathway_module',
    control.points = top_10_modules,
    colors = colors,
    mode = 'lookup',
    default.color = '#AA0000'
  )

  setNodeShapeRule(
    cw,
    node.attribute.name = 'attribute', c('module', 'interaction'),
    c('RECTANGLE', "ellipse")
  )


  layoutNetwork(cw, layout.name = 'grid')
  displayGraph(cw)

   setNodeFontSizeDirect(cw,
                         inter_nodes,
                         1)
  setNodeSizeDirect(cw,
                    inter_nodes,
                    30)

}





# Function to convert the traits table into a list
getPathwayName <- function(annotation.db) {
  traits <- list()
  for(i in 1:nrow(annotation.db)){
    traits[[ annotation.db[i, 2] ]] <- annotation.db[i, 1]
  }
  return(traits)
}
