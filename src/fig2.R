library(ggplot2)
library(magrittr)
library(tidyr)
library(WikidataQueryServiceR)
library(readr)

query_for_relations_used <- read_file("src/sparql/relations_used.rq")

df = query_wikidata(query_for_relations_used)
df$index <- "cell_lines"

df %>% 
  tidyr::pivot_longer(!index,names_to="relation", values_to = "count") ->
  long_df


clean_names <- function(a){
  clean_names = list(
    "cell_lines" = "instances of 'cell line'",
    "taxons"= "taxons",
    "diseases"= "diseases",
    "CLO_ids"= "Cell Line Ontology identifier",
    "cellosaurus_ids"= "Cellosaurus identifier",
    "sources"= "describing articles",
    "MESH_ids"= "MeSH identifiers",
    "wikipedia_articles"= "English Wikipedia articles",
    "exact_match_links"= "'exact match' crossreferences",
    "cell_lines_with_daughters"= "cell lines with daughters"	
  )
  
  return(clean_names[[a]])
  
}

long_df$clean_names <- as.character(lapply(long_df$relation, clean_names))


pl_relations_used <- ggplot(long_df, aes(x = reorder(clean_names, -count), y = count)) + 
  geom_point() + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 70, vjust =, hjust=1))+
  xlab("Relation on Wikidata") + 
  ylab("# entries on Wikidata") + scale_y_log10() 

ggsave("images/cellosaurus_count_per_category_2021_12_07.pdf", width = 7, height = 4)
plot(pl_relations_used)
dev.off()

ggsave("images/cellosaurus_count_per_category_2021_12_07.png", width = 7, height = 4)
plot(pl_relations_used)
dev.off()


##### Plot the MCF-7 cell line daughter lines 

library(network)
library(sna)
library(intergraph)
library(ggplot2)
library(igraph)
library(GGally)
library(ggrepel)

query_for_mcf7_daughters <- read_file("src/sparql/mcf7_daughters.rq")
df = query_wikidata(query_for_mcf7_daughters)
df = as.matrix(df)
igraph_object <- graph_from_edgelist(df, directed = TRUE)

nodes_to_label = 10
color = "purple"
name = "MCF-7 daughter lines"

degrees <- igraph::degree(igraph_object, normalized = FALSE)
igraph_object <-  igraph::set_vertex_attr(igraph_object, "degree", value = degrees)
net_obj <- intergraph::asNetwork(igraph_object)
graph_matrix <-  network::as.matrix.network.adjacency(net_obj)

# get coordinates from Fruchterman and Reingold's force-directed plafcent algorithm.
plotcord <- data.frame(sna::gplot.layout.fruchtermanreingold(graph_matrix, NULL))

colnames(plotcord) <- c("X1", "X2")

edglist <- network::as.matrix.network.edgelist(net_obj)

edge_coordinates <- data.frame(plotcord[edglist[, 1], ], plotcord[edglist[, 2], ])
plotcord$vertex.names <-  as.factor(network::get.vertex.attribute(net_obj, "vertex.names"))
plotcord$Degree <-   network::get.vertex.attribute(net_obj, "degree")

plotcord[, "shouldLabel"] <- FALSE

max_n <- min(nodes_to_label, length(degrees))

labeled_node_names <- names(sort(degrees, decreasing = TRUE))[seq_len(max_n)]
labeled_node_cooordinates <- plotcord[, "vertex.names"] %in% labeled_node_names

colnames(edge_coordinates) <- c("X1", "Y1", "X2", "Y2")

plotcord[which(plotcord[, "vertex.names"] %in% labeled_node_names), "shouldLabel"] <-
  TRUE

plotcord$Degree_cut <-
  cut(plotcord$Degree,
      breaks = 3,
      labels = FALSE
  )

plotcord$in_mod <- TRUE

# Adapted from https://github.com/csbl-usp/fcoex/blob/master/R/visualization.R
pl_mcf7 <- ggplot(plotcord) +
  geom_segment(
    data = edge_coordinates,
    aes_(
      x = ~X1,
      y = ~Y1,
      xend = ~X2,
      yend = ~Y2
    ),
    size = 0.5,
    alpha = 0.9,
    colour = "#9e9e9e"
  ) +
  geom_point(aes_(
    x = ~X1,
    y = ~X2,
    size = ~Degree,
    alpha = ~Degree
  ), color = color) +
  geom_label_repel(
    aes_(
      x = ~X1,
      y = ~X2,
      label = ~vertex.names
    ),
    box.padding = unit(1, "lines"),
    data = function(x) {
      x[x$shouldLabel, ]
    }
  ) +
  scale_colour_manual(values = c("#005E87")) +
  labs(title = name) +
  ggplot2::theme_bw(base_size = 12, base_family = "") +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    ),
    panel.border = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank()
  )
ggsave("images/mcf7_daughter_lines_2021_12_07.png", width = 5, height = 5)
plot(pl_mcf7)
dev.off()

ggsave("images/mcf7_daughter_lines_2021_12_07.pdf", width = 5, height = 5)
plot(pl_mcf7)
dev.off()

# Organize plots into a nice figure
ggsave("images/relations_and_mcf7_combined_2021_12_07.pdf", width = 10, height = 4)
pl_relations_used + pl_mcf7 + plot_layout(widths = c(5, 3))+ plot_annotation(tag_levels = 'A')
dev.off()

ggsave("images/relations_and_mcf7_combined_2021_12_07.png", width = 10, height = 4)
pl_relations_used + pl_mcf7 + plot_layout(widths = c(5, 3))+ plot_annotation(tag_levels = 'A')
dev.off()
