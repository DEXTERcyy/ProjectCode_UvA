library(dplyr)
library(ggplot2)
library(ggraph)
library(igraph)
library(RColorBrewer)
cir_plot <- function(network_table, taxa_table, labels, plot_path)
  {
    otu_tax_df <- taxa_table[,1:5] %>% # first 5 taxas
      rownames_to_column("OTU") %>%
        dplyr::select(-OTU, everything(), OTU)

    pairs <- list(
      c("Kingdom", "Phylum"),
      c("Phylum", "Class"),
      c("Class", "Order"),
      c("Order", "Family"),
      c("Family", "OTU"))
    # Function to create edges for a single pair
    create_edges <- function(pair, data)
      {
        from <- data[[pair[1]]]
        to <- data[[pair[2]]]
        data.frame(from = from, to = to)
      }

    # Apply the function to all pairs and combine the results
    edges <- unique(do.call(rbind, lapply(pairs, create_edges, data = otu_tax_df[otu_tax_df$OTU %in% labels, ])))

    # Extract lower triangular part
    lower_tri <- lower.tri(network_table, diag = FALSE)
    # Get non-zero elements and their indices
    non_zero <- which(lower_tri & network_table != 0, arr.ind = TRUE)
    # Create the new table
    connect <- data.frame(
      from = rownames(network_table)[non_zero[, 1]],
      to = colnames(network_table)[non_zero[, 2]],
      score = network_table[non_zero])

    # create a vertices data.frame. One line per object of our hierarchy
    vertices  <-  data.frame(
      name = unique(c(as.character(edges$from), as.character(edges$to))) , 
      value = runif(length(unique(c(as.character(edges$from), as.character(edges$to))))))
    vertices$group  <-  edges$from[ match( vertices$name, edges$to ) ]

    vertices$id <- NA
    myleaves <- which(is.na(match(vertices$name, edges$from)))
    vertices$value[myleaves] <- colSums(network_table) / sum(network_table)
    nleaves <- length(myleaves)
    vertices$id[myleaves] <- seq(1:nleaves)
    vertices$angle <- 90 - 360 * vertices$id / nleaves
    vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)
    vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)
    mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )
    from  <-  match( connect$from, vertices$name)
    to  <-  match( connect$to, vertices$name)

    ggraph::ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
      geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, width=0.1, aes(colour = after_stat(index))) +
      scale_edge_colour_gradient(low = "red", high = "blue") +
      geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=name, angle = angle, hjust=hjust, colour=group), size=2, alpha=1) +
      geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=value, alpha=0.2)) +
      scale_colour_manual(values= rep( brewer.pal(14,"Paired") , 30)) +
      scale_size_continuous(range = c(0.1,10) ) +
    # set width by edge cor value?
      theme_void() +
      theme(
        legend.position="none",
        plot.margin=unit(c(0,0,0,0),"cm"),
      ) +
      expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
    ggsave(filename=paste0(plot_path,"Filtered_plot_circularized.pdf"), width = 12, height = 12, units = "in")
  }