min_clusters <- function(nz_weights)
{
    n = max(nz_weights)
    edges = as.vector(t(nz_weights))
    g = igraph::make_empty_graph(n = n) %>%
        igraph::add_edges(edges)
    return(clusters(g)$no)
}
