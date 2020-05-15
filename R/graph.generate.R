
graph.generate <- function(p, graph = "hub", m = 1, neighb = 2) {
  ##############################################################################
  # - p          An integer giving the dimension. p should be >= 2.
  # - graph      character. Specify the network structure.
  # - m          integer. The number of blocks.
  # - neighb     interger. The neighborhood within which the vertices of the lattice
  #              will be connected if graph is "Watt-Strogatz" or "small-world".
  #              The number of off-diagonal bands used if graph is "banded".
  ##############################################################################

  graph <- match.arg(graph,c("star", "hub", "banded", "scale-free", "Barabasi-Albert",
                             "small-world", "Watts-Strogatz", "random","Erdos-Renyi"))

  if (graph == "hub" | graph == "star") {

    generator <- function(p) {
      adj <- diag(p)
      adj[1, ] <- adj[, 1] <- 1
      adj
    }

  } else if (graph == "banded") {

    generator <- function(p) {
      if (neighb > p) {
        stop("The number of bands cannot exceed the dimension of each block")
      }
      # adj <- 1/(abs(outer(1:p,1:p,"-"))+1)
      adj <- matrix(1, nrow=p, ncol=p)
      adj[abs(row(adj)-col(adj)) > neighb] <- 0
      adj
    }

  } else if (graph == "scale-free" | graph == "Barabasi-Albert") {

    generator <- function(p) {
      G <- sample_pa(p, power = 1, directed = FALSE)
      adj <- get.adjacency(G, sparse = FALSE) + diag(p)
      adj
    }

  } else if (graph == "small-world" | graph == "Watts-Strogatz") {

    generator <- function(p) {
      G <- sample_smallworld(1, p, neighb, 0.05)
      adj <- get.adjacency(G, sparse = FALSE) + diag(p)
      adj
    }

  } else if (graph == "random" | graph == "Erdos-Renyi") {

    generator <- function(p) {
      G <- sample_gnp(p, 1/p)
      adj <- get.adjacency(G, sparse = FALSE) + diag(p)
      adj
    }

  }

  # Construct the block split
  ind <- split(seq(p), ceiling(m*seq(p)/p))

  # Fill in blocks to construct full network structure
  net <- diag(p)
  for (i in ind) {
    net[i, i] <- generator(length(i))
  }
  net
}

