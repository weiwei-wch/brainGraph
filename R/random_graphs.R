#' Simulate N random graphs w/ same clustering and degree sequence as the input.
#'
#' \code{sim.rand.graph.par} simulates \code{N} simple random graphs with the
#' same clustering (optional) and degree sequence as the input. Essentially a
#' wrapper for \code{\link[igraph]{sample_degseq}} (or, if you want to match by
#' clustering, \code{\link{sim.rand.graph.clust}}) and
#' \code{\link{set_brainGraph_attr}}. It uses \code{\link[foreach]{foreach}} for
#' parallel processing.
#'
#' If you do not want to match by clustering, then simple rewiring of the input
#' graph is performed (the number of rewire's equaling the larger of \code{1e4}
#' and \eqn{10 \times m}, where \eqn{m} is the graph's edge count).
#'
#' @param g An \code{igraph} graph object
#' @param N Integer; the number of random graphs to simulate (default: 100)
#' @param clustering Logical; whether or not to control for clustering (default:
#'   \code{FALSE})
#' @param ... Other parameters (passed to \code{\link{sim.rand.graph.clust}})
#' @export
#'
#' @return \code{sim.rand.graph.par} - a \emph{list} of \emph{N} random graphs
#'   with some additional vertex and graph attributes
#'
#' @name RandomGraphs
#' @aliases sim.rand.graph.par
#' @rdname random_graphs
#'
#' @seealso \code{\link[igraph]{rewire}, \link[igraph]{sample_degseq},
#' \link[igraph]{keeping_degseq}}
#' @examples
#' \dontrun{
#' rand1 <- sim.rand.graph.par(g[[1]][[N]], N=1e3)
#' rand1.cl <- sim.rand.graph.par(g[[1]][[N]], N=1e2,
#'   clustering=T, max.iters=1e3)
#' }

sim.rand.graph.par <- function(g, N=100, clustering=FALSE, ...) {
  stopifnot(is_igraph(g))
  iters <- max(10*ecount(g), 1e4)
  if (isTRUE(clustering)) {
    r <- foreach(i=seq_len(N), .packages=c('igraph', 'brainGraph')) %dopar% {
      tmp <- sim.rand.graph.clust(g, rewire.iters=iters, ...)
      tmp <- set_brainGraph_attr(tmp, rand=TRUE)
      tmp
    }
  } else {
    if (is_connected(g)) {
      V(g)$degree <- degree(g)
      r <- foreach(i=seq_len(N)) %dopar% {
        tmp <- sample_degseq(V(g)$degree, method='vl')
        tmp <- set_brainGraph_attr(tmp, rand=TRUE)
        tmp
      }
    } else {
      r <- foreach(i=seq_len(N)) %dopar% {
        tmp <- rewire(g, keeping_degseq(loops=FALSE, iters))
        tmp <- set_brainGraph_attr(tmp, rand=TRUE)
        tmp
      }
    }
  }
}

#' Simulate a random graph with given degree sequence and clustering.
#'
#' \code{sim.rand.graph.clust} simulates a random graph with a given degree
#' sequence \emph{and} clustering coefficient. Increasing the \code{max.iters}
#' value will result in a closer match of clustering with the observed graph.
#'
#' @param rewire.iters Integer; number of rewiring iterations for the initial
#'   graph randomization (default: 1e4)
#' @param cl The clustering measure (default: transitivity)
#' @param max.iters The maximum number of iterations to perform; choosing a
#'   lower number may result in clustering that is further away from the
#'   observed graph's (default: 100)
#' @export
#'
#' @return \code{sim.rand.graph.clust} - A single \code{igraph} graph object
#'
#' @aliases sim.rand.graph.clust
#' @rdname random_graphs
#'
#' @seealso \code{\link[igraph]{transitivity}}
#'
#' @references Bansal S., Khandelwal S., Meyers L.A. (2009) \emph{Exploring
#' biological network structure with clustered random networks}. BMC
#' Bioinformatics, 10:405-421.

sim.rand.graph.clust <- function(g, rewire.iters=1e4, cl=g$transitivity, max.iters=100) {
  g <- rewire(g, keeping_degseq(loops=FALSE, rewire.iters))
  g.cand <- g
  A <- as_adj(g.cand, sparse=FALSE, names=FALSE)
  degs <- colSums(A)
  degs.large <- which(degs > 1)

  cur.iter <- 0
  while ((transitivity(g) < cl) & (cur.iter < max.iters)) {
    repeat {
      g.cand <- g
      A <- as_adj(g.cand, sparse=FALSE, names=FALSE)

      # If E(y1, y2) and E(z1, z2) don't exist, rewire 2 edges
      repeat {
        e <- choose.edges(A, degs.large)
        if ((A[e$y1, e$y2] == 0) && (A[e$z1, e$z2] == 0)) break
      }

      A[e$y1, e$z1] <- A[e$z1, e$y1] <- A[e$y2, e$z2] <- A[e$z2, e$y2] <- 0
      A[e$y1, e$y2] <- A[e$y2, e$y1] <- A[e$z1, e$z2] <- A[e$z2, e$z1] <- 1
      g.cand <- graph_from_adjacency_matrix(A, mode='undirected')

      if (transitivity(g.cand) > transitivity(g)) break
    }

    cur.iter <- cur.iter + 1
    g <- g.cand
  }

  return(g)
}

#' Select edges for re-wiring.
#'
#' \code{choose.edges} selects edges to be re-wired when simulating random
#' graphs while controlling for \emph{clustering}. It is based on the algorithm
#' by Bansal et al. (2009), BMC Bioinformatics.
#'
#' @param A Numeric (adjacency) matrix
#' @param degs.large Integer vector of vertex numbers with degree greater than
#'   one
#'
#' @return A data frame with four elements; two edges \code{(y1, z1)} and
#'   \code{(y2, z2)} will be removed, and two edges \code{(y1, y2)} and
#'   \code{(z1, z2)} will be added between the four vertices.
#'
#' @keywords internal
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Bansal S., Khandelwal S., Meyers L.A. (2009) \emph{Exploring
#'   biological network structure with clustered random networks}. BMC
#'   Bioinformatics, 10:405-421.

choose.edges <- function(A, degs.large) {
  # Uniformly select both a random node with degree > 1 (x), and 2 of its neighbors (y1 & y2)
  #-----------------------------------------------------------------------------
  get_random_v <- function(A, degs.large) {
    repeat {
      x <- degs.large[sample.int(length(degs.large), 1)]
      neighb <- intersect(which(A[x, ] == 1), degs.large)
      if (length(neighb) >= 2) return(list(x, neighb[sample.int(length(neighb), 2)]))
    }
  }
  #-----------------------------------------------------------------------------
  repeat {
    tmp <- get_random_v(A, degs.large)
    x <- tmp[[1]]
    y <- tmp[[2]]

    # Try to select random neighbors (z1 & z2) of y1 & y2 s.t. z1 != z2 != x
    y1.neighb <- which(A[y[1], ] == 1)
    choices1 <- setdiff(y1.neighb, c(x, y[2]))
    if (length(choices1) == 0) {
      next
    } else {
      z1 <- choices1[sample.int(length(choices1), 1)]
    }

    y2.neighb <- which(A[y[2], ] == 1)
    choices2 <- setdiff(y2.neighb, c(x, y[1], z1))
    if (length(choices2) > 0) break
  }

  z2 <- choices2[sample.int(length(choices2), 1)]
  return(data.frame(y1=y[1], y2=y[2], z1, z2))
}
