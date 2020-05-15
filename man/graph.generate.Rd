\name{graph.generate}
\alias{graph.generate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate graph
}
\description{
Generate one adjacency matrix according to graphical model.
}
\usage{
graph.generate(p, graph = "hub", m = 1, neighb = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
an integer giving the dimension. \code{p} should be >= 2.
}
  \item{graph}{
character. Specify the network structure. Possible network structure are \code{"star"} (default), \code{"hub", "banded", "scale-free", "Barabasi-Albert", "small-world", "Watts-Strogatz", "random", "Erdos-Renyi"}.
}
  \item{m}{
integer. The number of blocks.
}
  \item{neighb}{
interger. The neighborhood within which the vertices of the latticewill be connected if graph is "Watt-Strogatz" or "small-world".
The number of off-diagonal bands used if graph is "banded".
}
}
\details{
Generate a \code{"small-world"} network according to the Watts-Strogatz network model, \code{"scale-free"} network according to the Barabasi-Albert model, \code{"random"} network according to the Erdos-Renyi model.

}
\value{
A \code{p*p} matrix.
}

\examples{

##------ not run ------
adj <- graph.generate(20, graph = "scale-free", m = 2)
adj

}

