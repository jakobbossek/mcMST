% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcMSTEmoaZhou.R
\name{mcMSTEmoaZhou}
\alias{mcMSTEmoaZhou}
\title{Pruefer-EMOA for the multi-objective MST problem.}
\usage{
mcMSTEmoaZhou(
  instance,
  mu,
  lambda = mu,
  mut = mutUniformPruefer,
  selMating = ecr::selSimple,
  selSurvival = ecr::selNondom,
  ref.point = NULL,
  max.iter = 100L
)
}
\arguments{
\item{instance}{[\code{\link[grapherator]{grapherator}}]\cr
Multi-objective graph.}

\item{mu}{[\code{integer(1)}]\cr
Population size.}

\item{lambda}{[\code{integer(1)}]\cr
Number of offspring generated in each generation.
Default is \code{mu}.}

\item{mut}{[\code{ecr_mutator}]\cr
Mutation operator.
Defaults to \code{\link{mutUniformPruefer}}, i.e., each digit of the Pruefer encoding
is replaced with some probability with a random number from \eqn{V = \{1, \ldots, n\}}.}

\item{selMating}{[\code{ecr_selector}]\cr
Mating selector.
Default is \code{\link[ecr]{selSimple}}.}

\item{selSurvival}{[\code{ecr_selector}]\cr
Survival selector.
Default is \code{link[ecr]{selNondom}}.}

\item{ref.point}{[\code{numeric(n.objectives) | NULL}]\cr
Reference point for hypervolume computation used for logging.
If \code{NULL} the sum of the \eqn{n} largest edges in each objective
is used where \eqn{n} is the number of nodes of \code{instance}.
This is an upper bound for the size of each spanning tree
with \eqn{(n-1)} edges.}

\item{max.iter}{[\code{integer(1)}]\cr
Maximal number of iterations.
Default is \code{100}.}
}
\value{
[\code{\link[ecr]{ecr_result}}] List of type \code{\link[ecr]{ecr_result}}
 with the following components:
 \describe{
  \item{task}{The \code{ecr_optimization_task}.}
  \item{log}{Logger object.}
  \item{pareto.idx}{Indizes of the non-dominated solutions in the last population.}
  \item{pareto.front}{(n x d) matrix of the approximated non-dominated front where n
  is the number of non-dominated points and d is the number of objectives.}
  \item{pareto.set}{Matrix of decision space values resulting with objective values
  given in pareto.front.}
  \item{last.population}{Last population.}
  \item{message}{Character string describing the reason of termination.}
 }
}
\description{
Evolutionary multi-objective algorithm to solve the
multi-objective minimum spanning tree problem. The algorithm adopts the
so-called Pruefer-number as the encoding for spanning trees. A Pruefer-number
for a graph with nodes \eqn{V = \{1, \ldots, n\}} is a sequence of \eqn{n - 2}
numbers from \eqn{V}. Cayleys theorem states, that a complete graph width n nodes
has exactly \eqn{n^{n-2}} spanning trees.
The algorithm uses mutation only: each component of an individual is replaced
uniformly at random with another node number from the node set.
}
\references{
Zhou, G. and Gen, M. Genetic Algorithm Approach on Multi-Criteria
Minimum Spanning Tree Problem. In: European Journal of Operational Research (1999).
}
\seealso{
Mutator \code{\link{mutUniformPruefer}}

Other mcMST EMOAs: 
\code{\link{mcMSTEmoaBG}()}

Other mcMST algorithms: 
\code{\link{mcMSTEmoaBG}()},
\code{\link{mcMSTPrim}()}
}
\concept{mcMST EMOAs}
\concept{mcMST algorithms}
