#' mcMST: A Toolbox for the Multi-Criteria Minimum Spanning Tree Problem.
#'
#' The \pkg{mcMST} package provides a set of algorithms to
#' approximate the Pareto-optimal set/front of multi-criteria
#' minimum spanning tree (mcMST) problems.
#'
#' @section Algorithms:
#' Currently, the following algorithms are included:
#' \describe{
#'   \item{mcPrim}{A multi-criteria version of Prim's algorithm
#'   for the single-objective MST (see [1]).}
#'   \item{ZhouEmoa}{Evolutionary multi-objective algorithm operating
#'   on the Pruefer-encoding as proposed by Zhou and Gen [2].}
#'   \item{BGEmoa}{Evolutionary multi-objective algorithm operating
#'   on a direct edge list encoding. This algorithm applies a sub-tree
#'   based mutation operator as proposed by Bossek and Grimme [3].}
#'   \item{Exhaustive Enumeration}{A simple method to enumerate all Pareto-optimal
#'   solutions of a given combinatorial problem. This method is not limited to
#'   mcMST problems.}
#' }
#'
#' @references
#' [1] Knowles, J. D., and Corne, D. W. 2001. A Comparison of Encodings
#' and Algorithms for Multiobjective Minimum Spanning Tree Problems.
#' In Proceedings of the 2001 Congress on Evolutionary Computation (Ieee
#' Cat. No.01TH8546), 1:544–51 vol. 1. doi:10.1109/CEC.2001.934439.
#'
#' [2] Zhou, G., and Gen, M. 1999. Genetic Algorithm Approach on
#' Multi-Criteria Minimum Spanning Tree Problem. European Journal of
#' Operational Research 114 (1): 141–52.
#' doi:https://doi.org/10.1016/S0377- 2217(98)00016-2.
#'
#' [3] Bossek, J., and Grimme, C. 2017. A Pareto-Beneficial Sub-Tree Mutation
#' for the Multi-Criteria Minimum Spanning Tree Problem. In Proceedings of the
#' 2017 IEEE Symposium Series on Computational Intelligence. (accepted)
#'
#' @docType package
#' @name mcMST-package
#' @rdname mcMST-package
NULL
