# mcMST: A Toolbox for the Multi-Criteria Minimum Spanning Tree Problem

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00374/status.svg)](https://doi.org/10.21105/joss.00374)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/mcMST)](https://cran.r-project.org/package=mcMST)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/mcMST)](https://cran.r-project.org/package=mcMST)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/mcMST?color=orange)](https://cran.r-project.org/package=mcMST)
[![R-CMD-check](https://github.com/jakobbossek/smoof/workflows/R-CMD-check/badge.svg)](https://github.com/jakobbossek/mcMST/actions)
[![Coverage Status](https://coveralls.io/repos/github/jakobbossek/mcMST/badge.svg?branch=master)](https://coveralls.io/github/jakobbossek/mcMST?branch=master)

## Introduction

It is well known, that the single-objective spanning tree problem (MST) is solvable in polynomial time, e.g., by the Prim's algorithm. However, in real-world applications, e.g., in network design, often multiple conflicting objectives have to be considered simultaneously. The multi-criteria version of the MST is NP-hard. The **mcMST** package for the [R](https://www.r-project.org) programming language contains several methods for solving the multi-criteria spanning tree problem (mcMST).

Key features of the **mcMST** package are:

* A multi-objective version of Prim's algorithm.
* Evolutionary multi-objective algorithms (based on the Prüfer-encoding or  direct edge list representation) with several mutation operators.
* A modular approach for benchmark problem generation (**outsourced to R package [grapherator](https://github.com/jakobbossek/grapherator)**.)

## Example

Here we first generate a bi-criteria graph problem with n = 25 nodes with [grapherator](https://github.com/jakobbossek/grapherator). The first objective is the Euclidean distance of node coordinates in the euclidean plane [0, 10] x [0, 10]. The second objective follows a normal distribution (N(5, 1.5)). 
```r
library(grapherator)
set.seed(1)
g = graph(lower = 0, upper = 100)
g = addNodes(g, n = 25, generator = addNodesUniform)
g = addEdges(g, generator = addEdgesDelauney)
g = addWeights(g, generator = addWeightsDistance, method = "euclidean")
g = addWeights(g, generator = addWeightsRandom, method = rnorm, mean = 5, sd = 1.5)
print(g)
```

Next, we apply the multi-objective evolutionary algorithm proposed by Bossek & Grimme with population size `mu = 10` and number of offspring `lambda = 10` for `max.iter = 100` generations.
```r
library(ggplot2)
res = mcMSTEmoaBG(g, mu = 30L, max.iter = 300L)
ecr::plotFront(res$pareto.front)
```
See the package vignettes for more details.

## Installation Instructions

Install the [CRAN](https://cran.r-project.org) release version via:
```r
install.packages("mcMST")
```
If you are interested in trying out and playing around with the current development version use the [devtools](https://github.com/r-lib/devtools) package and install directly from GitHub:

```r
install.packages("devtools", dependencies = TRUE)
devtools::install_github("jakobbossek/mcMST")
```

## Contributing to mcMST

If you encounter problems using this software, e.g., bugs or insufficient/misleading documentation, or you simply have a question, feel free to open an issue in the [issue tracker](https://github.com/jakobbossek/mcMST/issues).
In order to reproduce potential problems, please provide a minimal and reproducible code example.

Contributions to this software package are welcome via [pull requests](https://help.github.com/articles/about-pull-requests/) and will be merged at the sole discretion of the author. 

## Publications

The following publications are strongly related to the package.

Bossek, J., & Grimme, C. (2017). A Pareto-Beneficial Sub-Tree Mutation for the Multi-Criteria Minimum Spanning Tree Problem. In Proceedings of the IEEE Symposium Series on Computational Intelligence, Honolulu, Hawai.

Bossek, J. (2017). mcMST: A Toolbox for the Multi-Criteria Minimum Spanning Tree Problem. The Journal of Open Source Software, 2017.

## Related work

The following packages provide methods to solve the __single-objective__ MST problem:

* [vegan: Community Ecology Package](https://cran.r-project.org/package=vegan)
* [igraph: Network Analysis and Visualization](https://cran.r-project.org/package=igraph)

Several packages implement methods to solve multi-criteria optimization problems in general:

* [ecr: Evolutionary Computation in R](https://cran.r-project.org/package=ecr)
* [mopsocd: Multi-Objective Particle Swarm Optimization with Crowding Distance](https://cran.r-project.org/package=mopsocd)
* [moko: Multi-Objective Kriging Optimization](https://cran.r-project.org/package=moko)
* [mco: Multiple Criteria Optimization Algorithms and Related Functions](https://cran.r-project.org/package=mco)



