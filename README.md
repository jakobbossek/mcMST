# mcMST: A Toolbox for the Multi-Criteria Minimum Spanning Tree Problem

[![CRAN Status Badge](http://www.r-pkg.org/badges/version/mcMST)](http://cran.r-project.org/web/packages/mcMST)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/mcMST)](http://cran.rstudio.com/web/packages/mcMST/index.html)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/mcMST?color=orange)](http://cran.rstudio.com/web/packages/mcMST/index.html)
[![Build Status](https://travis-ci.org/jakobbossek/mcMST.svg?branch=master)](https://travis-ci.org/jakobbossek/mcMST)
[![Build status](https://ci.appveyor.com/api/projects/status/f83u7suaxqxmtc80/branch/master?svg=true)](https://ci.appveyor.com/project/jakobbossek/mcmst/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/jakobbossek/mcMST/badge.svg?branch=master)](https://coveralls.io/github/jakobbossek/mcMST?branch=master)

## Introduction

It is well known, that the single-objective spanning tree problem (MST) is solvable in polynomial time, e.g., by the Prim's algorithm. However, in real-world applications, e.g., in network design, often multiple conflicting objectives have to be considered simultaneously. The multi-criteria version of the MST is NP-hard. The **mcMST** package for the [R](https://www.r-project.org) programming language contains methods for benchmark instance generation of multi-objective graph problems and methods for solving the multi-criteria spanning tree problem (mcMST).

Key features of the **mcMST** package are:

* A multi-objective version of Prim's algorithm.
* Evolutionary multi-objective algorithms (based on the Prüfer-encoding or  direct edge list representation) with several mutation operators.
* A modular approach for benchmark problem generation. 

## Example

Here we first generate a bi-criteria graph problem with n = 25 nodes. The first objective is the euclidean distance of node coordinates in the euclidean plane [0, 10] x [0, 10]. The second objective follows a normal distribution (N(5, 1.5)). 
```r
set.seed(1)
g = mcGP(lower = 0, upper = 10)
g = addCoordinates(g, n = 25, generator = coordUniform)
g = addWeights(g, method = "euclidean")
g = addWeights(g, method = "random", weight.fun = rnorm, mean = 5, sd = 1.5)
print(g)
```

Next, we apply the genetic algorithm proposed by Zhou & Gen with population size `mu = 10` and number of offspring `lambda = 10` for `max.iter = 100` generations.
```r
library(ggplot2)
res = mcMSTEmoaZhou(g, mu = 10L, lambda = 10L, max.iter = 100L)
ecr::plotFront(res$pareto.front)
```
See the package vignettes for more details.

## Installation Instructions

Install the [CRAN](http://cran.r-project.org) release version via:
```r
install.packages("mcMST")
```
If you are interested in trying out and playing around with the current development version use the [devtools](https://github.com/hadley/devtools) package and install directly from GitHub:

```r
install.packages("devtools", dependencies = TRUE)
devtools::install_github("jakobbossek/mcMST")
```

## Contributing to mcMST

If you encounter problems using this software, e.g., bugs or insufficient/misleading documentation, or you simply have a question, feel free to open an issue in the [issue tracker](https://github.com/jakobbossek/mcMST/issues).
In order to reproduce potential problems, please provide a minimal and reproducible code example.

Contributions to this software package are welcome via [pull requests](https://help.github.com/articles/about-pull-requests/) and will be merged at the sole discretion of the author. 

## Related work

The following packages provide methods to solve the __single-objective__ MST problem:

* [vegan: Community Ecology Package](https://cran.r-project.org/package=vegan)
* [igraph: Network Analysis and Visualization](https://cran.r-project.org/package=igraph)

Several packages implement methods to solve multi-criteria optimization problems in general:

* [ecr: Evolutionary Computation in R](https://cran.r-project.org/package=ecr)
* [mopsocd: Multi-Objective Particle Swarm Optimization with Crowding Distance](https://cran.r-project.org/package=mopsocd)
* [moko: Multi-Objective Kriging Optimization](https://cran.r-project.org/package=moko)
* [mco: Multiple Criteria Optimization Algorithms and Related Functions](https://cran.r-project.org/package=mco)



