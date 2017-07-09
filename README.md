# mcMST: Multi-Criteria Spanning-Tree-Problem

[![CRAN Status Badge](http://www.r-pkg.org/badges/version/mcMST)](http://cran.r-project.org/web/packages/mcMST)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/mcMST)](http://cran.rstudio.com/web/packages/mcMST/index.html)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/mcMST?color=orange)](http://cran.rstudio.com/web/packages/mcMST/index.html)
[![Build Status](https://travis-ci.org/jakobbossek/mcMST.svg?branch=master)](https://travis-ci.org/jakobbossek/mcMST)
[![Build status](https://ci.appveyor.com/api/projects/status/eu0nns2dsgocwntw/branch/master?svg=true)](https://ci.appveyor.com/project/jakobbossek/mcMST/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/jakobbossek/mcMST/badge.svg?branch=master)](https://coveralls.io/github/jakobbossek/mcMST?branch=master)

The **mcMST** package for the statistical programming language [R](https://www.r-project.org) contains methods for benchmark instance generation of multi-objective graph problems and methods for solving the multi-criteria spanning tree problem (mcMST).

## Example

Here we first generate a bi-criteria graph problem with n = 25 nodes. The first objective is the euclidean distance of node coordinates in [0, 10] x [0, 10] in the euclidean plane. The second objective follows a normal distribution (N(5, 1.5)). The instance generation process is modular and thus highly flexible.
```r
set.seed(1)
g = mcGP(lower = 0, upper = 10)
g = addCoordinates(g, n = 25, generator = coordUniform)
g = addWeights(g, method = "euclidean")
g = addWeights(g, method = "random", weight.fun = rnorm, mean = 5, sd = 1.5)
print(g)
```

Next, we apply the Genetic Algorithm proposed by Zhou & Gen with population size `mu = 10`, number of offspring `lambda = 10` for `max.iter = 100` generations.
```r
library(ggplot2)
res = mcMSTEmoaZhou(g, mu = 10L, lambda = 10L, max.iter = 1000L)
ecr::plotFront(res$pareto.front)
```

## Installation Instructions

The package will be available at [CRAN](http://cran.r-project.org) soon. Install the release version via:
```r
install.packages("mcMST")
```
If you are interested in trying out and playing around with the current github developer version use the [devtools](https://github.com/hadley/devtools) package and type the following command in R:

```r
install.packages("devtools", dependencies = TRUE)
devtools::install_github("jakobbossek/mcMST")
```

## Contributing to mcMST

If you encounter problems using this software, e.g., bugs or insufficient/misleading documentation, or you simply have a question, feel free to open an issue in the [issue tracker](https://github.com/jakobbossek/mcMST/issues).
In order to reproduce potential problems, please provide a minimal and reproducible code example.

Contributions to this software package are welcome via [pull requests](https://help.github.com/articles/about-pull-requests/) and will be merged at the sole discretion of the author. 



