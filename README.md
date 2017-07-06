# rmoco: Multi-Objective Combinatorial Optimization

[![CRAN Status Badge](http://www.r-pkg.org/badges/version/rmoco)](http://cran.r-project.org/web/packages/rmoco)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/rmoco)](http://cran.rstudio.com/web/packages/rmoco/index.html)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/rmoco?color=orange)](http://cran.rstudio.com/web/packages/rmoco/index.html)
[![Build Status](https://travis-ci.org/jakobbossek/rmoco.svg?branch=master)](https://travis-ci.org/jakobbossek/rmoco)
[![Build status](https://ci.appveyor.com/api/projects/status/eu0nns2dsgocwntw/branch/master?svg=true)](https://ci.appveyor.com/project/jakobbossek/rmoco/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/jakobbossek/rmoco/badge.svg?branch=master)](https://coveralls.io/github/jakobbossek/rmoco?branch=master)
[![Research software impact](http://depsy.org/api/package/cran/rmoco/badge.svg)](http://depsy.org/package/r/rmoco)

The **rmoco** package for the statistical programming language [R](https://www.r-project.org) contains methods for benchmark instance generation of multi-objective graph problems and methods for solving the multi-criteria Spanning Tree problem (mcMST).

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

Next, we apply the Genetic Algorithm proposed by Zhou & Gen with population size `mu = 10`, number of offspring `lambda = 10` for `max.iter = 100` generations. Additionally, a reference point `ref.point` is used for tracking the dominated hypervolume.
```r
library(ggplot2)
res = emoaMST_Zhou(g, mu = 10L, lambda = 10L, max.iter = 100L, ref.point = c(1000, 1000))
ecr::plotFront(res$pareto.front)
```

## Installation Instructions

The package will be available at [CRAN](http://cran.r-project.org) soon. Install the release version via:
```r
install.packages("rmoco")
```
If you are interested in trying out and playing around with the current github developer version use the [devtools](https://github.com/hadley/devtools) package and type the following command in R:

```r
devtools::install_github("jakobbossek/rmoco")
```

## Contact

Please address questions and missing features about the **rmoco** to the author Jakob Bossek <j.bossek@gmail.com>. Found some nasty bugs? Please use the [issue tracker](https://github.com/jakobbossek/rmoco/issues) for this. Pay attention to explain the problem as good as possible. At its best you provide an example, so I can reproduce your problem quickly.



