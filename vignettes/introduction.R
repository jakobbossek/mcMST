## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>", warning = FALSE, message = FALSE)
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ---- fig.width=8, out.width='100%', fig.cap='Plot of the bi-objective sample graph `g.'----
library(mcMST)
library(grapherator)

set.seed(1) # reproducability
g = graph(0, 10)
g = addNodes(g, n = 25, generator = addNodesUniform)
g = addWeights(g, generator = addWeightsDistance, method = "euclidean")
g = addWeights(g, generator = addWeightsRandom, method = rnorm, mean = 5, sd = 1.5)
print(g)
do.call(gridExtra::grid.arrange, c(plot(g), list(nrow = 1L)))

