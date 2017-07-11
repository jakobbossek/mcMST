---
title: 'mcMST: Multi-criteria minimum spanning tree'
authors:
- affiliation: 1
  name: Jakob Bossek
  orcid: 0000-0002-4121-4668
date: "11 July 2017"
output: pdf_document
bibliography: paper.bib
tags:
- R
- combinatorial optimization
- evolutionary algorithms
- multi-objective optimization
affiliations:
- index: 1
  name: University of Münster
---

# Summary

The single-objective minimum spanning tree (MST) problem is a combinatorial optimization problem known to be polynomial-time solvable, e.g., using the algorithm of Prim [@P57]. However, in real-world applications one is frequently confronted with multiple objectives which need to be minimized simultaneously. Since the objectives are usually conflicting, there is no single optimal solution to this problem. Instead the goal is the approximate the so-called Pareto-front, i.e., the set of nondominated solutions (see @D01). 

The `R` [@R] package [`mcMST`](https://github.com/jakobbossek/mcMST) provides algorithms to approximate the Pareto-front of mcMST problems. Beside a multi-objective version of Prims algorithm (see, e.g., @KC01) the package offers two evolutionary multi-objective algorithms (EMOAs). The first one is based on the Prüfer-encoding [@P18] and is very similar to the genetic algorithm proposed by @ZG99. The second EMOA operates on the edge list directly. Further, a simple and generic enumeration algorithm is included, which is useful to compute the exact front of graph problems of small instance size by exhaustive enumeration. This is particularly useful to investigate properties of Pareto-optimal solutions.

Additionally, a modular toolbox to generate multi-objective benchmark graph problems is included in the `mcMST` package. This allows for easy generation of a diverse benchmark set.

# References
