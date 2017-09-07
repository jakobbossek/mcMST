---
title: 'mcMST: A Toolbox for the Multi-Criteria Minimum Spanning Tree Problem'
authors:
- affiliation: 1
  name: Jakob Bossek
  orcid: 0000-0002-4121-4668
date: "07 September 2017"
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

The `R` [@R] package [`mcMST`](https://github.com/jakobbossek/mcMST) provides algorithms to approximate the Pareto-front of mcMST problems. Beside a multi-objective version of Prims algorithm (see, e.g., @KC01) the package offers evolutionary multi-objective algorithms (EMOAs), e.g., an Prüfer-encoding [@P18] based EMOA as proposed by @ZG99 and an EMOA based on direct encoding and a Pareto-beneficial subtree mutation operator [@BG17]. Further, a simple and generic enumeration algorithm is included, which is useful to compute the exact front of graph problems (not limited to mcMST) of small instance size by exhaustive enumeration. This is particularly useful to investigate properties of Pareto-optimal solutions.

Additionally, a modular toolbox to generate multi-objective benchmark graph problems is included in the `mcMST` package. This allows for easy generation of diverse benchmark sets.

# References
