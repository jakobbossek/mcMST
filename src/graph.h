#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <stack>
#include <limits>
#include <random>
#include <chrono>
#include <assert.h>
#include "debug.h"
#include "UnionFind.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifndef GRAPH
#define GRAPH

typedef std::pair< int, std::vector<double> > Edge;
typedef std::pair< std::pair<int, int>, std::vector<double> > Edge2;
typedef std::vector< std::vector <Edge> > AdjacencyList;

std::vector<int> getNondominatedPoints(std::vector<std::vector<double>> points);

class Graph {
public:
  Graph(int V, int W) {
    this->V = V;
    this->E = 0;
    this->W = W;
    // Pay attention: we add a dummy element at position 0
    for (int i = 0; i <= V; ++i) {
      degrees.push_back(0);
      adjList.push_back(std::vector<Edge>());
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rngGenerator(seed);
    this->rngGenerator = rngGenerator;
  }

  int getV() const {
    return this->V;
  }

  int getE() const {
    return this->E;
  }

  int getW() const {
    return this->W;
  }

  unsigned int getDegree(int v) const {
    assert(v >= 1 && v <= this->V);
    return this->degrees[v];
  }

  std::vector<unsigned int> getDegrees() const {
    return this->degrees;
  }

  unsigned int getMaximumDegree() const {
    int maxDegree = this->degrees[1];
    for (int i = 2; i <= this->getV(); ++i) {
      if (this->degrees[i] > maxDegree) {
        maxDegree = this->degrees[i];
      }
    }
    return maxDegree;
  }

  unsigned int getNumberOfLeafs() const {
    int nLeafs = 0;
    for (int i = 1; i <= this->getV(); ++i) {
      if (this->degrees[i] == 1) {
        nLeafs += 1;
      }
    }
    return nLeafs;
  }

  unsigned int getDiameter() const {
    int V = this->getV();
    std::vector<std::vector<int>> distances(V + 1, std::vector<int>(V + 1));

    // for each node v in V ...
    for (int v = 1; v <= V; ++v) {
      // ... start BFS from v
      std::vector<bool> visited(V + 1);
      std::vector<int> queue = {v};
      distances[v][v] = 0;
      visited[v] = TRUE;

      while (!queue.empty()) {
        // get topmost element from queue
        int w = queue.back();
        queue.pop_back();

        // go through adjacency list and eventually add neighbours to queue
        for (auto edge: this->adjList[w]) {
          int u = edge.first;
          if (!visited[u]) {
            visited[u] = TRUE;
            distances[v][u] = distances[v][w] + 1L;
            queue.push_back(u);
          }
        }
      }
    } // BFS

    // now search for maximum distance
    int diameter = -1;
    for (int v = 1; v <= V; ++v) {
      for (int w = 1; w <= V; ++w) {
        if (distances[v][w] > diameter) {
          diameter = distances[v][w];
        }
      }
    }

    return diameter;
  }

  unsigned int getRandomNumber(unsigned int maxInt) {
    std::uniform_int_distribution<int> distribution(1, maxInt);
    return(distribution(this->rngGenerator));
  }

  unsigned int getRandomEdgeId() {
    return(this->edgeDistribution(this->rngGenerator));
  }

  void addEdge(int u, int v, std::vector<double> weights, int checkExists = 1) {
    assert(u >= 1 && u <= this->V);
    assert(v >= 1 && u <= this->V);
    assert(checkExists == 0 || checkExists == 1);

    if ((checkExists == 0) || (!this->hasEdge(u, v))) {
      this->adjList[u].push_back({v, weights});
      this->adjList[v].push_back({u, weights});
      this->degrees[u] += 1;
      this->degrees[v] += 1;
      this->E += 1;
    }
  }

  bool hasEdge(int u, int v) const {
    for (int i = 0; i < this->adjList[u].size(); ++i) {
      if (this->adjList[u][i].first == v)
        return true;
    }
    return false;
  }

  void saveVectorOfEdges() {
    this->edgeVector = this->getEdges();
    //std::cout << "Saved vector of edges." << std::endl;
  }

  void setEdgeProbabilities(std::vector<double> probs) {
    // assert that there is a weight for each edge
    assert(this->getE() == probs.size());

    std::discrete_distribution<int> edgeDistribution (probs.begin(), probs.end());
    this->edgeDistribution = edgeDistribution;
  }

  std::vector<double> getEdgeProbabilities() {
    return this->edgeDistribution.probabilities();
  }

  Edge2 getEdge(int u, int v) const {
    Edge2 e;
    for (Edge edge: this->adjList[u]) {
      if (edge.first == v) {
        e = {{u, v}, edge.second};
        break;
      }
    }
    return e;
  }

  std::vector<Edge2> getEdges() const {
    std::vector<Edge2> edges(this->getE());
    int edgeCounter = 0;
    for (int u = 1; u <= this->getV(); ++u) {
      for (Edge edge : this->adjList[u]) {
        int v = edge.first;
        // if (u,v) with u < v already added, then (v, u) can be skipped
        if (u < v) {
          edges[edgeCounter] = {{u, v}, edge.second};
          edgeCounter += 1;
        }
      }
    }
    return edges;
  }

  std::vector<double> getSumOfEdgeWeights() const {
    int W = this->W;
    std::vector<double> sum(W);

    for (std::vector<Edge> alist: this->adjList) {
      for (Edge edge: alist) {
        for (int i = 0; i < W; ++i) {
          sum[i] += edge.second[i];
        }
      }
    }

    for (int i = 0; i < W; ++i) {
      sum[i] /= 2;
    }

    return sum;
  }

  static Graph getIntersectionGraph(const Graph &g1, const Graph &g2) {
    // sanity checks
    assert(g1.getV() == g2.getV());
    assert(g1.getW() == g2.getW());

    int V = g1.getV();

    // std::cout << "Intersection of graphs with " << V << " nodes" << std::endl;
    // bare skeleton
    Graph g(V, g1.getW());

    // std::cout << "so far" << std::endl;

    // now iterate over all edges
    // Should be possible in O(|E|)
    unsigned int u;
    for (u = 1; u <= V; ++u) {
      // std::cout << "u = " << u << std::endl;

      for (int j = 0; j < g1.adjList[u].size(); ++j) {
        // std::cout << "j = " << j << std::endl;

        unsigned int v = g1.adjList[u][j].first;
        if (g2.hasEdge(u, v)) {
          std::vector<double> weight = g1.adjList[u][j].second;
          g.addEdge(u, v, weight);
        }
      }
    }
    // std::cout << "Finalized" << std::endl;
    return g;
  }

  static Graph getUnionGraph(const Graph &g1, const Graph &g2, unsigned int maxDegree = INT_MAX) {
    // sanity checks
    assert(g1.getV() == g2.getV());
    assert(g1.getW() == g2.getW());

    int V = g1.getV();

    // std::cout << "Intersection of graphs with " << V << " nodes" << std::endl;
    // bare skeleton
    Graph g(g1);
    // std::cout << "so far" << std::endl;

    // now iterate over all edges
    // Should be possible in O(|E|)
    unsigned int u;
    for (u = 1; u <= V; ++u) {
      // std::cout << "u = " << u << std::endl;
      for (int j = 0; j < g2.adjList[u].size(); ++j) {
        // std::cout << "j = " << j << std::endl;
        unsigned int v = g2.adjList[u][j].first;
        if ((
          !g.hasEdge(u, v))
          && ((g.getDegree(v) + 1) <= maxDegree)
          && ((g.getDegree(u) + 1) <= maxDegree)) {
          std::vector<double> weight = g1.adjList[u][j].second;
          g.addEdge(u, v, weight);
        }
      }
    }
    // std::cout << "Finalized" << std::endl;
    return g;
  }

  // static Graph importFromGrapheratorFile(std::string pathToFile) {
  //   std::ifstream infile(pathToFile, std::ios::in);

  //   char sep = ',';

  //   // read first line (nnodes, nedges, nclusters, nweights)
  //   int V, E, C, W;
  //   infile >> V >> sep >> E >> sep >> C >> sep >> W;

  //   // init graph
  //   Graph g(V, W);

  //   // now go to first character and ignore first 4 + |V| lines
  //   // (meta data + node coordinates)
  //   //FIXME: have a look ad ignore function
  //   infile.seekg(0);
  //   unsigned int linesToSkip = 4 + V;
  //   if (C > 0)
  //     linesToSkip += 1; // skip cluster membership
  //   unsigned int curLine = 1;
  //   std::string line;
  //   while (curLine <= linesToSkip) {
  //     std::getline(infile, line); // discard
  //     curLine++;
  //   }

  //   // now we are ready to read edge costs section
  //   //FIXME: generalize to >= 2 objectives
  //   int u, v;
  //   double c1, c2;

  //   while (infile.good()) {
  //     infile >> u >> sep >> v >> sep >> c1 >> sep >> c2;
  //     g.addEdge(u, v, c1, c2);
  //   }

  //   infile.close();
  //   return g;
  // }


  void removeEdge(const int u, const int v) {
    assert(u >= 1 && u <= this->V);
    assert(v >= 1 && u <= this->V);

    for (int j = 0; j < this->adjList[u].size(); ++j) {
      if (this->adjList[u][j].first == v) {
        // this is ugly!
        this->adjList[u].erase(this->adjList[u].begin() + j);
        this->E -= 1;
        this->degrees[u] -= 1;
        break;
      }
    }

    for (int j = 0; j < this->adjList[v].size(); ++j) {
      if (this->adjList[v][j].first == u) {
        // this is ugly!
        this->adjList[v].erase(this->adjList[v].begin() + j);
        this->degrees[v] -= 1;
        break;
      }
    }
  }

  bool isSpanningTree() {
    if (this->getE() != this->getV() - 1) {
      // std::cout << "Number of edges is wrong!" << std::endl;
      return false;
    }

    std::vector<std::vector<int>> comps = this->getConnectedComponents();
    if (comps.size() != 1) {
      // std::cout << "ST must have 1 conn component, but has " << comps.size() << std::endl;
    }

    return (comps.size() == 1);
  }

  void print(bool detailed = false) const {
    std::cout << "Weighted graph: n = "
              << V << ", m = "
              << E << ", p = "
              << W << std::endl;
    if (detailed) {
      // iterate over nodes
      for (unsigned int u = 1; u <= this->getV(); ++u) {
        // iterate over adjacency list
        for (unsigned int j = 0; j < this->adjList[u].size(); ++j) {
          std::cout << "c(" << u << ", " << this->adjList[u][j].first << ") = (";
          for (int i = 0; i < W; ++i) {
            std::cout << this->adjList[u][j].second[i];
          }
          std::cout << ")";
          if (j == (this->adjList[u].size() - 1)) {
            std::cout << std::endl;
          } else {
            std::cout << ", ";
          }
        }
      }
      // std::vector<int> sum = this->getSumOfEdgeWeights();
      // std::cout << "Sum of edge weights: c(" << sum[1] << ", " << sum[2] << ")" << std::endl;
    }
  }

  std::vector<int> getConnectedSubtree(int startNode, unsigned int maxNodes) {
    assert(startNode >= 1 && startNode <= this->getV());
    assert(maxNodes >= 1 && maxNodes <= this->getV());

    // we need to store node and its predecessor in BFS tree
    std::vector<int> queue;
    queue.push_back(startNode);

    std::vector<bool> done(this->getV() + 1);
    for (int i = 0; i <= this->getV(); ++i) {
      done[i] = false;
    }

    //FIXME: we know the size (maxNodes)
    std::vector<int> output;

    unsigned int nsel = 0;

    unsigned int curNode = -1;
    while (nsel < maxNodes && queue.size() > 0) {
      // get node
      curNode = queue.back();
      queue.pop_back();
      output.push_back(curNode);
      done[curNode] = true;

      // access adjList and put all neighbours into queue
      for (Edge edge: this->adjList[curNode]) {
        // skip nodes already added
        if (!done[edge.first]) {
          queue.push_back(edge.first);
        }
      }
      nsel += 1;
    }

    return output;
  }

  std::vector<std::vector<int>> getConnectedComponents() const {
    std::vector<std::vector<int>> components;

    unsigned int V = this->getV();
    std::vector<bool> visited(V + 1);

    for (int node = 1; node <= V; ++node) {
      // already visited, i.e., in some component?
      if (visited[node])
        continue;

      // otherwise perform BFS
      std::vector<int> queue = {node};
      std::vector<int> component;

      while (!queue.empty()) {
        // get topmost element from queue
        int curNode = queue.back();
        queue.pop_back();

        // mark as visited
        visited[curNode] = true;

        // add to current component
        component.push_back(curNode);

        // go through adjacency list and eventually add neighbours to queue
        for (auto edge: this->adjList[curNode]) {
          if (!visited[edge.first]) {
            queue.push_back(edge.first);
          }
        }
      }
      components.push_back(component);
    }

    return components;
  }

  static Graph getInducedSubgraph(const Graph &g, std::vector<int> nodes, int direction) {
    // copy constructor
    Graph gind(g);

    // which nodes should be kept
    std::vector<bool> keep(g.getV() + 1);
    for (auto node: nodes) {
      keep[node] = true;
    }

    //FIXME: here crap happens!
    // now go through all edges and drop, if not both endpoints should be kept
    for (unsigned int u = 1; u <= g.getV(); ++u) {
      for (unsigned int j = 0; j < g.adjList[u].size(); ++j) {
        unsigned int v = g.adjList[u][j].first;
        if (direction == 0) {
          if (!keep[u] || !keep[v]) {
            // std::cout << "Removing edge (" << u << ", " << v << ")" << std::endl;
            gind.removeEdge(u, v);
          }
        } else {
          if (!keep[u] && !keep[v]) {
            // std::cout << "Removing edge (" << u << ", " << v << ")" << std::endl;
            gind.removeEdge(u, v);
          }
        }
      }
    }

    return gind;
  }

  Graph getRandomMSTKruskal() {
    Graph initialTree(this->getV(), this->getW());
    UnionFind UF(this->V);

    std::vector<std::pair<double, std::pair<std::pair<int, int>, std::vector<double>>>> edgelist;
    for (int u = 1; u <= this->V; ++u) {
      for (int j = 0; j < this->adjList[u].size(); ++j) {
        int v = this->adjList[u][j].first;
        if (u < v) {
          double costs = (double)rand() / (double)RAND_MAX;
          edgelist.push_back({costs, {{u, v}, this->adjList[u][j].second}});
        }
      }
    }

    // sort edges in increasing order
    sort(edgelist.begin(), edgelist.end());

    // Edge iterator
    std::vector<std::pair<double, std::pair<std::pair<int, int>, std::vector<double>>>>::iterator it;
    for (it = edgelist.begin(); it != edgelist.end(); it++) {
      // get end nodes
      int u = it->second.first.first;
      int v = it->second.first.second;
      std::vector<double> weights = it->second.second;

      if (!UF.find(u, v)) {
        // link components
        initialTree.addEdge(u, v, weights);
        // merge components
        UF.unite(u, v);
      }

      // found spanning tree if number of edge is |V| - 1
      if (initialTree.getE() == (this->getV() - 1)) {
        //DEBUG("Tree has |V| - 1 edges, i.e., it is a spanning tree.");
        break;
      }
    }
    return initialTree;
  }

  Graph getMSTKruskal(std::vector<double> weights) {
    // assert(weight >= 0 && weight <= 1);

    // represent minimum spanning tree with graph object
    Graph initialTree(this->getV(), this->getW());
    UnionFind UF(this->V);

    return this->getMSTKruskal(weights, initialTree, UF);
  }

  Graph getMSTKruskal(std::vector<double> weights, Graph &initialTree) {
    // assert(weight >= 0 && weight <= 1);

    std::vector<std::vector<int>> components = initialTree.getConnectedComponents();
    //std::cout << "Building MST from " << components.size() << " components" << std::endl;
    UnionFind UF(this->V, components);
    //UF.print();

    return this->getMSTKruskal(weights, initialTree, UF);
  }

  double scalarize(std::vector<double> values, std::vector<double> lambdas) {
    double out = 0.0;
    unsigned int n = values.size();
    for (int i  = 0; i < n; ++i) {
      out += (values[i] * lambdas[i]);
    }
    return out;
  }

  Graph getMSTKruskal(std::vector<double> weights, Graph &initialTree, UnionFind &UF) {
    // assert(weight >= 0 && weight <= 1);

    // // represent minimum spanning tree with graph object
    // Graph mst(this->getV(), this->getW());

    // // init efficient set data structure
    // UnionFind UF(this->V);

    // now we need to transform graph to list of edges which can be
    // sorted by weights

    std::vector<std::pair<double, std::pair<std::pair<int, int>, std::vector<double>>>> edgelist;
    for (int u = 1; u <= this->V; ++u) {
      for (int j = 0; j < this->adjList[u].size(); ++j) {
        int v = this->adjList[u][j].first;
        if (u < v) {
          std::vector<double> costs = this->adjList[u][j].second;
          double scalcosts = this->scalarize(costs, weights);
          edgelist.push_back({scalcosts, {{u, v}, costs}});
        }
      }
    }

    // sort edges in increasing order
    sort(edgelist.begin(), edgelist.end());

    // Edge iterator
    std::vector<std::pair<double, std::pair<std::pair<int, int>, std::vector<double>>>>::iterator it;
    for (it = edgelist.begin(); it != edgelist.end(); it++) {
      // get end nodes
      int u = it->second.first.first;
      int v = it->second.first.second;
      std::vector<double> costs = it->second.second;

      if (!UF.find(u, v)) {
        // link components
        initialTree.addEdge(u, v, costs);
        // merge components
        UF.unite(u, v);
      }

      // found spanning tree if number of edge is |V| - 1
      if (initialTree.getE() == (this->getV() - 1)) {
        //DEBUG("Tree has |V| - 1 edges, i.e., it is a spanning tree.");
        break;
      }
    }

    if (!initialTree.isSpanningTree()) {
      Rcout << "Massive fail" << std::endl;
    }
    return initialTree;
  }

  std::vector<int> toBitstring(Graph &subgraph) {
    std::vector<int> bitstring(this->getE());

    std::vector<Edge2> edges = this->getEdges();
    for (int i = 0; i < edges.size(); ++i) {
      int from = edges[i].first.first;
      int to = edges[i].first.second;
      if (subgraph.hasEdge(from, to)) {
        bitstring[i] = 1;
      } else {
        bitstring[i] = 0;
      }
    }
    return bitstring;
  }

  //FIXME: sample between {1, ..., maxDrop}. At the moment we just
  // use maxDrop
  Graph getMSTBySubforestMutation(Graph &mst, unsigned int maxDrop, std::vector<double> weights) {
    assert(maxDrop <= this->getE());

    Graph forest(mst);

    // get IDs of ALL edges
    std::vector<int> selEdges(forest.getE());
    for (int i = 0; i < forest.getE(); ++i) {
      selEdges[i] = i;
    }

    // now shuffle edges
    std::random_shuffle(selEdges.begin(), selEdges.end());

    // Remove edges from tree and build forest
    //FIXME: this is very inefficient for dense graphs
    /*
    Write method removeEdges(std::vector<int> edgeIDs) which
    expects ids and removes these edges
    */
    std::vector<Edge2> edges = forest.getEdges();
    for (int i = 0; i < maxDrop; ++i) {
      int selectedEdgeID = selEdges[i];
      Edge2 edge = edges[selectedEdgeID];
      //std::cout << "Removing " << (i+1) << "-th edge" << std::endl;
      forest.removeEdge(edge.first.first, edge.first.second);
    }

    // build new spanning tree by reconnecting forest
    Graph mst2 = this->getMSTKruskal(weights, forest);

    return mst2;
  }


  std::vector<Edge2> getCircleEdges(int u) const {
    unsigned int V = this->getV();

    // circle edges
    std::vector<Edge2> circleEdges;

    // bookkeeping of visited nodes
    std::vector<bool> visited(V + 1);

    // implicit bookkeeping of DFS tree, i.e., pi[v] = w means (w, v) is DFS tree edge
    std::vector<int> pi(V + 1);

    // costs of DFS tree edges
    std::vector<std::vector<double>> wi(V + 1);

    // create stack for iterative DFS
    std::stack<int> stack;

    // put start node on stack
    stack.push(u);
    pi[u] = u;

    while (!stack.empty()) {
      u = stack.top();
      stack.pop();

      if (visited[u]) {
        continue;
      }

      visited[u] = true;

      // go through adjacency list and eventually add neighbours to queue
      for (auto edge: this->adjList[u]) {
        int v = edge.first;

        // put on stack if not already visited
        if (!visited[v]) {
          // set predecessor
          pi[v] = u;
          wi[v] = edge.second;
          stack.push(v);
        }

        // found circle! Reconstruct edges on circle
        if (visited[v] && !(v == pi[u])) {
          circleEdges.push_back({{u, v}, edge.second});
          int w;
          while (pi[u] != u) {
            w = pi[u];
            circleEdges.push_back({{w, u}, wi[u]});
            u = w;
          }
          return(circleEdges);
        }
      }
    } // while
    return circleEdges;
  }

  std::vector<Graph> filterNondominatedTrees(std::vector<Graph> graphs) const {
    std::vector<std::vector<double>> costs;
    for (Graph g: graphs) {
      costs.push_back(g.getSumOfEdgeWeights());
    }

    std::vector<int> nondomIndizes = getNondominatedPoints(costs);

    std::vector<Graph> nondomGraphs;
    for (int i: nondomIndizes) {
      nondomGraphs.push_back(graphs[i]);
    }
    return nondomGraphs;
  }

  std::vector<Edge2> getNonDominatedEdges(std::vector<Edge2> edges) const {
    std::vector<std::vector<double>> costs;
    for (Edge2 edge: edges) {
      std::vector<double> edgeCosts = edge.second;
      costs.push_back(edgeCosts);
    }

    std::vector<int> nondomIndizes = getNondominatedPoints(costs);

    std::vector<Edge2> nondomEdges;
    for (int i: nondomIndizes) {
      nondomEdges.push_back(edges[i]);
    }
    return nondomEdges;
  }

  Graph getMSTByEdgeExchange(Graph &mst, unsigned int repls, bool deleteLargest = false) {
    unsigned int V = this->getV();

    // sanity checks
    assert(repls <= V);
    assert(mst.getE() == (V - 1));

    Graph mst2(mst);

    for (unsigned int repl = 0; repl < repls; ++repl) {
      //std::cout << "Performing edge exchange " << repl + 1 << std::endl;
      int edgeToAddId = this->getRandomEdgeId();
      Edge2 edgeToAdd = this->edgeVector[edgeToAddId];
      int u = edgeToAdd.first.first;
      int v = edgeToAdd.first.second;
      std::vector<double> weights = edgeToAdd.second;
      // do nothing if edge already in tree
      if (mst2.hasEdge(u, v)) {
        // std::cout << "Selected edge already in tree." << std::endl;
        continue;
      }
      // otherwise add edge ...
      mst2.addEdge(u, v, weights);
      // and get rid of random edge on introduced circle
      std::vector<Edge2> edgesOnCircle = mst2.getCircleEdges(u);
      // std::cout << "Added edge: (" << u << ", " << v << "). There is a circle of length " << edgesOnCircle.size() << " edges now!" << std::endl;
      if (!deleteLargest) {
        int edgeToDeleteIdx = getRandomNumber(edgesOnCircle.size()) - 1;
        std::pair<int, int> edgeToDelete = edgesOnCircle[edgeToDeleteIdx].first;
        mst2.removeEdge(edgeToDelete.first, edgeToDelete.second);
      } else {
        int largestIdx = 0;
        double largestWeight = -1;
        // sample random weight
        //FIXTHIS
        double rndWeight = (double)rand() / (double)RAND_MAX;
        for (int i = 1; i < edgesOnCircle.size(); ++i) {
          std::vector<double> edgeWeight = edgesOnCircle[i].second;
          double scalWeight = this->scalarize(edgeWeight, {rndWeight, 1 - rndWeight});
          if (scalWeight > largestWeight) {
            largestWeight = scalWeight;
            largestIdx = i;
          }
        }
        std::pair<int, int> edgeToDelete = edgesOnCircle[largestIdx].first;
        mst2.removeEdge(edgeToDelete.first, edgeToDelete.second);
      }
      assert(mst2.isSpanningTree());
    }

    return(mst2);
  }

  std::vector<Graph> doMCPrim() {
    int V = this->getV();
    int W = this->getW();

    std::vector<Graph> trees;
    std::vector<Graph> trees2;

    // wlog the start node is not 1
    int startNode = 1;

    // now get all edges (1, v) which are non-dominated
    std::vector<Edge2> firstEdges;
    for (Edge2 edge : this->edgeVector) {
      if (edge.first.first == startNode || edge.first.second == startNode) {
        firstEdges.push_back(edge);
      }
    }
    std::vector<Edge2> nonDomEdges = this->getNonDominatedEdges(firstEdges);

    // now build initial partial trees, each one for each non-dominated edge
    for (Edge2 edge: nonDomEdges) {
      Graph tree(V, W);
      tree.addEdge(edge.first.first, edge.first.second, edge.second);
      trees.push_back(tree);
    }

    // now loop until partial trees are actually spanning trees
    unsigned int i = 1;
    while (i < V - 1) {
      trees2.clear();
      Rcout << "Adding " << i + 1 << "-th edge" << std::endl;
      // for each partial tree, append non-dominated edges and check
      for (Graph partialTree: trees) {
        // go through neighbors of partial tree, i.e., determine candidate edges
        std::vector<Edge2> candidateEdges;

        for (Edge2 edge: this->edgeVector) {
          int v = edge.first.first;
          int w = edge.first.second;
          if ((partialTree.getDegree(v) == 0 && partialTree.getDegree(w) > 0) ||
            (partialTree.getDegree(v) > 0 && partialTree.getDegree(w) == 0)) {
            candidateEdges.push_back(edge);
          }
        }

        // now search for non dominated edges among those selected
        nonDomEdges = getNonDominatedEdges(candidateEdges);
        //Rcout << "Found " << nonDomEdges.size() << " nondominated edges!" << std::endl;
        for (Edge2 edge: nonDomEdges) {
          // make copy of partial tree
          Graph partialTreeCopy(partialTree);
          // add nondominated edge
          partialTreeCopy.addEdge(edge.first.first, edge.first.second, edge.second);
          // check if we have a duplicate, i.e., this paertial tree was already created
          bool alreadyExists = FALSE;
          std::vector<Edge2> partialTreeEdges = partialTreeCopy.getEdges();
          for (Graph otherTree: trees2) {
            unsigned int commonEdges = 0;
            for (Edge2 partialEdge: partialTreeEdges) {
              if (!otherTree.hasEdge(partialEdge.first.first, partialEdge.first.second)) {
                break;
              }
              commonEdges += 1;
            }
            if (commonEdges == partialTreeCopy.getE()) {
              alreadyExists = TRUE;
              break;
            }
            // BACK-UP
            // Graph intersectionTree = Graph::getIntersectionGraph(otherTree, partialTreeCopy);
            // // if intersection contains the same edges -> break
            // if (intersectionTree.getE() == partialTreeCopy.getE()) {
            //   alreadyExists = TRUE;
            //   break;
            // }
          }
          // save partial tree
          if (!alreadyExists) {
            trees2.push_back(partialTreeCopy);
          }
        }
      }

      // unsigned int ntrees = trees2.size();
      // unsigned int treesize = trees2[0].getE();


      // // all (partial) trees are considered not to be dups
      // std::vector<bool> dups(ntrees, FALSE);
      // // now iterate over all trees
      // for (int i = 0; i < ntrees; ++i) {
      //   // if already marked, just skip
      //   if (dups[i]) {
      //     continue;
      //   }
      //   // otherwise iterate over all other trees
      //   for (int j = i + 1; j < ntrees; ++j) {
      //     // dups do not need to be considered again
      //     if (dups[j])
      //       continue;
      //     std::vector<Edge2> t1edges = trees2[i].getEdges();
      //     int commonEdges = 0;
      //     for (Edge2 e: t1edges) {
      //        if (!trees2[j].hasEdge(e.first.first, e.first.second)) {
      //          break;
      //        }
      //        commonEdges += 1;
      //     }
      //     // found dup
      //     if (commonEdges == treesize) {
      //       dups[j] = TRUE;
      //     }
      //   }
      // }
      // // finally filter all duplicates
      // std::vector<Graph> trees3;
      // for (int i = 0; i < ntrees; ++i) {
      //   if (!dups[i]) {
      //     trees3.push_back(trees2[i]);
      //   }
      // }
      // trees = trees3;
      trees = trees2;

      Rcout << "Now we have " << trees.size() << " partial trees!" << std::endl;
      i += 1;
    } // while

    trees = filterNondominatedTrees(trees);

    for (Graph mst: trees) {
      if (!mst.isSpanningTree()) {
        Rcout << "FAIL" << std::endl;
      }
    }
    return trees;
  }

  std::vector<Graph> doMCPrim2() {
    int V = this->getV();
    int W = this->getW();

    std::vector<Graph> trees;

    // get non-dominated edges
    std::vector<Edge2> nonDomEdges = this->getNonDominatedEdges(this->edgeVector);
    // now build initial partial trees, each one for each non-dominated edge
    for (Edge2 edge: nonDomEdges) {
      Graph tree(V, W);
      tree.addEdge(edge.first.first, edge.first.second, edge.second);
      trees.push_back(tree);
    }

    // now loop until partial trees are actually spanning trees
    unsigned int i = 1;
    while (i < V - 1) {
      Rcout << "Adding " << i + 1 << "-th edge" << std::endl;
      std::vector<Graph> trees2;
      // for each partial tree, append edges and check
      for (Graph partialTree: trees) {
        // go through neighbors of partial tree, i.e., determine candidate edges
        std::vector<Edge2> candidateEdges;

        for (Edge2 edge: this->edgeVector) {
          int v = edge.first.first;
          int w = edge.first.second;
          if ((partialTree.getDegree(v) == 0 && partialTree.getDegree(w) > 0) ||
            (partialTree.getDegree(v) > 0 && partialTree.getDegree(w) == 0)) {
            candidateEdges.push_back(edge);
          }
        }

        // now search for non dominated edges among those selected
        nonDomEdges = getNonDominatedEdges(candidateEdges);
        //Rcout << "Found " << nonDomEdges.size() << " nondominated edges!" << std::endl;
        for (Edge2 edge: nonDomEdges) {
          // make copy of partial tree
          Graph partialTreeCopy(partialTree);
          // add nondominated edge
          partialTreeCopy.addEdge(edge.first.first, edge.first.second, edge.second);
          // save partial tree
          trees2.push_back(partialTreeCopy);
        }
      }
      trees = trees2;

      // finally filter dominated partial trees
      trees = this->filterNondominatedTrees(trees);

      Rcout << "Now we have " << trees.size() << " partial trees!" << std::endl;
      i += 1;
    }

    for (Graph mst: trees) {
      if (!mst.isSpanningTree()) {
        Rcout << "FAIL" << std::endl;
      }
    }
    return trees;
  }

  // Graph getMSTBySubgraphMutation(Graph &mst, unsigned int maxSelect, bool scalarize = true) {
  //   assert(maxSelect <= this->getV());

  //   unsigned int V = this->getV();

  //   // first make a copy of input
  //   Graph forest(mst);

  //   // get random start node
  //   //FIXME: write helper getRandomInteger(max)
  //   int rndNode = (int)(((double)rand() / (double)RAND_MAX) * V) + 1;
  //   // prevent memory not mapped error
  //   if (rndNode > V) {
  //     rndNode = V;
  //   }

  //   // now we need to extract connected subtree and delete edges
  //   // we need to store node and its predecessor in BFS tree
  //   std::vector<int> queue;
  //   queue.push_back(rndNode);

  //   std::vector<bool> done(V + 1);
  //   for (int i = 0; i <= V; ++i) {
  //     done[i] = false;
  //   }

  //   //FIXME: we know the size (maxNodes)
  //   std::vector<int> nodesintree;

  //   unsigned int nsel = 0;

  //   /*
  //   IDEA:
  //   - copy mst O(V)
  //   - BFS until W nodes selected (keep boolean array
  //     w with w[i] = true if i-th node selected) O(W)
  //   - go through nodes selected an go through their
  //     adj. list in source graph, check whether neighbour
  //     selected (O(1)) and eventually add to forest.
  //   - Apply Kruskal to forest: O(W^2 log(W))


  //   */

  //   unsigned int curNode = -1;
  //   while (nsel < maxSelect && queue.size() > 0) {
  //     // get node
  //     curNode = queue.back();
  //     queue.pop_back();
  //     nodesintree.push_back(curNode);
  //     done[curNode] = true;

  //     // access adjList and put all neighbours into queue
  //     for (Edge edge: forest.adjList[curNode]) {
  //       // skip nodes already added
  //       if (!done[edge.first]) {
  //         queue.push_back(edge.first);
  //         // here we remove the edge
  //         //forest.removeEdge(curNode, edge.first);
  //       }
  //     }
  //     nsel += 1;
  //   } // while

  //   // now add all existing edges between nodes in nodesintree
  //   // to forest
  //   for (int i = 0; i < nodesintree.size(); ++i) {
  //     int u = nodesintree[i];
  //     for (Edge edge: this->adjList[u]) {
  //       int v = edge.first;
  //       if (done[v]) {
  //         forest.addEdge(u, v, edge.second.first, edge.second.second);
  //       }
  //     }
  //   }

  //   // now we need to rewire the edges
  //   double rndWeight = (double)rand() / (double)RAND_MAX;
  //   if (!scalarize) {
  //     rndWeight = round(rndWeight);
  //   }

  //   // build new spanning tree by reconnecting forest
  //   Graph mst2 = forest.getMSTKruskal(rndWeight);

  //   return mst2;
  // }

  Graph getMSTByInsertionFirstSubgraphMutation(Graph &mst, unsigned int maxAdd, std::vector<double> weights, unsigned int maxDegree) {
    // FIXME: magic number 2 (think of reasonable maximum value for maxAdd)
    assert(maxAdd >= 1);

    Graph uniong(mst);

    // Now sample edges:
    // 1) Sample adjacency list entry in O(1)
    // 2) Sample and access edge in O(1) (we know the number of edges)
    // NOTE: I have no idea why it u is always (n+1) if I pass this->getDegrees().begin() and .end() to discrete_distribution :(
    std::vector<int> probs(this->getV());
    for (int i = 0; i < this->getV(); ++i) {
      probs[i] = this->getDegree(i + 1);
    }
    std::discrete_distribution<int> distribution_adjlist(probs.begin(), probs.end());

    // for (int i = 0; i <= 3; ++i) {
    //   std::cout << "degree[" << i << "] = " << this->getDegrees()[i] << std::endl;
    //   std::cout << "RND: " << distribution_adjlist(this->rngGenerator) << std::endl;
    // }

    for (int j = 0; j < maxAdd; ++j) {
      // sample first node (+1 since for some weird reason I decided to start at 1, so 0 is a dummy)
      int u = distribution_adjlist(this->rngGenerator) + 1;
      // std::cout << "u=" << u << std::endl;

      // workaround since |this->getDegrees()|=n+1
      if (u > this->getV()) {
        u = this->getV();
      }
      std::uniform_int_distribution<int> distribution_edge(0, this->getDegree(u) - 1);
      // sample second node
      int vID = distribution_edge(this->rngGenerator);
      Edge e = this->adjList[u][vID];
      int v = e.first;
      // std::cout << "v=" << v << std::endl;

      if (!uniong.hasEdge(u, v)) {
        if ((uniong.getDegree(u) < maxDegree) && (uniong.getDegree(v) < maxDegree)) {
          uniong.addEdge(u, v, e.second);
        }
      }
    }

    Graph mst2 = uniong.getMSTKruskal(weights);

    return mst2;
  }

  Graph getMSTByUnionCrossover(Graph &mst1, Graph &mst2, std::vector<double> weights, unsigned int maxDegree  = INT_MAX) {
    // FIXME: throws "memory not mapped error". I have no nclue why
    // Graph uniong = Graph::getUnionGraph(mst1, mst2, maxDegree);

    // FIXME: copy & paste from Graph::getUnionGraph
    Graph uniong(mst1);
    int V = uniong.getV();

    // now iterate over all edges
    // Should be possible in O(|E|)
    unsigned int u;
    for (u = 1; u <= V; ++u) {
      // std::cout << "u = " << u << std::endl;
      for (int j = 0; j < mst2.adjList[u].size(); ++j) {
        // std::cout << "j = " << j << std::endl;
        unsigned int v = mst2.adjList[u][j].first;
        if ((
          !uniong.hasEdge(u, v))
          && (uniong.getDegree(v) < maxDegree)
          && (uniong.getDegree(u) < maxDegree)) {
          std::vector<double> weight = mst2.adjList[u][j].second;
          uniong.addEdge(u, v, weight);
        }
      }
    }

    Graph mst = uniong.getMSTKruskal(weights);
    return mst;
  }

  Graph getMSTBySubgraphMutation(Graph &mst, unsigned int maxSelect, std::vector<double> weights) {
    assert(maxSelect >= 2 && maxSelect <= this->getV());

    // Rcout << "maxSelect: " << maxSelect << std::endl;
    // Rcout << "scalarize: " << scalarize << std::endl;

    // get number of nodes
    unsigned int V = this->getV();

    // first make a copy of input: Theta(|V|)
    Graph mst2(mst);

    // get random start node
    int rndNode = (int)(((double)rand() / (double)RAND_MAX) * V) + 1;
    // prevent memory not mapped error
    if (rndNode > V) {
      rndNode = V;
    }

    // now we need to extract connected subtree and delete edges
    // we need to store node and its predecessor in BFS tree
    std::vector<int> queue;
    queue.push_back(rndNode);

    // initialize which nodes have already been visited
    std::vector<bool> done(V + 1);
    for (int i = 0; i <= V; ++i) {
      done[i] = false;
    }

    //FIXME: we know the size (maxNodes)
    std::vector<int> nodesintree;

    // mapping from nodesintree to 1, ..., maxSelect
    std::vector<int> revmapping(V + 1);
    for (int i = 0; i <= V; ++i) {
      revmapping[i] = 0;
    }

    unsigned int nsel = 0;

    // IDEA:
    // - copy mst O(V)
    // - BFS until W nodes selected (keep boolean array
    //   w with w[i] = true if i-th node selected) O(W)
    // - go through nodes selected an go through their
    //   adj. list in source graph, check whether neighbour
    //   selected (O(1)) and eventually add to forest.
    // - Apply Kruskal to forest: O(W^2 log(W))


    // BFS: Theta(sigma * |V|)
    unsigned int curNode = -1;
    while (nsel < maxSelect && queue.size() > 0) {
      // get node
      curNode = queue.back();
      queue.pop_back();
      nodesintree.push_back(curNode);
      done[curNode] = true;

      // access adjList and put all neighbours into queue
      for (Edge edge: mst2.adjList[curNode]) {
        // skip nodes already added
        if (!done[edge.first]) {
          queue.push_back(edge.first);
          // here we remove the edge
          //mst2.removeEdge(curNode, edge.first);
        }
      }
      nsel += 1;
    } // while

    // remove edges from copy
    for (Edge2 edge: mst2.getEdges()) {
      if (done[edge.first.first] && done[edge.first.second]) {
        mst2.removeEdge(edge.first.first, edge.first.second);
      }
    }

    // Rcout << "        T has nodes: " << mst.getV() << ", edges: "  << mst.getE() << std::endl;
    // Rcout << "Copy of T has nodes: " << mst2.getV() << ", edges: "  << mst2.getE() << std::endl;

    // build rev-mapping: (O(|V|))
    for (int i = 0; i < nodesintree.size(); ++i) {
      revmapping[nodesintree[i]] = i + 1;
    }

    // Rcout << "Nodes in tree: " << nodesintree.size() << std::endl;

    // initialize new sub-graph: O(sigma * |V|)
    Graph subgraph(maxSelect, this->getW());
    for (int i = 0; i < nodesintree.size(); ++i) {
      int curNode = nodesintree[i];
      for (Edge edge: this->adjList[curNode]) {
        // if target is also selected
        if (done[edge.first]) {
          // add mapped node
          int sourceNode = revmapping[curNode];
          int targetNode = revmapping[edge.first];
          subgraph.addEdge(sourceNode, targetNode, edge.second);
        }
      }
    }

    // Graph subgraph(V, this->getW());
    // // now add all nodes and corresponding edges from source graph into sub-graph
    // for (int i = 0; i < nodesintree.size(); ++i) {
    //   int curNode = nodesintree[i];
    //   for (Edge edge: this->adjList[curNode]) {
    //     // if the target node is also selected, add it to. subgraph
    //     if (done[edge.first]) {
    //       subgraph.addEdge(curNode, edge.first, edge.second.first, edge.second.second);
    //     }
    //   }
    // }

    // Now calculate MST on sub-graph: O((sigma + sigma^2) * log(sigma)) = O(sigma^2 log(sigma)) = O(log(|V|)^2 log(log(|V|))) if sigma = log(|V|)
    Graph subgraphMST = subgraph.getMSTKruskal(weights);

    // Rcout << "Genertated subgraph with " << subgraph.getV() << " nodes and " << subgraph.getE() << " edges " << std::endl;


    // Finally go through edges of subgraphMST and add to our copy: O(sigma)
    std::vector<Edge2> subgraphMSTEdges = subgraphMST.getEdges();
    for (Edge2 edge: subgraphMSTEdges) {
      int sourceNode = nodesintree[edge.first.first - 1];
      int targetNode = nodesintree[edge.first.second - 1];
      mst2.addEdge(sourceNode, targetNode, edge.second);
    }

    // build new spanning tree by reconnecting forest
    //Graph mst2 = this->getMSTKruskal(rndWeight, forest);

    return mst2;
  }

  // Graph getMSTBySubgraphMutation(Graph &mst, unsigned int maxSelect, bool scalarize = true) {
  //   assert(maxSelect <= this->getV());

  //   unsigned int V = this->getV();

  //   // first make a copy of input
  //   Graph forest(mst);

  //   // get random start node
  //   //FIXME: write helper getRandomInteger(max)
  //   int rndNode = (int)(((double)rand() / (double)RAND_MAX) * V) + 1;
  //   // prevent memory not mapped error
  //   if (rndNode > V) {
  //     rndNode = V;
  //   }

  //   // now we need to extract connected subtree and delete edges
  //   // we need to store node and its predecessor in BFS tree
  //   std::vector<int> queue;
  //   queue.push_back(rndNode);

  //   std::vector<bool> done(V + 1);
  //   for (int i = 0; i <= V; ++i) {
  //     done[i] = false;
  //   }

  //   //FIXME: we know the size (maxNodes)
  //   std::vector<int> nodesintree;

  //   unsigned int nsel = 0;


  //   // IDEA:
  //   // - copy mst O(V)
  //   // - BFS until W nodes selected (keep boolean array
  //   //   w with w[i] = true if i-th node selected) O(W)
  //   // - go through nodes selected an go through their
  //   //   adj. list in source graph, check whether neighbour
  //   //   selected (O(1)) and eventually add to forest.
  //   // - Apply Kruskal to forest: O(W^2 log(W))


  //   unsigned int curNode = -1;
  //   while (nsel < maxSelect && queue.size() > 0) {
  //     // get node
  //     curNode = queue.back();
  //     queue.pop_back();
  //     nodesintree.push_back(curNode);
  //     done[curNode] = true;

  //     // access adjList and put all neighbours into queue
  //     for (Edge edge: forest.adjList[curNode]) {
  //       // skip nodes already added
  //       if (!done[edge.first]) {
  //         queue.push_back(edge.first);
  //         // here we remove the edge
  //         forest.removeEdge(curNode, edge.first);
  //       }
  //     }
  //     nsel += 1;
  //   } // while

  //   // now we need to rewire the edges
  //   double rndWeight = (double)rand() / (double)RAND_MAX;
  //   if (!scalarize) {
  //     rndWeight = round(rndWeight);
  //   }

  //   // build new spanning tree by reconnecting forest
  //   Graph mst2 = this->getMSTKruskal(rndWeight, forest);

  //   return mst2;
  // }

private:
  int V;
  int E;
  int W;
  std::vector<unsigned int> degrees;
  AdjacencyList adjList;
  std::vector<Edge2> edgeVector;

  std::default_random_engine rngGenerator;
  std::discrete_distribution<int> edgeDistribution;
};
#endif

class GraphImporter {
public:
  Graph importFromGrapheratorFileCPP(CharacterVector pathToFileR) {
    std::string pathToFile = Rcpp::as<std::string>(pathToFileR);

    std::ifstream infile(pathToFile, std::ios::in);

    char sep = ',';

    // read first line (nnodes, nedges, nclusters, nweights)
    int V, E, C, W;
    infile >> V >> sep >> E >> sep >> C >> sep >> W;

    // init graph
    Graph g(V, W);

    // now go to first character and ignore first 4 + |V| lines
    // (meta data + node coordinates)
    //FIXME: have a look ad ignore function
    infile.seekg(0);
    unsigned int linesToSkip = 4 + V;
    if (C > 0)
      linesToSkip += 1; // skip cluster membership
    unsigned int curLine = 1;
    std::string line;
    while (curLine <= linesToSkip) {
      std::getline(infile, line); // discard
      curLine++;
    }

    // now we are ready to read edge costs section
    //FIXME: generalize to >= 2 objectives
    int u, v;
    std::vector<double> ws(W);

    while (infile.good()) {
      // Old "two weights" special case
      // infile >> u >> sep >> v >> sep >> c1 >> sep >> c2;
      infile >> u >> sep >> v;
      for (int i = 0; i < W; ++i) {
       infile >> sep >> ws[i];
      }
      g.addEdge(u, v, ws);
    }

    infile.close();
    return g;
  }
};

class RepresentationConverter {
public:
  Graph prueferCodeToGraph(Graph *g, IntegerVector pcode) {
    int n = pcode.size();
    int V = n + 2;

    Graph tree(V, g->getW());

    std::vector<int> degrees(V + 1);

    // initialize
    for (int i = 1; i <= V; ++i) {
      degrees[i] = 0;
    }

    // update node which are in Pruefer Code
    for (int i = 0; i < n; ++i) {
      degrees[pcode[i]] += 1;
    }

    // now build the tree
    int v = 1;
    for (int i = 0; i < n; ++i) {
      int u = pcode[i];
      for (v = 1; v <= V; ++v) {
        if (degrees[v] == 0) {
          // add edge
          std::vector<double> weight = g->getEdge(u, v).second;
          tree.addEdge(u, v, weight);

          // reduce degrees
          degrees[u] -= 1;
          degrees[v] -= 1;
          break;
        }
      }
    }

    // now two more nodes with degree 1 are left
    int u = 0;
    v = 0;
    for (int i = 1; i <= V; ++i) {
      if (degrees[i] == 0) {
        if (u == 0) {
          u = i;
        } else {
          v = i;
          break;
        }
      }
    }
    std::vector<double> weight = g->getEdge(u, v).second;
    tree.addEdge(u, v, weight);

    return tree;
  }

  Graph edgeListToGraph(Graph *g, NumericMatrix edgeList) {
    Graph g2(g->getV(), g->getW());

    for (int i = 0; i < edgeList.ncol(); ++i) {
      int u = edgeList(0, i);
      int v = edgeList(1, i);
      std::vector<double> weight = g->getEdge(u, v).second;
      g2.addEdge(u, v, weight);
    }

    return g2;
  }
};

std::vector<int> getNondominatedPoints(std::vector<std::vector<double>> points) {
  int n = points.size();
  //int m = points[1].size();

  // Ugly: make new vector and store indizes so they are ordered as well
  std::vector<std::pair<int, std::vector<double>>> points2;
  for (int i = 0; i < n; ++i) {
    points2.push_back({i, points[i]});
  }

  // now sort regarding first dimension in ascending order
  std::sort(points2.begin(), points2.end(),
    [](const std::pair<int, std::vector<double>>& x, const std::pair<int, std::vector<double>>& y) {
      // assure strict weak order (i.e., comp(x,x) == FALSE for all x)
      // Otherwise you get some awesome segfaults
      // See, e.g. https://stackoverflow.com/questions/18291620/why-will-stdsort-crash-if-the-comparison-function-is-not-as-operator
      // if (x.second[0] == y.second[0] && x.second[1] == y.second[1])
      //   return false;
      // if (x.second[0] == y.second[0]) {
      //    return (x.second[1] < y.second[1]);
      // } else {
        return (x.second[0] < y.second[0]);
      //}
    });


  std::vector<int> nondomIndizes;
  double currentX1 = points2[0].second[0];
  double bestX2 = points2[0].second[1];

  int tokeep = 0;

  for (int i = 1; i < n; ++i) {
    double X1 = points2[i].second[0];
    double X2 = points2[i].second[1];

    if (X1 == currentX1 && X2 < bestX2) {
      tokeep = i;
      bestX2 = X2;
    // keep identical points - start
    } else if (X1 == currentX1 && X2 == bestX2 && i < (n - 1) && points2[i+1].second[0] > currentX1) {
      // Solution with equal quality
      nondomIndizes.push_back(points2[tokeep].first);
      tokeep = i;
    } else if (X1 == currentX1 && X2 == bestX2 && i == (n - 1)) {
      nondomIndizes.push_back(points2[tokeep].first);
      tokeep = i;
    // keep identical points - end
    } else if (X1 > currentX1 && X2 < bestX2) {
      nondomIndizes.push_back(points2[tokeep].first);
      tokeep = i;
      bestX2 = X2;
      currentX1 = X1;
    }
  }
  nondomIndizes.push_back(points2[tokeep].first);


  // std::vector<bool> nondominated(n);
  // for (int i = 0; i < n; ++i) {
  //   nondominated[i] = true;
  // }

  // // FIXME: inefficient
  // for (int i = 0; i < n; ++i) {
  //   for (int j = 0; j < n; ++j) {
  //     if (i == j) {
  //       continue;
  //     }
  //     std::vector<double> p1 = points[i];
  //     std::vector<double> p2 = points[j];
  //     if ((p2[0] < p1[0] && p2[1] <= p1[1]) || (p2[0] <= p1[0] && p2[1] < p1[1])) { //|| ((j < i) && (p1[0] == p2[0]) && (p1[1] == p2[1]))) {
  //       nondominated[i] = false;
  //       break; // break nested loop
  //     }
  //   }
  // }

  // std::vector<int> nondomIndizes;
  // for (int i = 0; i < n; ++i) {
  //   if (nondominated[i]) {
  //     nondomIndizes.push_back(i);
  //   }
  // }
  return nondomIndizes;
}


// @deprecated (basically an alias for getRandomMST)
Graph getMST(Graph * mst) {
  Graph g2 = mst->getRandomMSTKruskal();
  return g2;
}

Graph getRandomMST(Graph* mst) {
  Graph g2 = mst->getRandomMSTKruskal();
  return g2;
}

Graph getMSTByWeightedSumScalarization(Graph* g, std::vector<double> weights) {
  Graph mst = g->getMSTKruskal(weights);
  return mst;
}

Graph getMSTBySubforestMutationR(Graph* g, Graph* mst, int maxDrop, std::vector<double> weights) {
  Graph mstnew = g->getMSTBySubforestMutation(*mst, maxDrop, weights);
  return mstnew;
}

Graph getMSTBySubgraphMutationR(Graph* g, Graph* mst, int maxSelect, std::vector<double> weights) {
  Graph mstnew = g->getMSTBySubgraphMutation(*mst, maxSelect, weights);
  return mstnew;
}

Graph getMSTByInsertionFirstSubgraphMutationR(Graph* g, Graph* mst, int maxAdd, std::vector<double> weights, int maxDegree) {
  Graph mstnew = g->getMSTByInsertionFirstSubgraphMutation(*mst, maxAdd, weights, maxDegree);
  return mstnew;
}

Graph getMSTByUnionCrossoverR(Graph* g, Graph* mst1, Graph* mst2, std::vector<double> weights, int maxDegree) {
  Graph mstnew = g->getMSTByUnionCrossover(*mst1, *mst2, weights, maxDegree);
  return mstnew;
}

Graph getMSTByEdgeExchangeR(Graph* g, Graph* mst, int repls, bool dropLargest) {
  Graph mstnew = g->getMSTByEdgeExchange(*mst, repls, dropLargest);
  return mstnew;
}

NumericVector getSumOfEdgeWeightsR(Graph* g) {
  std::vector<double> sums = g->getSumOfEdgeWeights();
  int n = g->getW();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = sums[i];
  }
  return out;
}

NumericVector getMaxWeight(Graph* g) {
  NumericVector mw(g->getW());

  std::vector<Edge2> edges = g->getEdges();
  // for each weight
  for (Edge2 edge: edges) {
    for (int i = 0; i < g->getW(); ++i) {
      if (edge.second[i] > mw[i]) {
        mw[i] = edge.second[i];
      }
    }
  }

  return mw;
}

int getNumberOfCommonEdges(Graph* g1, Graph* g2) {
  Graph intersection = Graph::getIntersectionGraph(*g1, *g2);
  return intersection.getE();
}

int getNumberOfCommonComponents(Graph* g1, Graph* g2) {
  Graph intersection = Graph::getIntersectionGraph(*g1, *g2);
  std::vector<std::vector<int>> components = intersection.getConnectedComponents();
  return components.size();
}

int getSizeOfLargestCommonComponent(Graph* g1, Graph* g2) {
  Graph intersection = Graph::getIntersectionGraph(*g1, *g2);
  std::vector<std::vector<int>> components = intersection.getConnectedComponents();
  int largest = components[0].size();
  for (int i = 1; i < components.size(); ++i) {
    if (components[i].size() > largest) {
      largest = components[i].size();
    }
  }
  return largest;
}

bool setEdgeProbabilitiesR(Graph* g, NumericVector probs) {
  std::vector<double> probs2(probs.begin(), probs.end());
  g->setEdgeProbabilities(probs2);
  return true;
}

NumericVector getEdgeProbabilitiesR(Graph* g) {
  NumericVector out(g->getE());
  std::vector<double> probs = g->getEdgeProbabilities();
  for (int i = 0; i < probs.size(); ++i) {
    out[i] = probs[i];
  }
  return out;
}

NumericMatrix getWeightsAsMatrix(Graph *g) {
  unsigned int E = g->getE();
  unsigned int W = g->getW();
  NumericMatrix weightMatrix(W, E);

  std::vector<Edge2> edges = g->getEdges();
  for (int i = 0; i < E; ++i) {
    for (int j = 0; j < W; ++j) {
      weightMatrix(j, i) = edges[i].second[j];
    }
  }
  return weightMatrix;
}

NumericMatrix toEdgeList(Graph *g) {
  unsigned int E = g->getE();
  NumericMatrix edgeMatrix(2, E);

  std::vector<Edge2> edges = g->getEdges();
  for (int i = 0; i < E; ++i) {
    edgeMatrix(0, i) = edges[i].first.first;
    edgeMatrix(1, i) = edges[i].first.second;
  }
  return edgeMatrix;
}

NumericVector toBitstringR(Graph *g, Graph *subgraph) {
  std::vector<int> bitstring = g->toBitstring(*subgraph);
  int m = g->getE();
  NumericVector bitstringR(m);
  for (int i = 0; i < m; ++i) {
    bitstringR[i] = bitstring[i];
  }
  return bitstringR;
}

List doMCPrim(Graph *g) {
  int V = g->getV();

  std::vector<Graph> trees = g->doMCPrim();
  unsigned int ntrees = trees.size();
  //FIXME: generalize
  unsigned int W = g->getW();

  NumericMatrix costMatrix(W, ntrees);
  //FIXME: this is absolutely ultra-ugly!!!
  NumericMatrix edgeMatrix(2 * ntrees, V - 1);

  for (int i = 0; i < ntrees; ++i) {
    // Fixme
    Graph tree = trees[i];
    std::vector<double> costs = trees[i].getSumOfEdgeWeights();
    for (int j = 0; j < W; ++j) {
      costMatrix(j, i) = costs[j];
    }

    std::vector<Edge2> treeEdges = tree.getEdges();
    for (int j = 0; j < treeEdges.size(); ++j) {
      Edge2 edge = treeEdges[j];
      // store each two rows for edges
      edgeMatrix(2 * i, j) = edge.first.first;
      edgeMatrix(2 * i + 1, j) = edge.first.second;
    }
  }

  List output;
  output["weights"] = costMatrix;
  output["edges"] = edgeMatrix;
  return output;
}


RCPP_EXPOSED_CLASS(Graph)
RCPP_EXPOSED_CLASS(GraphImporter)
RCPP_EXPOSED_CLASS(RepresentationConverter)

RCPP_MODULE(graph_module) {
  using namespace Rcpp;
  class_<Graph>("Graph")
    .constructor<int, int>()
    .method("getV", &Graph::getV)
    .method("getE", &Graph::getE)
    .method("getW", &Graph::getW)
    .method("getDegree", &Graph::getDegree)
    .method("getMaximumDegree", &Graph::getMaximumDegree)
    .method("getNumberOfLeafs", &Graph::getNumberOfLeafs)
    .method("getDiameter", &Graph::getDiameter)
    .method("saveVectorOfEdges", &Graph::saveVectorOfEdges)
    .method("addEdge", &Graph::addEdge)
    .method("isSpanningTree", &Graph::isSpanningTree)
    //.method("getMSTKruskal", &Graph::getMSTKruskal)
    //.method("importFromGrapheratorFile", &importFromGrapheratorFileR)
    .method("getMST", &getMST)
    .method("getRandomMST", &getRandomMST)
    .method("getMSTByWeightedSumScalarization", &getMSTByWeightedSumScalarization)
    .method("getSumOfEdgeWeights", &getSumOfEdgeWeightsR)
    .method("getMaxWeight", &getMaxWeight)
    .method("getNumberOfCommonEdges", &getNumberOfCommonEdges)
    .method("getNumberOfCommonComponents", &getNumberOfCommonComponents)
    .method("getSizeOfLargestCommonComponent", &getSizeOfLargestCommonComponent)
    .method("getMSTBySubforestMutation", &getMSTBySubforestMutationR)
    .method("getMSTBySubgraphMutation", &getMSTBySubgraphMutationR)
    .method("getMSTByInsertionFirstSubgraphMutation", &getMSTByInsertionFirstSubgraphMutationR)
    .method("getMSTByUnionCrossover", &getMSTByUnionCrossoverR)
    .method("getMSTByEdgeExchange", &getMSTByEdgeExchangeR)
    .method("getWeightsAsMatrix", &getWeightsAsMatrix)
    .method("toEdgeList", &toEdgeList)
    .method("toBitstring", &toBitstringR)
    .method("getEdgeProbabilities", &getEdgeProbabilitiesR)
    .method("setEdgeProbabilities", &setEdgeProbabilitiesR)
    .method("doMCPrim", &doMCPrim)
  ;
  class_<GraphImporter>("GraphImporter")
    .constructor()
    .method("importFromGrapheratorFile", &GraphImporter::importFromGrapheratorFileCPP)
  ;
  class_<RepresentationConverter>("RepresentationConverter")
    .constructor()
    .method("prueferCodeToGraph", &RepresentationConverter::prueferCodeToGraph)
    .method("edgeListToGraph", &RepresentationConverter::edgeListToGraph)
  ;
}
