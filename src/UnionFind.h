#include <iostream>
#include <cstdlib>
#include <vector>
#include <assert.h>

#ifndef UNION_FIND_DATASTRUCTURE
#define UNION_FIND_DATASTRUCTURE
class UnionFind {
public:
  //FIXME: code optimization?
  //FIXME: add tests
  /**
   * General constructor.
   *
   * @param[in] n Number of elements.
   * @param[in] pc Activate path compression?
   * @return Object of type UnionFind.
   */
  UnionFind(unsigned int n, bool pc = true) {
    assert(n >= 2);

    this->n = n;
    this->nsets = n;
    this->pc = pc;
    // vector is zero based, i.e., we add a dummy element here
    // for simplification reasons
    this->root.reserve(n + 1);
    this->size.reserve(n + 1);
    for (unsigned int i = 1; i <= n; ++i) {
      this->root[i] = i;
      this->size[i] = 1;
    }
  }

  /**
   * Constructor based on non-isolated sets.
   *
   * @param[in] n Number of elements.
   * @param[in] sets Sets.
   * @param[in] pc Activate path compression?
   * @return Object of type UnionFind.
   */
  UnionFind(unsigned int n, std::vector<std::vector<int>> sets, bool pc = true) {
    assert(n >= 2);

    this->n = n;
    this->nsets = sets.size();
    this->pc = pc;

    this->root.reserve(n + 1);
    this->size.reserve(n + 1);

    // now go through all sets
    for (auto set: sets) {
      int setRoot;
      int setSize = set.size();
      // now go through set
      for (int i = 0; i < setSize; ++i) {
        // get element
        int element = set[i];
        // first element in set is always the root
        if (i == 0) {
          setRoot = element;
        }
        this->root[element] = setRoot;
        this->size[element] = setSize;
      }
    }
  }

  /**
   * Get root element / representative.
   *
   * @param[in] i Element.
   * @return Root element.
   */
  int getRoot(unsigned int i) {
    assert(i >= 1 & i <= this->n);

    // propagate parent
    unsigned int j = i;
    while (j != root[j]) {
      j = root[j];
    }

    // iterate again. However, set pointer to root directly.
    if (this->pc) {
      unsigned int k = i;
      while (k != root[k]) {
        unsigned int tmp = k;
        k = root[k];
        root[tmp] = j;
      }
    }

    return j;
  }

  /**
   * Check if two elements are in the same set.
   *
   * @param i, j Elements.
   * @return Boolean indicating whether i and j are in the same set.
   */
  bool find(unsigned int i, unsigned int j) {
    return getRoot(i) == getRoot(j);
  }

  /**
   * Set-union operation.
   *
   * @param i, j Set elements.
   */
  void unite(unsigned int i, unsigned int j) {
    // sanity checks
    assert(i >= 1 & i <= this->n);
    assert(j >= 1 & j <= this->n);

    unsigned int root_i = getRoot(i);
    unsigned int root_j = getRoot(j);

    // only if not already in same set
    if (root_i != root_j) {
      if (size[i] < size[j]) {
        root[root_j] = i;
        size[i] += size[root_j];
      } else {
        root[root_i] = j;
        size[j] += size[root_i];
      }
      this->nsets -= 1;
    }
  }

  /// Printer for debugging.
  void print() {
    std::cout << "UnionFind DS: n = " << this->n << std::endl;
    std::vector<bool> output(this->n);
    for (unsigned int i = 1; i <= this->n; ++i) {
      if (output[i]) {
        continue;
      }
      for (unsigned int j = 1; j <= this->n; ++j) {
        if (getRoot(j) == i) {
          std::cout << j << ", ";
          output[j] = true;
        }
      }
      std::cout << std::endl;
      output[i] = true;
    }
  }

private:
  /// number of elements initially added
  unsigned int n;
  /// number of sets
  unsigned int nsets;
  /// path compression flag
  bool pc;
  /// pointers to root elements
  std::vector<unsigned int> root;
  /// sizes of the set of each element
  std::vector<unsigned int> size;
};
#endif
