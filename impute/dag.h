/* Copyright notice.... */
#ifndef DAG_H
#define DAG_H

#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"
#include "impute/haplotypepair.h"

class Dag
{
public:
  /**
   * Returns the number of edges at the specified level of the DAG.
   *
   * @param level a level of the DAG
   */
  virtual int nEdges(int level) const = 0;

  /**
   * Returns the number of parent nodes at the specified level of the DAG.
   *
   * @param level a level of the DAG
   */
  virtual int nParentNodes(int level) const = 0;

  /**
   * Returns the number of child nodes at the specified level of the DAG.
   *
   * @param level a level of the DAG
   */
  virtual int nChildNodes(int level) const = 0;

  /**
   * Returns the index of the specified parent node in the DAG.
   *
   * @param level a level of the DAG.
   * @param edge the index of an edge at the specified level of the DAG
   */
  virtual int parentNode(int level, int edge) const = 0;

  /**
   * Returns the index of the specified child node in the DAG.
   *
   * @param level a level of the DAG
   * @param edge the index of an edge at the specified level of the DAG
   */
  virtual int childNode(int level, int edge) const = 0;

  /**
   * Returns the symbol labeling the specified edge of the DAG.
   *
   * @param level a level of the DAG
   * @param edge the index of an edge at the specified level of the DAG
   */
  virtual int symbol(int level, int edge) const = 0;

  /**
   * Returns the sum of the weights of the sequences that pass
   * through the specified edge of the DAG.
   *
   * @param level a level of the DAG
   * @param edge the index of an edge at the specified level of the DAG
   */
  virtual double edgeWeight(int level, int edge) const = 0;

  /**
   * Returns the sum of the weights of the sequences that pass
   * through the specified node of the DAG.
   *
   * @param level a level of the DAG
   * @param parentNode the index of a parent node at the specified level
   * of the DAG
   */
  virtual double parentWeight(int level, int parentNode) const = 0;

  /**
   * Returns the ratio of the sum of the weights of the sequences that pass
   * through the specified edge of the DAG and
   * the sum of the weights of the sequences that pass through the parent
   * node of the specified edge of the DAG.
   *
   * @param level a level of the DAG
   */
  virtual double condEdgeProb(int level, int edge) const = 0;

  /**
   * Returns the ratio of the sum of the weights of the sequences that pass
   * through the specified edge of the DAG and the sum of the weights of all
   * sequences.
   *
   * @param level a level of the DAG
   * @param edge the index of an edge at the specified level of the DAG
   */
  virtual double edgeProb(int level, int edge) const = 0;

  /**
   * Returns the ratio of the sum of the weights of the sequences that pass
   * through the specified parent node of the DAG and the sum of the weights
   * of all sequences.
   *
   * @param level a level of the DAG
   * @param parentNode the index of a parent node at the specified level
   * of the DAG
   */
  virtual double parentProb(int level, int parentNode) const = 0;

  /**
   * Returns the number of markers.
   */
  virtual int nLevels() const = 0;

  /**
   * Returns the markers represented by this DAG.
   */
  virtual Markers markers() const = 0;

  /**
   * Returns the number of nodes in the DAG.
   */
  virtual int nNodes() const = 0;

  /**
   * Returns the number of edges in the DAG.
   */
  virtual int nEdges() const = 0;

  /**
   * Returns the maximum number of parent nodes at any level of the DAG.
   */
  virtual int maxNodes() const = 0;

  /**
   * Returns the maximum number of edges at any level of the DAG.
   */
  virtual int maxEdges() const = 0;

  /**
   * Returns the number of outgoing edges for the specified node of the DAG.
   *
   * @param level a level of the DAG
   * @param parentNode the index of a parent node at the specified
   * level of the DAG
   */
  virtual int nOutEdges(int level, int parentNode) const = 0;

  /**
   * Returns the index of the specified edge in the DAG.
   *
   * @param level a level of the DAG
   * @param parentNode the index of a parent node at the specified
   * level of the DAG
   * @param outEdge the index of an outgoing edge of the specified
   * parent node
   */
  virtual int outEdge(int level, int parentNode, int outEdge) const = 0;

  /**
   * Returns the index of the specified edge at the specified level of the
   * DAG or {@code -1} if no such edge exists.
   *
   * @param level a level of the DAG
   * @param parentNode the index of a parent node at the specified
   * level of the DAG
   * @param symbol a symbol labeling an outgoing edge of the specified
   * parent node of the DAG
   */
  virtual int outEdgeBySymbol(int level, int parentNode, int symbol) const = 0;

  /**
   * Returns the number of ingoing edges for the specified node of the DAG.
   *
   * @param level a level of the DAG
   * @param childNode the index of a child node at the specified
   * level of the DAG
   */
  virtual int nInEdges(int level, int childNode) const = 0;

  /**
   * Returns the index of the specified edge in the DAG.
   *
   * @param level a level of the DAG
   * @param childNode the index of a child node at the specified
   * level of the DAG
   * @param inEdge the index of an ingoing edge of the specified
   * child node in the DAG
   */
  virtual int inEdge(int level, int childNode, int inEdge) const = 0;

  /**
   * Returns {@code true} if the child node of the specified parent
   * edge equals the parent node of the specified child edge and
   * returns {@code false} otherwise.
   *
   * @param parentLevel a level of the DAG
   * @param parentEdge the index of an edge at the specified level
   * of the DAG
   * @param childEdge the index of an edge at level {@code (parentLevel + 1)}
   * of the DAG
   */
  virtual bool isChildOf(int parentLevel, int parentEdge, int childEdge) const = 0;

  /**
   * Returns a QList of length {@code this.nMarkers()} whose {@code j}-th
   * element is the distance from the root node to
   * the child node at level {@code j} of the DAG.
   * The distance from parent node to child node at level {@code lev}
   * equals {@code -Math.log10(P)} where {@code P} is the weighted conditional
   * edge probability at level {@code lev}, when each edge {@code e} is
   * weighted by {@code this.counts(lev, e)}.
   */
  virtual QList<double> posArray() const = 0;
};

class LinkageEquilibriumDag : public Dag
{
public:
  /**
   * Constructs a new {@code LinkageEquilibriumDag} instance that represents
   * markers in linkage equilibrium, with one level per marker,
   * one parent node per level, one edge per allele at each level,
   * and edge count equal to the estimated allele frequency.
   * @param gl the genotype emission probabilities which determine
   * the estimated allele frequencies
   * @param minFreq the minimum allele frequency that will be used
   */
  LinkageEquilibriumDag(const SplicedGL &gl, double minFreq);

  int nEdges(int level) const { return _alleleFreq[level].length(); }
  int nParentNodes(int level) const
  {
    checkLevel(level);
    return 1;
  }

  int nChildNodes(int level) const
  {
    checkLevel(level);
    return 1;
  }

  int parentNode(int level, int edge) const
  {
    checkEdge(level, edge);
    return 0;
  }

  int childNode(int level, int edge) const
  {
    checkEdge(level, edge);
    return 0;
  }

  int symbol(int level, int edge) const
  {
    checkEdge(level, edge);
    return edge;
  }

  double edgeWeight(int level, int edge) const { return _alleleFreq[level][edge]; }
  double parentWeight(int level, int parentNode) const
  {
    checkParentNode(level, parentNode);
    return 1.0;
  }

  double condEdgeProb(int level, int edge) const { return _alleleFreq[level][edge]; }
  double edgeProb(int level, int edge) const { return _alleleFreq[level][edge]; }
  double parentProb(int level, int node) const
  {
    checkParentNode(level, node);
    return 1.0;
  }

  int nLevels() const { return _alleleFreq.length(); }
  Markers markers() const { return _markers; }
  int nNodes() const { return (_alleleFreq.length() + 1); }
  int nEdges() const { return _sumAlleles; }
  int maxNodes() const { return 1; }
  int maxEdges() const { return _maxAlleles; }
  int nOutEdges(int level, int parentNode) const { return _alleleFreq[level].length(); }
  int outEdge(int level, int parentNode, int outEdge) const;

  int outEdgeBySymbol(int level, int parentNode, int symbol) const { return symbol; }
  int nInEdges(int level, int childNode) const;

  int inEdge(int level, int childNode, int inEdge) const;

  bool isChildOf(int parentLevel, int parentEdge, int childEdge) const;

  QList<double> posArray() const;

private:
  QList<double> alleleFrequencies(const SplicedGL &gl, int marker, double minFreq);
  void divideEntriesBySum(QList<double> &fa);
  void enforceMinFrequency(QList<double> &alleleFreq, double minFreq);
  void checkLevel(int level) const;
  void checkEdge(int level, int edge) const;
  void checkParentNode(int level, int node) const;

  Markers _markers;
  QList<QList<double> > _alleleFreq;
  int _maxAlleles;
  int _sumAlleles;
};

#endif
