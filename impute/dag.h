/* Copyright notice.... */
#ifndef DAG_H
#define DAG_H

#include <QVector>

#include "impute/haplotypepair.h"
#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"

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
  virtual float edgeWeight(int level, int edge) const = 0;

  /**
   * Returns the sum of the weights of the sequences that pass
   * through the specified node of the DAG.
   *
   * @param level a level of the DAG
   * @param parentNode the index of a parent node at the specified level
   * of the DAG
   */
  virtual float parentWeight(int level, int parentNode) const = 0;

  /**
   * Returns the ratio of the sum of the weights of the sequences that pass
   * through the specified edge of the DAG and
   * the sum of the weights of the sequences that pass through the parent
   * node of the specified edge of the DAG.
   *
   * @param level a level of the DAG
   */
  virtual float condEdgeProb(int level, int edge) const = 0;

  /**
   * Returns the ratio of the sum of the weights of the sequences that pass
   * through the specified edge of the DAG and the sum of the weights of all
   * sequences.
   *
   * @param level a level of the DAG
   * @param edge the index of an edge at the specified level of the DAG
   */
  virtual float edgeProb(int level, int edge) const = 0;

  /**
   * Returns the ratio of the sum of the weights of the sequences that pass
   * through the specified parent node of the DAG and the sum of the weights
   * of all sequences.
   *
   * @param level a level of the DAG
   * @param parentNode the index of a parent node at the specified level
   * of the DAG
   */
  virtual float parentProb(int level, int parentNode) const = 0;

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
  virtual long nNodes() const = 0;

  /**
   * Returns the number of edges in the DAG.
   */
  virtual long nEdges() const = 0;

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
  LinkageEquilibriumDag(const SplicedGL &gl, float minFreq);

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

  float edgeWeight(int level, int edge) const { return _alleleFreq[level][edge]; }
  float parentWeight(int level, int parentNode) const
  {
    checkParentNode(level, parentNode);
    return 1.0;
  }

  float condEdgeProb(int level, int edge) const { return _alleleFreq[level][edge]; }
  float edgeProb(int level, int edge) const { return _alleleFreq[level][edge]; }
  float parentProb(int level, int node) const
  {
    checkParentNode(level, node);
    return 1.0;
  }

  int nLevels() const { return _alleleFreq.length(); }
  Markers markers() const { return _markers; }
  long nNodes() const { return (_alleleFreq.length() + 1); }
  long nEdges() const { return _sumAlleles; }
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
  QList<float> alleleFrequencies(const SplicedGL &gl, int marker, float minFreq);
  void divideEntriesBySum(QList<float> &fa);
  void enforceMinFrequency(QList<float> &alleleFreq, float minFreq);
  void checkLevel(int level) const;
  void checkEdge(int level, int edge) const;
  void checkParentNode(int level, int node) const;

  Markers _markers;
  QList<QList<float> > _alleleFreq;
  int _maxAlleles;
  int _sumAlleles;
};

/**
 * Class {@code DagLevel} represents a level (within class {@code
 * ImmutableDag}) of a leveled directed acyclic graph (DAG). This
 * implementation can contain up to 65535 edges. ("Medium capacity".)
 *
 * "All instances of {@code DagLevel} are required to be immutable."
 */
class DagLevel
{
public:
  /**
   * Constructs a new {@code DagLevel} instance from the specified
   * data.
   *
   * @param parentNodes an array mapping edge index to parent node index
   * @param childNodes an array mapping edge index to child node index
   * @param symbols an array mapping edge index to the symbol labeling the
   * edge
   * @param counts an array mapping edge index to edge count
   */
  DagLevel(QVector<quint16> parentNodes, QVector<quint16> childNodes,
           QVector<quint16> symbols, QVector<float> counts);

  /**
   * Returns the number of edges at this level of the DAG.
   */
  int nEdges() const { return _condEdgeProbs.length(); }
  /**
   * Returns the number of parent nodes at this level of the DAG.
   */
  int nParentNodes() const { return _parentIndices.length() - 1; }
  /**
   * Returns the number of child nodes at this level of the DAG.
   */
  int nChildNodes() const { return _childIndices.length() - 1; }
  /**
   * Returns the index of the parent node of the specified edge
   * at this level of the DAG.
   *
   * @param edge an edge index
   */
  int parentNode(int edge) const { return _parentNodes[edge]; }
  /**
   * Returns the index of the child node of the specified edge
   * at this level of the DAG.
   *
   * @param edge an edge index.
   */
  int childNode(int edge) const { return _childNodes[edge]; }
  /**
   * Returns the symbol labeling the specified edge at this level
   * of the DAG.
   *
   * @param edge an edge index
   */
  int symbol(int edge) const { return _symbols[edge]; }
  /**
   * Returns the sum of weights for the sequences that pass
   * through the specified edge at this level of the DAG.
   *
   * @param edge an edge index
   */
  float edgeWeight(int edge) const { return _edgeCounts[edge]; }
  /**
   * Returns the sum of weights for the sequences that pass
   * through the specified node at this level of the DAG.
   *
   * @param parentNode a parent node index
   */
  float parentWeight(int parentNode) const { return _parentCounts[parentNode]; }
  /**
   * Returns the conditional edge probability, which is defined to be
   * the ratio of the sum of the weights of the sequences that pass
   * through the specified edge at this level of the DAG and
   * the sum of the weights of the sequences that pass through the parent
   * node of the specified edge.
   *
   * @param edge an edge index
   */
  float condEdgeProb(int edge) const { return _condEdgeProbs[edge]; }
  /**
   * Returns the edge probability, which is defined to be the ratio of the
   * sum of the weights of the sequences that pass through the specified
   * edge at this level of the DAG and the sum of the weights of the
   * sequences that pass through any edge at this level of the DAG.
   *
   * @param edge an edge index
   */
  float edgeProb(int edge) const { return (_edgeCounts[edge] / _count); }
  /**
   * Returns the parent node probability, which is defined to be the
   * ratio of the sum of the weights of the sequences that pass through
   * the specified parent node at this level of the DAG and the sum of
   * the weights of the sequences that pass through any parent node at this
   * level of the DAG.
   *
   * @param parentNode a parent node index
   */
  float parentProb(int node) const { return (_parentCounts[node] / _count); }
  /**
   * Returns the number of outgoing edges of the specified parent node
   * at this level of the DAG.
   *
   * @param parentNode a parent node index
   */
  int nOutEdges(int parentNode) const
  {
    return (_parentIndices[parentNode + 1] - _parentIndices[parentNode]);
  }

  /**
   * Returns the index of the specified edge at this level of the DAG.
   *
   * @param parentNode a parent node index
   * @param outEdge the index of the outgoing edge of the specified
   * parent node
   */
  int outEdge(int parentNode, int outEdgeIndex) const;

  /**
   * Returns the index of the specified edge at this level of the
   * DAG or {@code -1} if no such edge exists.
   *
   * @param parentNode a parent node index
   * @param symbol a symbol labeling an outgoing edge of the specified
   * parent node
   */
  int outEdgeBySymbol(int parentNode, int symbol) const;

  /**
   * Returns the number of ingoing edges for the specified child node
   * at this level of the DAG.
   *
   * @param childNode a child node index
   */
  int nInEdges(int childNode) const
  {
    return (_childIndices[childNode + 1] - _childIndices[childNode]);
  }

  /**
   * Returns the index of the specified edge at this level of the DAG.
   *
   * @param childNode index of the child node
   * @param inEdge index of an ingoing edge of the specified child node
   */
  int inEdge(int childNode, int inEdgeIndex) const;

private:
  void obtainParentCounts(QVector<quint16> parentNodes, QVector<float> counts, int nNodes);

  /*
   * The k-th edge parent node index is stored in {@code _parentNodes[k]}.
   * The k-th edge child node index is stored in {@code _childNodes[k]}.
   * The k-th edge symbol is stored in {@code _symbols[k]}.
   * The k-th edge count is stored in {@code _edgeCounts[k]}.
   * The k-th edge conditional edge probability is stored in
   * {@code _condEdgeProbs[k]}, and is defined to be the
   * k-th edge count divided by the k-th edge's parent node count.
   * The k-th node count is stored in {@code _nodeCounts[k]}.
   *
   * The outgoing edges indices of the k-th parent node are stored in consecutive
   * entries of {@code _parents} beginning with
   * {@code _parentIndices[k]} (inclusive) and ending with
   * {@code _parentIndices[k+1]} (exclusive).
   *
   * The ingoing edges indices of the k-th child node are stored in consecutive
   * entries of {@code _children} beginning with
   * {@code _childIndices[k]} (inclusive) and ending with
   * {@code _childIndices[k+1]} (exclusive).
   */
  float _count;
  QVector<quint16> _parentNodes;
  QVector<quint16> _childNodes;
  QVector<quint16> _parentIndices;
  QVector<quint16> _parents;
  QVector<quint16> _childIndices;
  QVector<quint16> _children;
  QVector<quint16> _symbols;
  QVector<float> _edgeCounts;
  QVector<float> _condEdgeProbs;
  QVector<float> _parentCounts;
};

/**
 * Class {@code ImmutableDag} represents a leveled Directed Acyclic Graph
 * (DAG).
 *
 * "Instances of class {@code ImmutableDag} are immutable."
 *
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
class ImmutableDag : public Dag
{
public:
  /**
   * Constructs a new {@code ImmutableDag} instance from the specified data.
   * @param hapPairs the sequence data
   * @param weights an array whose {@code j}-th element is the
   * weight for the {@code j}-th haplotype
   * @param scale a parameter that multiplicatively scales the node
   * similarity threshold
   * @param nInitLevels the number of initial levels to read
   *
   * (Internally constructs a new {@code MergeableDagFactory}
   *  instance from the above specified data, then uses that to
   *  finish constructing an {@code ImmutableDag} instance by
   *  obtaining the following data from hapPairs and from the
   *  factory instance:
   *  @param markers the markers
   *  @param levels the levels of the leveled DAG)
   */
  ImmutableDag(const HapPairs &hapPairs, const QVector<float> &weights, float scale,
               int nInitLevels);

  int nEdges(int level) const { return _dagLevels[level].nEdges(); }
  int nParentNodes(int level) const { return _dagLevels[level].nParentNodes(); }
  int nChildNodes(int level) const { return _dagLevels[level].nChildNodes(); }
  int parentNode(int level, int edge) const { return _dagLevels[level].parentNode(edge); }
  int childNode(int level, int edge) const { return _dagLevels[level].childNode(edge); }
  int symbol(int level, int edge) const { return _dagLevels[level].symbol(edge); }
  float edgeWeight(int level, int edge) const { return _dagLevels[level].edgeWeight(edge); }
  float parentWeight(int level, int parentNode) const
  {
    return _dagLevels[level].parentWeight(parentNode);
  }

  float condEdgeProb(int level, int edge) const { return _dagLevels[level].condEdgeProb(edge); }
  float edgeProb(int level, int edge) const { return _dagLevels[level].edgeProb(edge); }
  float parentProb(int level, int node) const { return _dagLevels[level].parentProb(node); }
  int nLevels() const { return _dagLevels.length(); }
  Markers markers() const { return _markers; }
  long nNodes() const { return _nNodes; }
  long nEdges() const { return _nEdges; }
  int maxNodes() const { return _maxNodes; }
  int maxEdges() const { return _maxEdges; }
  int nOutEdges(int level, int parentNode) const { return _dagLevels[level].nOutEdges(parentNode); }
  int outEdge(int level, int parentNode, int outEdge) const
  {
    return _dagLevels[level].outEdge(parentNode, outEdge);
  }

  int outEdgeBySymbol(int level, int parentNode, int symbol) const
  {
    return _dagLevels[level].outEdgeBySymbol(parentNode, symbol);
  }

  int nInEdges(int level, int childNode) const { return _dagLevels[level].nInEdges(childNode); }
  int inEdge(int level, int childNode, int inEdge) const
  {
    return _dagLevels[level].inEdge(childNode, inEdge);
  }

  bool isChildOf(int parentLevel, int parentEdge, int childEdge) const
  {
    int nodeA = _dagLevels[parentLevel + 1].parentNode(childEdge);
    int nodeB = _dagLevels[parentLevel].childNode(parentEdge);
    return nodeA == nodeB;
  }

  QList<double> posArray() const { return _posArray; }

private:
  Markers _markers;

  long _nNodes;   // total number of nodes
  long _nEdges;   // total number of edges
  int _maxNodes;  // maximum number of nodes on any level
  int _maxEdges;  // maximum number of edges on any level

  QList<DagLevel> _dagLevels;
  QList<double> _posArray;
};

#endif
