#include "impute/dag.h"

#include <float.h>
#include <math.h>

#include <QSet>
#include <QLinkedList>
#include <QMutableLinkedListIterator>

#define UINT16_MAX 65535

#define MAX_PROP_UNMERGED 0.01
#define MIN_DEPTH 10
#define MAX_THRESHOLD_RATIO 1.4

LinkageEquilibriumDag::LinkageEquilibriumDag(const SplicedGL &gl, float minFreq)
{
  Q_ASSERT_X(minFreq > 0.0 || minFreq < 0.5,
             "LinkageEquilibriumDag::LinkageEquilibriumDag",
             "bad value for minFreq");

  int nMarkers = gl.nMarkers();
  int localMaxAlleles = 0;
  _markers = gl.markers();

  for (int marker = 0; marker < nMarkers; ++marker) {
    _alleleFreq.append(alleleFrequencies(gl, marker, minFreq));
    if (_alleleFreq[marker].length() > localMaxAlleles) {
      localMaxAlleles = _alleleFreq[marker].length();
    }
  }
  _maxAlleles = localMaxAlleles;
  _sumAlleles = gl.markers().sumAlleles();
}

QList<float> LinkageEquilibriumDag::alleleFrequencies(const SplicedGL &gl,
                                                       int marker,
                                                       float minFreq)
{
  int nSamples = gl.nSamples();
  int nAlleles = gl.marker(marker).nAlleles();
  QList<float> alleleFreq;
  QList<float> scaledFreq;
  for (int a1 = 0; a1 < nAlleles; a1++) {
    alleleFreq.append(0.0);
    scaledFreq.append(0.0);
  }
  for (int sample = 0; sample < nSamples; ++sample) {
    for (int a1 = 0; a1 < nAlleles; ++a1) {
      for (int a2 = 0; a2 < nAlleles; ++a2) {
        float likelihood = gl.gl(marker, sample, a1, a2);
        scaledFreq[a1] += likelihood;
        scaledFreq[a2] += likelihood;
      }
    }
    divideEntriesBySum(scaledFreq);
    for (int j = 0; j < scaledFreq.length(); ++j) {
      alleleFreq[j] += scaledFreq[j];
      scaledFreq[j] = 0.0f;
    }
  }
  divideEntriesBySum(alleleFreq);
  enforceMinFrequency(alleleFreq, minFreq);
  return alleleFreq;
}

void LinkageEquilibriumDag::divideEntriesBySum(QList<float> &fa)
{
  float sum = 0.0f;
  for (int j = 0; j < fa.length(); ++j)
    sum += fa[j];

  for (int j = 0; j < fa.length(); ++j)
    fa[j] /= sum;
}

void LinkageEquilibriumDag::enforceMinFrequency(QList<float> &alleleFreq, float minFreq)
{
  bool changedFreq = false;
  for (int j = 0; j < alleleFreq.length(); ++j) {
    if (alleleFreq[j] < minFreq) {
      alleleFreq[j] = minFreq;
      changedFreq = true;
    }
  }
  if (changedFreq)
    divideEntriesBySum(alleleFreq);
}

int LinkageEquilibriumDag::outEdge(int level, int parentNode, int outEdge) const
{
  checkParentNode(level, parentNode);
  checkEdge(level, outEdge);
  return outEdge;
}

int LinkageEquilibriumDag::nInEdges(int level, int childNode) const
{
  checkLevel(level);
  return _alleleFreq[level].length();
}

int LinkageEquilibriumDag::inEdge(int level, int childNode, int inEdge) const
{
  checkEdge(level, inEdge);
  Q_ASSERT_X(childNode == 0, "LinkageEquilibriumDag::inEdge", "childNode != 0");
  return inEdge;
}

bool LinkageEquilibriumDag::isChildOf(int parentLevel, int parentEdge, int childEdge) const
{
  checkEdge(parentLevel, parentEdge);
  checkEdge(parentLevel + 1, childEdge);
  return true;
}

QList<double> LinkageEquilibriumDag::posArray() const
{
  QList<double> pos;
  for (int j = 0; j < _alleleFreq.length(); ++j) {
    double condEdgeProb = 0.0;
    for (int a = 0; a < _alleleFreq[j].length(); ++a)
      condEdgeProb += _alleleFreq[j][a] * _alleleFreq[j][a];

    if (j == 0)
      pos.append(-log10(condEdgeProb));
    else
      pos.append(pos[j - 1] - log10(condEdgeProb));
  }
  return pos;
}

void LinkageEquilibriumDag::checkLevel(int level) const
{
  Q_ASSERT_X(level >= 0 && level < _alleleFreq.length(),
             "LinkageEquilibriumDag::checkLevel",
             "invalid level");
}

void LinkageEquilibriumDag::checkEdge(int level, int edge) const
{
  Q_ASSERT_X(edge >= 0 && edge < _alleleFreq[level].length(),
             "LinkageEquilibriumDag::checkEdge",
             "invalid edge for this level");
}

void LinkageEquilibriumDag::checkParentNode(int level, int node) const
{
  checkLevel(level);
  Q_ASSERT_X(node == 0, "LinkageEquilibriumDag::checkParentNode", "invalid parent node (!= 0)");
}

/**
 * Iterator for the HapPairs object--except that instead of iterating
 * over the haplotypes, we iterate over the markers.
 */
class HapsMarkerIterator
{
public:
  HapsMarkerIterator(HapPairs hp) : _hp(hp), _pos(-1) {}
  bool hasNext() const { return _pos < (_hp.nMarkers() - 1); }
  void next() { _pos++; }
  int nHaps() const { return _hp.nHaps(); }
  Marker marker() const { return _hp.marker(_pos); }
  int allele(int hap) const { return _hp.allele(_pos, hap); }
private:
  HapPairs _hp;
  int _pos;
};

/**
 * Class {@code MergeableDagLevel} represents a level of a leveled
 * directed acyclic graph (DAG). The class includes a public method
 * for merging parent nodes, and is designed to be used by {@code
 * MergeableDagFactory}.
 *
 * "Instances of class {@code MergebleDagLevel} are not thread-safe."
 */
class MergeableDagLevel
{
public:
  /**
   * Constructs a new {@code MergeableDagLevel} instance from the specified
   * phased genotype data and haplotype weights.  The {@code previous()}
   * method of the constructed instance will return zero.
   * @param data the phased genotype data
   * @param weights an array mapping haplotype indices to non-negative
   * weights
   */
  MergeableDagLevel(const HapsMarkerIterator &data, QVector<float> weights);

  /**
   * Constructs a new {@code MergeableDagLevel} instance with the
   * specified previous {@code MergeableDagLevel} and the
   * specified phased genotype data.  This constructor does not alter
   * any field of the specified {@code prevLevel} object.
   * @param prevLevel (a pointer to) the previous {@code MergeableDagLevel}
   * @param data the phased genotype data
   */
  MergeableDagLevel(MergeableDagLevel *prevLevel, const HapsMarkerIterator &data);

  /**
   * Sets the previous DAG level to {@code null}, and returns
   * the previous DAG level that existed immediately prior to the invocation
   * of this method.
   */
  MergeableDagLevel *takePrevious();

  /**
   * Sets the next level to the specified {@code MergeableDagLevel}.
   * @param nextLevel the next level
   */
  void setNextLevel(MergeableDagLevel *nextLevel);

  /**
   * Returns the previous DAG level or {@code null} if no previous level
   * is stored.
   */
  MergeableDagLevel *previous() const { return _prevLevel; }
  /**
   * Returns the next DAG level or {@code null} if no next level is stored.
   */
  MergeableDagLevel *next() const { return _nextLevel; }
  /**
   * Returns {@code true} if the specified parent node has a
   * sibling  and returns {@code false} otherwise.
   * Two parent nodes are siblings if they are connected by an
   * edge to the same parent node at the previous level of the DAG.
   *
   * @param parentNode a parent node index
   */
  bool hasSibling(int parentNode) const { return _prevLevel->hasSiblingInMe(parentNode); }
  /**
   * Manufactures an immutable {@code DagLevel} corresponding to
   * this object. The parent node, edge, and child node indices
   * in the returned {@code DagLevel} are the ranks of the
   * parent node, edge, and child node indices for this object,
   * with rank 0 corresponding to the smallest index.
   */
  DagLevel toDagLevel();

  /**
   * Merges the two specified parent nodes and assigns the specified
   * {@code retainedNode} index to the merged node.
   *
   * @param retainedNode a parent node which will receive ingoing and
   * outgoing edges of {@code removedNode}
   * @param removedNode a parent node that will be deleted after merging.
   */
  void mergeParentNodes(int retainedNode, int removedNode);

  /**
   * Returns the marker index.
   */
  int index() const { return _levelIndex; }
  /**
   * Returns the number of sequences used to construct the DAG.
   */
  int nHaps() const { return _nHaps; }
  /**
   * Returns the number of alleles.
   */
  int nAlleles() const { return _nAlleles; }
  /**
   * Returns the sum of weights for the sequences that pass
   * through the specified edge or 0 if the edge does not exist.
   *
   * @param edge index of the edge
   */
  float edgeCount(int edge) const { return _counts[edge]; }
  /**
   * Returns the sum of weights for the sequences that pass
   * through the specified parent node or 0 if the parent node
   * does not exist.
   *
   * @param parentNode index of the parent node
   */
  float nodeCount(int parentNode) const;

  /**
   * Returns an array of parent node indices.
   */
  QList<int> parentNodeArray() const;

  /**
   * Returns the parent node of the specified edge or -1 if the edge does
   * not exist.
   *
   * @param edge index of the edge
   */
  int parentNode(int edge) const { return _parentNodes[edge]; }
  /**
   * Returns the child node of the specified edge or -1 if the edge does
   * not exist
   *
   * @param edge the edge
   */
  int childNode(int edge) const { return _childNodes[edge]; }
  /**
   * Returns the edge that is the outgoing edge of the specified
   * parent parent node having the specified symbol, or
   * returns -1 if no such edge exists.
   *
   * @param parentNode the parent node
   * @param symbol symbol labeling the outgoing edge
   */
  int outEdge(int parentNode, int symbol) const { return _outEdges[symbol][parentNode]; }
private:
  void checkParameters(const HapsMarkerIterator &data, const QVector<float> weights) const;
  void checkParameters(MergeableDagLevel *parent, const HapsMarkerIterator &data) const;
  void initializeArrays();
  void fillArraysForFirst(const HapsMarkerIterator &data);
  void fillArraysWithPrev(const HapsMarkerIterator &data);
  void addEdge(int parentNode, int symbol, float weight, int edge, int haplotype);
  void reduceEdgeArrayLengths(int newLength);
  void removeHaplotypeIndices();  // Removes the haplotype index data from this object.
  bool hasSiblingInMe(int parentNodeOfNextLevel) const;
  void mergeParentNodesForward(int retainedNode, int removedNode);
  bool isParentNode(int node) const;

  /*
   * Merges the two specified child nodes and assigns the merged
   * node to the specified {@code retainedNode} index.  Ingoing edges
   * to {@code removedNode} are redirected to be ingoing edges
   * to {@code retainedNode}.
   *
   * @param retainedNode a child node which will receive ingoing edges of
   * {@code removedNode}
   * @param removedNode a child node that will be deleted after merging
   */
  void mergeChildNodes(int retainedNode, int removedNode);

  void changeParent(int edge, int newParent);
  void mergeEdges(int retainedEdge, int removedEdge);
  void mergeHaplotypes(int retainedChild, int removedChild);

  MergeableDagLevel *_nextLevel;  // Pointer to the next element of the doubly-linked list.
  MergeableDagLevel *_prevLevel;  // Pointer to the previous element of the doubly-linked list.

  int _levelIndex;
  int _nAlleles;
  int _nHaps;
  QVector<float> _weights;

  QList<QVector<int> > _outEdges;  // [allele][parent node]
  QVector<int> _child2FirstInEdge;
  QVector<int> _inEdge2NextInEdge;

  QVector<int> _parentNodes;  // edge -> parent node
  QVector<int> _childNodes;   // edge -> child node
  QVector<int> _symbols;      // edge -> symbol
  QVector<float> _counts;     // edge -> weight

  QVector<int> _child2FirstHap;  // child node -> first hap index
  QVector<int> _hap2NextHap;     // current hap index -> next hap index
};

/**
 * Class {@code Score} represents a similarity score for a pair
 * of trees.
 *
 * "Instances of class {@code Score} are immutable."
 */
class Score
{
public:
  /**
   * Constructs a new {@code Score} instance.  Smaller similarity scores
   * correspond to greater similarity.
   * @param nodeA root node index for the first tree
   * @param nodeB root node index for the second tree
   * @param score the a non-negative similarity score for the two specified
   * trees
   * @param isMergeable {@code true} if the two trees may be
   * merged, and {@code false} otherwise
   */
  Score(int nodeA, int nodeB, float score, bool isMergeable);

  /**
   * Default constructor for {@code Score} is used to create the
   * maximum score value.
   */
  Score();

  /**
   * Returns the root node index for the first tree.
   */
  int nodeA() const { return _nodeA; }
  /**
   * Returns the root node index for the second tree.
   */
  int nodeB() const { return _nodeB; }
  /**
   * Returns the similarity score for the two trees.
   */
  float score() const { return (_score < 0) ? -_score : _score; }
  /**
   * Returns {@code true} if the two trees may be merged, and
   * returns {@code false} otherwise.
   */
  bool isMergeable() const { return _score >= 0; }
private:
  int _nodeA;
  int _nodeB;
  float _score;
};

/**
 * Class {@code MergeableDagFactory} contains a static, thread-safe factory
 * method that constructs a Directed Acyclic Graph (DAG) from sequence data.
 *
 * References:
 * <br>
 * Ron D, Singer Y, and Tishby N (1998) On the Learnability and
 * usage of acyclic probabilistic finite automata.  Journal of Computer
 * and SystemSciences 56:133-152.
 * <br>
 * Browning S (2006) Multi-locus association mapping using variable length
 * Markov chains.  Am J Hum Genet. 78:903-913.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
class MergeableDagFactory
{
public:
  /**
   * From the specified data, constructs the (remaining) information
   * (specifically, the required list of ImmutableDagLevel objects)
   * needed to construct a new {@code ImmutableDag} instance.
   *
   * @param hapPairs the sequence data
   * @param weights an array whose {@code j}-th element is the
   * weight for the {@code j}-th haplotype
   * @param scale a parameter that multiplicatively scales the node
   * similarity threshold
   * @param nInitLevels the number of initial levels to read
   */
  MergeableDagFactory(const HapPairs &hapPairs, const QVector<float> &weights, float scale, int nInitLevels);

  const QList<DagLevel> &levels() const { return _mergedLevels; }
private:
  void checkParameters(const HapPairs &hapPairs, const QVector<float> &weights, float scale, int nInitLevels) const;
  MergeableDagLevel *readFirstLevel(HapsMarkerIterator &it, const QVector<float> &weights) const;
  float maxUnmergedAtLeaf(const HapPairs &hapPairs, const QVector<float> &weights) const;
  int nextReadDepth(float unmergedRatio, int depth, int lastDepth) const;
  MergeableDagLevel *readLevels(HapsMarkerIterator &it, int nLevelsToRead,
                                MergeableDagLevel *leafLevel) const;
  void mergeParentNodes(MergeableDagLevel *level);
  void getPairwiseScores(Score &scoreObj, MergeableDagLevel *level, QLinkedList<Score> &scores);
  void findSortedParents(QList<int> &sortedParentNodes, int &nParentsWithSibs, MergeableDagLevel *level) const;
  bool score(Score &scoreObj, MergeableDagLevel *level, int nodeA, int nodeB);

  /*
   * Returns a similarity-score (lower scores correspond to higher
   * similarity).
   *
   * @param marker marker DAG marker containing the specified nodes
   * @param nodeA a parent node at the specified DAG level in tree A
   * @param nodeB a parent node at the specified DAG level in tree B
   * @param nodeCntA the count of the parent node in tree A
   * @param nodeCntB the count of the parent node in tree B
   * @param baseMarker the marker index at the root of trees A and B
   * @param nA the node count of the root of tree A
   * @param nB the node count of the root of tree B
   * @param maxDiff the current maximum difference in proportions in
   * the counts of corresponding tree branches
   * @param threshold the maximum permitted node similarity
   */
  float similar(MergeableDagLevel *level, int nodeA, int nodeB, float nodeCntA, float nodeCntB,
                 int baseMarker, float nA, float nB, float maxDiff, float threshold);

  QList<DagLevel> _mergedLevels;
  float _scale;
  float _nUnmergedAtLeaf;
};

/**
 * Utilities for MergeableDagLevel....
 */

static QVector<quint16> toQuint16Array(const QVector<int> &ia)
{
  int n = ia.length();
  QVector<quint16> ca(n);

  for (int j = 0; j < n; ++j) {
    Q_ASSERT_X(ia[j] >= 0 && ia[j] <= UINT16_MAX, "toQuint16Array (dag.cpp)", "value out of range");

    ca[j] = (quint16)ia[j];
  }

  return ca;
}

/*
 * Returns an array obtained by replacing each array value with its
 * rank when the set of array values is ordered: the smallest value
 * is replaced by 0, the next smallest value is replaced by 1, etc.
 */
static QVector<quint16> rankedQuint16Values(const QVector<int> &array)
{
  Q_ASSERT_X(array.length() > 0, "rankedQuint16Values (dag.cpp)", "array.length() == 0");

  QVector<int> sortedCopy = array;
  qSort(sortedCopy.begin(), sortedCopy.end());
  Q_ASSERT_X(
      sortedCopy[0] >= 0, "rankedQuint16Values (dag.cpp)", "sortedCopy has negative elements");

  int n = sortedCopy[sortedCopy.length() - 1] + 1;
  QVector<int> indexMap(n, 0);
  int index = 0;
  indexMap[sortedCopy[0]] = index++;
  for (int j = 1, n = sortedCopy.length(); j < n; ++j) {
    if (sortedCopy[j] != sortedCopy[j - 1])
      indexMap[sortedCopy[j]] = index++;
  }

  Q_ASSERT_X(index <= UINT16_MAX,
             "rankedQuint16Values (dag.cpp)",
             "Array has more than UINT16_MAX values");

  QVector<quint16> transformedArray(array.length());
  for (int j = 0, n = transformedArray.length(); j < n; ++j)
    transformedArray[j] = (quint16)indexMap[array[j]];

  return transformedArray;
}

MergeableDagLevel::MergeableDagLevel(const HapsMarkerIterator &data, QVector<float> weights)
{
  checkParameters(data, weights);
  _prevLevel = 0;
  _nextLevel = 0;
  _levelIndex = 0;
  _nAlleles = data.marker().nAlleles();
  _nHaps = data.nHaps();
  _weights = weights;
  initializeArrays();
  fillArraysForFirst(data);
}

MergeableDagLevel::MergeableDagLevel(MergeableDagLevel *prevLevel, const HapsMarkerIterator &data)
{
  checkParameters(prevLevel, data);
  _prevLevel = prevLevel;
  _nextLevel = 0;
  _levelIndex = prevLevel->index() + 1;
  _nAlleles = data.marker().nAlleles();
  _nHaps = data.nHaps();
  _weights = prevLevel->_weights;
  initializeArrays();
  fillArraysWithPrev(data);
}

void MergeableDagLevel::checkParameters(const HapsMarkerIterator &data, const QVector<float> weights) const
{
  Q_ASSERT_X(weights.length() == data.nHaps(),
             "MergeableDagLevel::checkParameters(data, weights)",
             "data.nHaps() != weights.length");
}

void MergeableDagLevel::checkParameters(MergeableDagLevel *parent, const HapsMarkerIterator &data) const
{
  Q_ASSERT_X(!parent->_nextLevel,
             "MergeableDagLevel::checkParameters(parent, data)",
             "parent->nextLevel!=null");

  Q_ASSERT_X(parent->nHaps() == data.nHaps(),
             "MergeableDagLevel::checkParameters(parent, data)",
             "inconsistent samples");

  // NB: the sequences of sample ID indices are not checked
}

void MergeableDagLevel::initializeArrays()
{
  QVector<int> oe(_nHaps, -1);
  for (int j = 0; j < _nAlleles; j++)
    _outEdges.append(oe);

  _child2FirstInEdge.fill(-1, _nHaps);
  _inEdge2NextInEdge.fill(-1, _nHaps);
  _parentNodes.fill(-1, _nHaps);
  _childNodes.fill(-1, _nHaps);
  _symbols.fill(-1, _nHaps);
  _counts.fill(0.0, _nHaps);
  _child2FirstHap.fill(-1, _nHaps);
  _hap2NextHap.fill(-1, _nHaps);
}

void MergeableDagLevel::fillArraysForFirst(const HapsMarkerIterator &data)
{
  int parentNode = 0;
  int nEdges = 0;
  for (int hap = 0, n = data.nHaps(); hap < n; ++hap) {
    int symbol = data.allele(hap);
    float count = _weights[hap];
    int edge = _outEdges[symbol][parentNode];
    if (edge == -1)
      addEdge(parentNode, symbol, count, nEdges++, hap);
    else {
      Q_ASSERT_X(
          edge == _childNodes[edge], "MergeableDagLevel::fillArrays", "edge!=childNodes[edge]");
      int child = _childNodes[edge];
      _counts[edge] += count;
      _hap2NextHap[hap] = _child2FirstHap[child];
      _child2FirstHap[child] = hap;
    }
  }

  if (nEdges < 0.75 * _nHaps)
    reduceEdgeArrayLengths(nEdges);
}

void MergeableDagLevel::fillArraysWithPrev(const HapsMarkerIterator &data)
{
  int nEdges = 0;
  for (int node = 0, n = _prevLevel->_child2FirstHap.length(); node < n; ++node) {
    if (_prevLevel->_child2FirstHap[node] >= 0) {
      int hap = _prevLevel->_child2FirstHap[node];
      while (hap != -1) {
        int symbol = data.allele(hap);
        float count = _weights[hap];
        int edge = _outEdges[symbol][node];
        if (edge == -1)
          addEdge(node, symbol, count, nEdges++, hap);
        else {
          Q_ASSERT_X(edge == _childNodes[edge],
                     "MergeableDagLevel::fillArraysWithPrev",
                     "edge != _childNodes[edge]");

          int child = _childNodes[edge];
          _counts[edge] += count;
          _hap2NextHap[hap] = _child2FirstHap[child];
          _child2FirstHap[child] = hap;
        }
        hap = _prevLevel->_hap2NextHap[hap];
      }
    }
  }

  if (nEdges < 0.75 * _nHaps)
    reduceEdgeArrayLengths(nEdges);

  _prevLevel->removeHaplotypeIndices();
}

void MergeableDagLevel::addEdge(int parentNode, int symbol, float weight, int edge, int haplotype)
{
  int childNode = edge;
  _outEdges[symbol][parentNode] = edge;
  _child2FirstInEdge[childNode] = edge;
  _parentNodes[edge] = parentNode;
  _childNodes[edge] = childNode;
  _symbols[edge] = symbol;
  _counts[edge] = weight;
  _child2FirstHap[childNode] = haplotype;
}

void MergeableDagLevel::reduceEdgeArrayLengths(int newLength)
{
  _child2FirstInEdge.resize(newLength);
  _inEdge2NextInEdge.resize(newLength);
  _parentNodes.resize(newLength);
  _childNodes.resize(newLength);
  _symbols.resize(newLength);
  _counts.resize(newLength);
}

void MergeableDagLevel::removeHaplotypeIndices()
{
  _child2FirstHap.clear();
  _hap2NextHap.clear();
}

MergeableDagLevel *MergeableDagLevel::takePrevious()
{
  MergeableDagLevel *prev = _prevLevel;
  _prevLevel = 0;
  return prev;
}

void MergeableDagLevel::setNextLevel(MergeableDagLevel *newNextLevel)
{
  Q_ASSERT_X(newNextLevel->_prevLevel == this,
             "MergeableDagLevel::setNextLevel",
             "nextLevel->previousLevel != this");

  _nextLevel = newNextLevel;
}

bool MergeableDagLevel::hasSiblingInMe(int parentNodeOfNextLevel) const
{
  int edge = _child2FirstInEdge[parentNodeOfNextLevel];

  while (edge >= 0) {
    int pn = _parentNodes[edge];
    int cnt = 0;

    for (int allele = 0; allele < _nAlleles; ++allele) {
      if (_outEdges[allele][pn] >= 0)
        ++cnt;
    }

    if (cnt > 1)
      return true;

    edge = _inEdge2NextInEdge[edge];
  }

  return false;
}

DagLevel MergeableDagLevel::toDagLevel()
{
  _counts.removeAll(0.0);
  _symbols.removeAll(-1);
  _parentNodes.removeAll(-1);
  _childNodes.removeAll(-1);

  Q_ASSERT_X(_counts.length() <= UINT16_MAX, "MergeableDagLevel::toDagLevel", "too many counts");

  DagLevel newdl(rankedQuint16Values(_parentNodes), rankedQuint16Values(_childNodes),
                 toQuint16Array(_symbols), _counts);
  return newdl;
}

void MergeableDagLevel::mergeParentNodes(int retainedNode, int removedNode)
{
  Q_ASSERT_X(isParentNode(retainedNode),
             "MergeableDagLevel::mergeParentNodes",
             "invalid retained parent node");

  Q_ASSERT_X(isParentNode(removedNode),
             "MergeableDagLevel::mergeParentNodes",
             "invalid removed parent node");

  _prevLevel->mergeChildNodes(retainedNode, removedNode);
  mergeParentNodesForward(retainedNode, removedNode);
}

void MergeableDagLevel::mergeParentNodesForward(int retainedNode, int removedNode)
{
  for (int j = 0; j < _nAlleles; ++j) {
    int retainedEdge = _outEdges[j][retainedNode];
    int removedEdge = _outEdges[j][removedNode];
    if (removedEdge >= 0) {
      if (retainedEdge == -1)
        changeParent(removedEdge, retainedNode);
      else {
        int retainedChild = childNode(retainedEdge);
        int removedChild = childNode(removedEdge);
        mergeEdges(retainedEdge, removedEdge);
        if (_nextLevel)
          _nextLevel->mergeParentNodesForward(retainedChild, removedChild);
      }
    }
  }
}

void MergeableDagLevel::mergeChildNodes(int retainedNode, int removedNode)
{
  int lastEdge = -1;
  int edge = _child2FirstInEdge[removedNode];

  while (edge != -1) {
    Q_ASSERT_X(_childNodes[edge] == removedNode,
               "MergeableDagLevel::mergeChildNodes",
               "_childNodes[edge] != removedNode");

    _childNodes[edge] = retainedNode;
    lastEdge = edge;
    edge = _inEdge2NextInEdge[edge];
  }

  if (lastEdge != -1) {
    _inEdge2NextInEdge[lastEdge] = _child2FirstInEdge[retainedNode];
    _child2FirstInEdge[retainedNode] = _child2FirstInEdge[removedNode];
    _child2FirstInEdge[removedNode] = -1;
  }
}

void MergeableDagLevel::changeParent(int edge, int newParent)
{
  int oldParent = _parentNodes[edge];
  int symbol = _symbols[edge];

  Q_ASSERT_X(_outEdges[symbol][oldParent] == edge,
             "MergeableDagLevel::changeParent",
             "_outEdges[symbol][oldParent] != edge");

  Q_ASSERT_X(_outEdges[symbol][newParent] == -1,
             "MergeableDagLevel::changeParent",
             "_outEdges[symbol][newParent] != -1");

  _outEdges[symbol][oldParent] = -1;
  _outEdges[symbol][newParent] = edge;
  _parentNodes[edge] = newParent;
}

void MergeableDagLevel::mergeEdges(int retainedEdge, int removedEdge)
{
  Q_ASSERT_X(_symbols[retainedEdge] == _symbols[removedEdge],
             "MergeableDagLevel::mergeEdges",
             "_symbols[retainedEdge] != _symbols[removedEdge]");

  Q_ASSERT_X(
      _counts[removedEdge] > 0.0, "MergeableDagLevel::mergeEdges", "_counts[removedEdge] <= 0.0");

  _counts[retainedEdge] += _counts[removedEdge];

  if (!_nextLevel)
    mergeHaplotypes(_childNodes[retainedEdge], _childNodes[removedEdge]);

  int parentNode = _parentNodes[removedEdge];
  int childNode = _childNodes[removedEdge];
  int symbol = _symbols[removedEdge];

  Q_ASSERT_X(_inEdge2NextInEdge[_child2FirstInEdge[childNode]] == -1,
             "MergeableDagLevel::mergeEdges",
             "_inEdge2NextInEdge[_child2FirstInEdge[childNode]] != -1");

  _outEdges[symbol][parentNode] = -1;
  _child2FirstInEdge[childNode] = -1;
  _counts[removedEdge] = 0.0;
  _parentNodes[removedEdge] = -1;
  _childNodes[removedEdge] = -1;
  _symbols[removedEdge] = -1;
}

void MergeableDagLevel::mergeHaplotypes(int retainedChild, int removedChild)
{
  int hap = _child2FirstHap[removedChild];

  while (_hap2NextHap[hap] != -1)
    hap = _hap2NextHap[hap];

  _hap2NextHap[hap] = _child2FirstHap[retainedChild];
  _child2FirstHap[retainedChild] = _child2FirstHap[removedChild];
  _child2FirstHap[removedChild] = -1;
}

float MergeableDagLevel::nodeCount(int parentNode) const
{
  float sum = 0.0;
  for (int symbol = 0; symbol < _nAlleles; ++symbol) {
    if (_outEdges[symbol][parentNode] >= 0)
      sum += edgeCount(_outEdges[symbol][parentNode]);
  }
  return sum;
}

QList<int> MergeableDagLevel::parentNodeArray() const
{
  QList<int> sortedReducedArray = _parentNodes.toList();
  sortedReducedArray.removeAll(-1);
  qSort(sortedReducedArray.begin(), sortedReducedArray.end());

  Q_ASSERT_X(sortedReducedArray.length() > 0,
             "MergeableDagLevel::parentNodeArray",
             "sortedReducedArray.length() == 0");

  QList<int> uniqueNodes;
  uniqueNodes.append(sortedReducedArray[0]);
  for (int j = 1, n = sortedReducedArray.length(); j < n; ++j) {
    if (sortedReducedArray[j] != sortedReducedArray[j - 1])
      uniqueNodes.append(sortedReducedArray[j]);
  }

  return uniqueNodes;
}

bool MergeableDagLevel::isParentNode(int node) const
{
  if (_prevLevel)
    return (_prevLevel->_child2FirstInEdge[node] >= 0);
  else {
    for (int j = 0; j < _nAlleles; ++j) {
      if (_outEdges[j][node] != -1)
        return true;
    }
    return false;
  }
}

/**
 * Utility for DagLevel and MergeableDagFactory
 */
static float sumArray(QVector<float> fa)
{
  float sum = 0.0;
  foreach (const float f, fa)
    sum += f;
  return sum;
}

Score::Score(int nodeA, int nodeB, float score, bool isMergeable) : _nodeA(nodeA), _nodeB(nodeB)
{
  Q_ASSERT_X(score >= 0  &&  (score != 0  ||  isMergeable != false),
             "Score::Score",
             "invalid value for score");

  _score = isMergeable ? score : -score;
}

Score::Score() : _nodeA(-1), _nodeB(-1), _score(-FLT_MAX)
{
}

MergeableDagFactory::MergeableDagFactory(const HapPairs &hapPairs, const QVector<float> &weights, float scale,
                                         int nInitLevels)
{
  checkParameters(hapPairs, weights, scale, nInitLevels);
  _scale = scale;

  float maxUnmerged = maxUnmergedAtLeaf(hapPairs, weights);
  int lastReadDepth = nInitLevels;

  // "current" and "leaf" point to a singly-linked list of
  // MergeableDagLevel elements.
  HapsMarkerIterator it(hapPairs);
  MergeableDagLevel *current = readFirstLevel(it, weights);
  MergeableDagLevel *leaf = readLevels(it, nInitLevels, current);

  while (current->next()) {
    _nUnmergedAtLeaf = 0.0;
    current = current->next();
    mergeParentNodes(current);

    MergeableDagLevel *previousLevel = current->takePrevious();
    _mergedLevels.append(previousLevel->toDagLevel());
    delete previousLevel;

    if (it.hasNext()) {
      float ratio = (_nUnmergedAtLeaf / maxUnmerged);
      int depth = (leaf->index() - current->index());
      int readDepth = nextReadDepth(ratio, depth, lastReadDepth);

      if (readDepth > depth) {
        leaf = readLevels(it, (readDepth - depth), leaf);
        lastReadDepth = readDepth;
      }
    }
  }

  _mergedLevels.append(current->toDagLevel());
  delete current;
}

void MergeableDagFactory::checkParameters(const HapPairs &hapPairs, const QVector<float> &weights, float scale,
                                          int nInitLevels) const
{
  Q_ASSERT_X(nInitLevels >= 1, "MergeableDagFactory::checkParameters", "nInitLevels < 1");

  Q_ASSERT_X(
      hapPairs.nMarkers() > 0, "MergeableDagFactory::checkParameters", "hapPairs.nMarkers()==0");

  for (int j = 0; j < weights.length(); ++j) {
    Q_ASSERT_X(weights[j] > 0.0, "MergeableDagFactory::checkParameters", "invalid weight");
  }

  Q_ASSERT_X(scale > 0, "MergeableDagFactory::checkParameters", "scale <= 0");
}

MergeableDagLevel *MergeableDagFactory::readFirstLevel(
    HapsMarkerIterator &it, const QVector<float> &weights) const
{
  it.next();
  return new MergeableDagLevel(it, weights);
}

float MergeableDagFactory::maxUnmergedAtLeaf(const HapPairs &hapPairs, const QVector<float> &weights) const
{
  float sum = sumArray(weights);
  return MAX_PROP_UNMERGED * sum;
}

int MergeableDagFactory::nextReadDepth(float unmergedRatio, int depth, int lastDepth) const
{
  if (unmergedRatio <= 1) {
    return MIN_DEPTH;
  } else if ((float)depth < (0.85 * lastDepth)) {
    return 1 + (int)round(0.95 * lastDepth + 0.5);
  } else if ((unmergedRatio > 2.0) && ((float)depth > (0.95 * lastDepth))) {
    return (int)ceil((1 + unmergedRatio / 20) * lastDepth);
  } else {
    return lastDepth;
  }
}

MergeableDagLevel *MergeableDagFactory::readLevels(HapsMarkerIterator &it,
                                                   int nLevelsToRead,
                                                   MergeableDagLevel *leafLevel) const
{
  for (int j = 0; it.hasNext() && j < nLevelsToRead; ++j)
  {
    it.next();
    MergeableDagLevel *newLeaf = new MergeableDagLevel(leafLevel, it);
    leafLevel->setNextLevel(newLeaf);
    leafLevel = newLeaf;
  }

  return leafLevel;
}

void MergeableDagFactory::mergeParentNodes(MergeableDagLevel *level)
{
  QLinkedList<Score> scores;

  Score minScore;  // Initialized by the call to "getPairwiseScores" just below.
  Score maxScore;  // Class "Score" default-constructs to the maximum value.
  Score newScore;  // Use "newScore" as an output variable for method "score".

  getPairwiseScores(minScore, level, scores);

  while (minScore.isMergeable()) {
    int retainedNode = minScore.nodeA();
    int removedNode = minScore.nodeB();
    if (level->hasSibling(retainedNode) == false) {
      // Ensure that no-sibling nodes are always removed
      retainedNode = minScore.nodeB();
      removedNode = minScore.nodeA();
      Q_ASSERT_X(level->hasSibling(retainedNode),
                 "MergeableDagFactory::mergeParentNodes",
                 "retainedNode has no sibling");
    } else if (level->hasSibling(removedNode) &&
               level->nodeCount(minScore.nodeA()) < level->nodeCount(minScore.nodeB())) {
      removedNode = minScore.nodeB();
      retainedNode = minScore.nodeA();
    }
    level->mergeParentNodes(retainedNode, removedNode);
    minScore = maxScore;

    QMutableLinkedListIterator<Score> scoreIt(scores);
    while (scoreIt.hasNext()) {
      Score &s = scoreIt.next();
      if (s.nodeA() == removedNode || s.nodeB() == removedNode)
        scoreIt.remove();
      else {
        if (s.nodeA() == retainedNode || s.nodeB() == retainedNode) {
          if (score(newScore, level, s.nodeA(), s.nodeB())) {
            if (newScore.score() < minScore.score() && newScore.isMergeable())
              minScore = newScore;

            scoreIt.setValue(newScore);
          } else
            scoreIt.remove();
        } else if (s.score() < minScore.score() && s.isMergeable())
          minScore = s;
      }
    }
  }
}

void MergeableDagFactory::getPairwiseScores(Score &scoreObj, MergeableDagLevel *level,
                                            QLinkedList<Score> &scores)
{
  Score minScore;  // Class "Score" default-constructs to the maximum value.
  Score newScore;  // Use "newScore" as an output variable for method "score".

  QList<int> sortedParentNodes;
  int nParentsWithSibs;
  findSortedParents(sortedParentNodes, nParentsWithSibs, level);

  int spnLength = sortedParentNodes.length();

  for (int j = 0; j < nParentsWithSibs; ++j) {
    int nodeA = sortedParentNodes[j];

    for (int k = j + 1; k < spnLength; ++k) {
      int nodeB = sortedParentNodes[k];

      if (score(newScore, level, nodeA, nodeB)) {
        if (newScore.score() < minScore.score() && newScore.isMergeable())
          minScore = newScore;

        scores.append(newScore);
      }
    }
  }

  scoreObj = minScore;
}

void MergeableDagFactory::findSortedParents(QList<int> &sortedParentNodes, int &nParentsWithSibs, MergeableDagLevel *level) const
{
  sortedParentNodes = level->parentNodeArray();
  int index1 = 0;
  int index2 = sortedParentNodes.length() - 1;

  while (index1 < index2) {
    if (level->hasSibling(sortedParentNodes[index1]))
      ++index1;
    else {
      int tmp = sortedParentNodes[index1];
      sortedParentNodes[index1] = sortedParentNodes[index2];
      sortedParentNodes[index2--] = tmp;
    }
  }

  if (level->hasSibling(sortedParentNodes[index1]))
    ++index1;

  nParentsWithSibs = index1;
}

bool MergeableDagFactory::score(Score &scoreObj, MergeableDagLevel *level, int nodeA, int nodeB)
{
  float maxDiff = 0.0;
  float nodeCntA = level->nodeCount(nodeA);
  float nodeCntB = level->nodeCount(nodeB);
  float threshold = (float)(_scale * sqrt((1.0 / nodeCntA) + (1.0 / nodeCntB)));
  maxDiff = similar(level, nodeA, nodeB, nodeCntA, nodeCntB, level->index(), nodeCntA, nodeCntB,
                    maxDiff, threshold);
  if (maxDiff > MAX_THRESHOLD_RATIO * threshold)
    return false;
  else {
    bool isMergeable = (maxDiff < threshold);
    Score newScore(nodeA, nodeB, maxDiff, isMergeable);
    scoreObj = newScore;
    return true;
  }
}

float MergeableDagFactory::similar(MergeableDagLevel *level, int nodeA, int nodeB, float nodeCntA,
                                   float nodeCntB, int baseMarker, float nA, float nB,
                                   float maxDiff, float threshold)
{
  float propA = nodeCntA / nA;
  float propB = nodeCntB / nB;
  float diff = abs(propA - propB);

  if (diff >= threshold)
    return diff;
  else if (propA <= maxDiff && propB <= maxDiff)
    return maxDiff;
  else if (diff > maxDiff)
    maxDiff = diff;

  if ((nodeA == -1) ^ (nodeB == -1))
    return maxDiff;
  else if (!level) {
    _nUnmergedAtLeaf += (nodeCntA + nodeCntB);
    return maxDiff;
  }

  for (int j = 0, n = level->nAlleles(); j < n; ++j) {
    int edgeA = level->outEdge(nodeA, j);
    int edgeB = level->outEdge(nodeB, j);
    int childA = (edgeA != -1) ? level->childNode(edgeA) : -1;
    int childB = (edgeB != -1) ? level->childNode(edgeB) : -1;
    nodeCntA = (edgeA != -1) ? level->edgeCount(edgeA) : 0.0;
    nodeCntB = (edgeB != -1) ? level->edgeCount(edgeB) : 0.0;
    float childMaxDiff = similar(level->next(), childA, childB, nodeCntA, nodeCntB, baseMarker, nA,
                                 nB, maxDiff, threshold);

    if (childMaxDiff > maxDiff) {
      if (childMaxDiff >= threshold)
        return childMaxDiff;
      else
        maxDiff = childMaxDiff;
    }
  }

  return maxDiff;
}

/**
 * Utility routines for class DagLevel....
 */

static int checkLengths(QVector<quint16> parentNodes, QVector<quint16> childNodes,
                        QVector<quint16> symbols, QVector<float> counts)
{
  Q_ASSERT_X(parentNodes.length() <= UINT16_MAX,
             "checkLengths (dag.cpp)",
             "parentNodes.length() > UINT16_MAX");

  Q_ASSERT_X(parentNodes.length() == childNodes.length() && parentNodes.length() == symbols.length() &&
                 parentNodes.length() == counts.length(),
             "checkLengths (dag.cpp)",
             "inconsistent arrays");

  return parentNodes.length();
}

static void checkForDuplicateOutEdges(QVector<quint16> parentIndices, QVector<quint16> parents,
                                      QVector<quint16> symbols)
{
  QSet<int> indexSet;
  indexSet.reserve(symbols.length());

  for (int j = 1, n = parentIndices.length(); j < n; ++j)
  {
    indexSet.clear();
    for (int k = parentIndices[j - 1], n = parentIndices[j]; k < n; ++k) {
      int edge = parents[k];
      Q_ASSERT_X(!indexSet.contains(symbols[edge]),
                 "checkForDuplicateOutEdges (dag.cpp)",
                 "duplicate edge");
      indexSet.insert(symbols[edge]);
    }
  }
}

static int maxElement(QVector<quint16> ca)
{
  quint16 max = 0;
  foreach (const quint16 c, ca) {
    if (c > max) {
      max = c;
    }
  }
  return max;
}

/*
 * Returns an array of length {@code max(nodes) + 1}
 * whose {@code j}-th element is the number of
 * elements of the specified array that have value {@code j}.
 *
 * @param nodeCounts reference to the output array.
 * @param nodes an array of non-negative values.
 * "@throws IllegalArgumenException if set of elements of the
 * specified array is not equal to {@code {0, 1, 2, ..., k}} for some
 * {@code k}."
 */
static void elementCounts(QVector<int> &nodeCounts, const QVector<quint16> array)
{
  int maxNode = maxElement(array);

  nodeCounts.fill(0, maxNode + 1);

  foreach (const quint16 c, array)
    ++nodeCounts[c];

  for (int j = 0, n = nodeCounts.length(); j < n; ++j) {
    Q_ASSERT_X(nodeCounts[j] > 0, "elementCounts (dag.cpp)", "nodeCounts[j] == 0");
  }
}

static void getIndicesArray(QVector<quint16> &indicesArray, const QVector<quint16> nodes)
{
  QVector<int> countArray;
  elementCounts(countArray, nodes);

  indicesArray.fill(0, countArray.length() + 1);
  for (int j = 1, n = indicesArray.length(); j < n; ++j) {
    Q_ASSERT_X(countArray[j - 1] > 0, "getIndicesArray (dag.cpp)", "countArray[j-1] <= 0");
    int x = indicesArray[j - 1] + countArray[j - 1];
    Q_ASSERT_X(x <= UINT16_MAX, "getIndicesArray (dag.cpp)", "x > UINT16_MAX");
    indicesArray[j] = (quint16)x;
  }
}

DagLevel::DagLevel(QVector<quint16> parentNodes, QVector<quint16> childNodes,
                   QVector<quint16> symbols, QVector<float> counts)
{
  int nEdges = checkLengths(parentNodes, childNodes, symbols, counts);

  getIndicesArray(_parentIndices, parentNodes);
  getIndicesArray(_childIndices, childNodes);
  _parentNodes = parentNodes;
  _childNodes = childNodes;
  _symbols = symbols;
  _edgeCounts = counts;
  _condEdgeProbs.fill(0.0, nEdges);
  _parents.fill(0, nEdges);
  _children.fill(0, nEdges);

  QVector<quint16> pIndices = _parentIndices.mid(0, _parentIndices.length() - 1);
  QVector<quint16> cIndices = _childIndices.mid(0, _childIndices.length() - 1);
  obtainParentCounts(parentNodes, counts, pIndices.length());
  _count = sumArray(_parentCounts);

  for (quint16 j = 0; j < nEdges; ++j) {
    quint16 p = parentNodes[j];
    quint16 c = childNodes[j];
    _parents[pIndices[p]++] = j;
    _children[cIndices[c]++] = j;
    _condEdgeProbs[j] = counts[j] / _parentCounts[p];
  }
  checkForDuplicateOutEdges(_parentIndices, _parents, _symbols);
}

void DagLevel::obtainParentCounts(QVector<quint16> parentNodes, QVector<float> counts, int nNodes)
{
  _parentCounts.fill(0.0, nNodes);

  for (int j = 0; j < _condEdgeProbs.length(); ++j) {
    quint16 p = parentNodes[j];
    _parentCounts[p] += counts[j];
  }
}

int DagLevel::outEdge(int parentNode, int outEdgeIndex) const
{
  Q_ASSERT_X(outEdgeIndex >= 0 && outEdgeIndex < nOutEdges(parentNode),
             "DagLevel::outEdge",
             "outEdgeIndex out of bounds");
  return _parents[_parentIndices[parentNode] + outEdgeIndex];
}

int DagLevel::outEdgeBySymbol(int parentNode, int symbol) const
{
  int start = _parentIndices[parentNode];
  int end = _parentIndices[parentNode + 1];
  for (int j = start; j < end; ++j) {
    quint16 edgeIndex = _parents[j];
    if (_symbols[edgeIndex] == symbol)
      return (int)edgeIndex;
  }
  return -1;
}

int DagLevel::inEdge(int childNode, int inEdgeIndex) const
{
  Q_ASSERT_X(inEdgeIndex >= 0 && inEdgeIndex < nInEdges(childNode),
             "DagLevel::inEdge",
             "inEdgeIndex out of bounds");
  return _children[_childIndices[childNode] + inEdgeIndex];
}

/**
 * Utility routine for ImmutableDag:
 */
static double minusLog10CondEdgeProb(const DagLevel &level)
{
  float meanScore = 0.0f;
  for (int e = 0, n = level.nEdges(); e < n; ++e) {
    meanScore += level.edgeProb(e) * level.condEdgeProb(e);
  }
  double d = -log10(meanScore);
  return (d < 0) ? 0.0 : d;
}

ImmutableDag::ImmutableDag(const HapPairs &hapPairs, const QVector<float> &weights, float scale,
                           int nInitLevels)
{
  MergeableDagFactory mdf(hapPairs, weights, scale, nInitLevels);

  _markers = hapPairs.markers();
  _dagLevels = mdf.levels();

  Q_ASSERT_X(_dagLevels.length() > 0, "ImmutableDag::ImmutableDag", "_dagLevels.length()==0");

  Q_ASSERT_X(_dagLevels[0].nParentNodes() == 1,
             "ImmutableDag::ImmutableDag",
             "_dagLevels[0].nParentNodes()!=1");

  _nNodes = 1;
  _nEdges = 0;
  _maxNodes = 0;
  _maxEdges = 0;

  for (int j = 0; j < _dagLevels.length(); ++j)
  {
    Q_ASSERT_X(j == 0  ||  _dagLevels[j - 1].nChildNodes() == _dagLevels[j].nParentNodes(),
               "ImmutableDag::ImmutableDag",
               "j > 0  &&  _dagLevels[j - 1].nChildNodes() != _dagLevels[j].nParentNodes()");

    _nNodes += _dagLevels[j].nChildNodes();
    _nEdges += _dagLevels[j].nEdges();

    if (_dagLevels[j].nChildNodes() > _maxNodes)
      _maxNodes = _dagLevels[j].nChildNodes();

    if (_dagLevels[j].nEdges() > _maxEdges)
      _maxEdges = _dagLevels[j].nEdges();

    double d = minusLog10CondEdgeProb(_dagLevels[j]);

    _posArray.append((j == 0) ? d : (_posArray[j - 1] + d));
  }
}

