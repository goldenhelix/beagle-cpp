#ifndef BAUMHMM_H
#define BAUMHMM_H

#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"
#include "impute/haplotypepair.h"
#include "impute/dag.h"

#include <QMap>
#include <QMapIterator>

#define KEEP_TRACK_OF_ORDER_IN_SINGLENODES

class IntPair
{
public:
  IntPair(int a, int b) { _a=a; _b=b; }
  IntPair(const IntPair &other) { _a=other._a; _b=other._b; }

  bool operator<(const IntPair &other) const
  { return (_a == other._a) ? (_b < other._b) : (_a < other._a); }

  int firstInt() const { return _a; }
  int secondInt() const { return _b; }

private:
  int _a;
  int _b;
};

/**
 * Class {@code SingleNodes} stores ordered node pairs and associated values.
 */
class SingleNodes
{
public:

  /**
   * Creates a new instance of {@code SingleNodes} that is
   * considered to have an initial value of 0 for each ordered node
   * pair.
   */
  SingleNodes() {}

  /**
   * Adds the specified positive value to the stored value of the specified
   * node pair.
   *
   * @param node1 the first node
   * @param node2 the second node
   * @param value the value
   */
  void sumUpdate(int node1, int node2, double value);

#ifdef KEEP_TRACK_OF_ORDER_IN_SINGLENODES
  int node1(int order) const { return _orderedPairs[order].firstInt(); }
  int node2(int order) const { return _orderedPairs[order].secondInt(); }
  double value(int order) const { return _nodes[_orderedPairs[order]]; }
  int size() const { return _orderedPairs.length(); }
#else
  /**
   * Returns a Java-style iterator for the underlying QMap inside of
   * this ({@code SingleNodes}) object.
   */
  QMapIterator<IntPair, double> nodeIterator() const;
#endif

  /**
   * Returns the value of the specified node pair.
   *
   * @param node1 the first node
   * @param node2 the second node
   */
  double value(int node1, int node2) const;

  /**
   * Sets the value of each ordered node pair to 0.
   */
  void clear();

private:
  // NOTE: Only nodes with non-zero values are actually kept in this
  // map.
  QMap<IntPair, double> _nodes;

#ifdef KEEP_TRACK_OF_ORDER_IN_SINGLENODES
  QList<IntPair> _orderedPairs;
#endif
};

/**
 * Class {@code SingleBaumLevel} computes forward and backward Baum
 * values at a level of a hidden Markov model (HMM) whose states are
 * ordered edge pairs of a leveled directed acyclic graph (DAG).
 *
 * "Instances of class {@code SingleBaumLevel} are not thread-safe."
 */
class SingleBaumLevel
{
public:

  /**
   * Constructs a new {@code SingleBaumLevel} instance from the specified
   * data.
   * @param dag the directed acyclic graph that the determines transition
   * probabilities
   * @param gl the emission probabilities
   */
  SingleBaumLevel(Dag *dag, SplicedGL *gl);

  /**
   * Sets the Baum forward algorithm values for this level of the HMM
   * and records the child node pair values in the specified
   * {@code nodes} parameter. When the method call returns, the {@code nodes}
   * parameter will be reset to the child node pair values for this level of
   * the HMM.
   *
   * @param nodes child node pair values at the previous level of the HMM
   * @param marker the level of the HMM at which the Baum forward algorithm
   * values will be computed
   * @param sample a sample index
   */
  void setForwardValues(SingleNodes &nodes, int marker, int sample);

  /**
   * Stores the Baum forward algorithm child node pair values for this
   * level of the HMM in the specified {@code SingleNodes} object.
   *
   * @param nodes the node pair values that will be set
   */
  void setChildNodes(SingleNodes &nodes);

  /**
   * Sets the Baum backward algorithm values for this level of the HMM
   * and stores the parent node pair values in the specified
   * {@code nodes} parameter.  When the method call returns, the
   * {@code nodes} parameter will be reset to the parent
   * node pair values for this level of the HMM.
   *
   * @param nodes parent node pair values at the next level of HMM
   */
  // void setBackwardValues(SingleNodes &nodes);

  /**
   * Returns the directed acyclic graph that determines the transition
   * probabilities.
   */
  Dag &dag() const { return *_dag; }

  /**
   * Returns the emission probabilities.
   */
  SplicedGL &gl() const { return *_gl; }

  /**
   * Return the level of the HMM.
   */
  int marker() const { return _marker; }

  /**
   * Return the number of possible genotypes at this level of the HMM.
   */
  // int nGenotypes() const { return _nGenotypes; }

  /**
   * Return the number of states with nonzero forward probability at
   * this level of the HMM.
   */
  int size() const { return _size; }

  /**
   * Returns the first edge of the specified HMM state with nonzero forward
   * probability.
   * @param state an index of a HMM state with nonzero forward probability
   */
  int edge1(int state) const {
    checkIndex(state);
    return _edges1[state];
  }

  /**
   * Returns the second edge of the specified HMM state with nonzero forward
   * probability.
   * @param state an index of a HMM state with nonzero forward probability
   */
  int edge2(int state) const {
    checkIndex(state);
    return _edges2[state];
  }

  /**
   * Returns the parent node of the first edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int parentNode1(int state) const {
    checkIndex(state);
    return _dag->parentNode(_marker, _edges1[state]);
  }

  /**
   * Returns the parent node of the second edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int parentNode2(int state) const {
    checkIndex(state);
    return _dag->parentNode(_marker, _edges2[state]);
  }

  /**
   * Returns the child node of the first edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int childNode1(int state) const {
    checkIndex(state);
    return _dag->childNode(_marker, _edges1[state]);
  }

  /**
   * Returns the child node of the second edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int childNode2(int state) const {
    checkIndex(state);
    return _dag->childNode(_marker, _edges2[state]);
  }

  /**
   * Returns the symbol for the first edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int symbol1(int state) const { return _dag->symbol(_marker, edge1(state)); }

  /**
   * Returns the symbol for the second edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int symbol2(int state) const { return _dag->symbol(_marker, edge2(state)); }

  /**
   * Returns the normalized forward value for the specified HMM state
   * with nonzero forward probability.
   * The normalized forward value is obtained by dividing the
   * forward value by the sum of the forward values at this level
   * of the HMM.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  double forwardValue(int state) const {
    checkIndex(state);
    return _fwdValues[state];
  }

  /**
   * Returns the normalized backward value for the specified HMM state
   * with nonzero forward probability.
   * The normalized backward value is obtained by dividing the
   * backward value by the sum of the backward values at this level
   * of the HMM.
   *
   * @param state an index of a state with nonzero forward probability
   */
  double backwardValue(int state) const {
    checkIndex(state);
    return _bwdValues[state];
  }

  /**
   * Returns the sum of the forward values at this level of the HMM
   * when the forward values are computed using forward values
   * from the previous level that are normalized to sum to 1.
   */
  double forwardValuesSum() const {
    return _fwdValueSum;
  }

  /**
   * Returns the sum of the backward values at this level of the HMM
   * when the backward values are computed using backward
   * values from the next level that are normalized to sum to 1.
   */
  double backwardValuesSum() const {
    return _bwdValueSum;
  }

private:
  void setStates(const SingleNodes &nodes);
  void checkIndex(int state) const;

  Dag *_dag;
  SplicedGL *_gl;

  int _marker;
  int _sample;
  int _size;

  QList<int> _edges1;
  QList<int> _edges2;
  QList<double> _fwdValues;
  QList<double> _bwdValues;
  double _fwdValueSum;
  double _bwdValueSum;
};


/**
 * <p>Class {@code SingleBaum} implements the Baum forward and backward
 * algorithms for a hidden Markov model (HMM) of an individual's genotype data.
 * </p>
 * "Instances of class {@code SingleBaum} are not thread-safe."
 */
class SingleBaum
{
public:

  /**
   * Creates a new {@code SingleBaum} instance from the specified data.
   *
   * @param dag the directed acyclic graph that determines the
   * transition probabilities
   * @param gl the emission probabilities
   * @param seed the random seed
   * @param nSamplingsPerIndividual the number of haplotype pairs that
   * will be sampled for each individual
   * @param lowMem {@code true} if a low memory algorithm should be used, and
   * {@code false} otherwise
   */
  SingleBaum(Dag &dag, SplicedGL &gl, int seed, int nSamplingsPerIndividual,
             bool lowMem);

  QList<HapPair> randomSample(int sample);

  // QList<HapPair> randomSample(int sample, QList<double> gtProbs);

  Dag &dag() const {
    return *_dag;
  }

  SplicedGL &gl() const {
    return *_gl;
  }

  int nSamplingsPerIndividual() const {
    return _nSamplingsPerIndividual;
  }

  int seed() const {
    return _seed;
  }

private:
  QList<HapPair> hapList(int sample) const;
  void initSampleAlleles(const SingleBaumLevel &level, int sample);
  int initialRandomState(const SingleBaumLevel &level, int copy);
  double parentSum(const SingleBaumLevel &level, int sample, int state) const;
  void sampleAlleles(const SingleBaumLevel &level, int sample);
  int randomPreviousState(const SingleBaumLevel &level, int node1,
			  int node2, double nodeValue, int copy);
  SingleBaumLevel &nextLevel();

  SingleBaumLevel &currentLevel() { return _levels[_arrayIndex]; }

  SingleBaumLevel &previousLevel(int sample);
  void forwardAlgorithm(int sample);

  Dag *_dag;
  SplicedGL *_gl;
  int _nMarkers;
  int _nSamplingsPerIndividual;
  long _seed;
  // Random _random;

  QList<int> _node1;
  QList<int> _node2;
  QList<double> _nodeValue;

  QList< QList<int> > _alleles1;
  QList< QList<int> > _alleles2;

  QList<SingleBaumLevel> _levels;
  SingleNodes _fwdNodes;

  int _windowIndex;
  int _arrayIndex;
};

#endif
