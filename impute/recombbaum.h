#ifndef RECOMBBAUM_H
#define RECOMBBAUM_H

#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"
#include "impute/haplotypepair.h"
#include "impute/dag.h"

/**
 * Class {@code RecombSingleNodes} stores ordered node pairs and
 * associated values.
 *
 * "Instances of class {@code RecombSingleNodes} are not thread safe."
 */
class RecombSingleNodes
{
public:

  /**
   * Creates a new instance of {@code RecombSingleNodes} that has an
   * initial value of 0 for each ordered node pair. The first node
   * has index 0.
   */
  RecombSingleNodes();

  /**
   * Initialize the number of nodes this object will handle.
   *
   * @param nNodes the maximum number of distinct nodes
   * which will be paired to form ordered node pairs
   */
  void initialize(int nNodes);

  /**
   * Adds the specified positive value to the stored value of the specified
   * node pair.
   *
   * @param node1 the first node
   * @param node2 the second node
   * @param value the value
   */
  void sumUpdate(int node1, int node2, float value);

  /**
   * Returns the number of node pairs with non-zero value.
   */
  int size() const { return _size; }

  /**
   * Returns the number of nodes.
   */
  int nNodes() { return _nNodes; }

  /**
   * Returns the value of the specified node pair.
   *
   * @param node1 the first node
   * @param node2 the second node
   * @return the value of the specified node pair
   */
  float value(int node1, int node2) const;

  /**
   * Returns the sum of the values of the node pairs that have the specified
   * first node
   *
   * @param node1 a node
   */
  float sumNode1Value(int node1) const { return _sumNode1Value[node1]; }

  /**
   * Returns the sum of the values of the node pairs that have the specified
   * second node.
   *
   * @param node2 a node
   */
  float sumNode2Value(int node2) const { return _sumNode2Value[node2]; }

  /**
   * Returns the sum of the values of all node pairs.
   */
  float sumValue() const { return _sumValue; }

  /**
   * Sets the value of each ordered node pair to 0.
   */
  void clear();

private:
  /*
   * Return the storage index for specified node pair.  If the key is not
   * currently stored in the hash table, the index at which the value
   * should be stored is returned.
   */
  int index(int node1, int node2) const;

  /*
   * Increases the capacity of the internal hash table.
   */
  void rehash();
  void checkSize(int index) const;

  int _nNodes;

  int _size;
  int _capacity; // required to be a power of 2
  int _rehashThreshold;

  QVector<int> _index;
  QVector<int> _node1;
  QVector<int> _node2;
  QVector<float> _value;
  QVector<float> _sumNode1Value;
  QVector<float> _sumNode2Value;
  float _sumValue;
};

class SamplerData;
class SinglePermittedStates;

/**
 * Class {@code RecombSingleBaumLevel} computes forward and backward Baum
 * values at a level of a hidden Markov model (HMM) whose states are
 * ordered edge pairs of a leveled directed acyclic graph (DAG).
 *
 * "Instances of class {@code RecombSingleBaumLevel} are not thread-safe."
 */
class RecombSingleBaumLevel
{
public:

  /**
   * Constructs a new {@code RecombSingleBaumLevel} instance from the
   * specified data.
   * @param samplerData the analysis data
   */
  RecombSingleBaumLevel(const SamplerData &samplerData);


  /**
   * Resets the size of this level to 0 and clears the QLists in
   * this level.
   */
  void reset();

  /**
   * Sets the Baum forward algorithm values for this level of the HMM
   * and records the child node pair values in the specified
   * {@code nodes} parameter. When the method call returns, the {@code nodes}
   * parameter will be reset to the child node pair values for this level of
   * the HMM.
   *
   * @param nodes child node pair values at the previous level of the HMM
   * @param permittedStates the permitted diploid model states
   * @param marker the level of the HMM at which the Baum forward algorithm
   * values will be computed
   * @param sample a sample index
   */
  void setForwardValues(RecombSingleNodes &nodes,
                        SinglePermittedStates &permittedStates,
                        int marker, int sample);

  /**
   * Stores the Baum forward algorithm child node pair values for this
   * level of the HMM in the specified {@code SingleNodes} object.
   *
   * @param nodes the node pair values that will be set
   */
  void setChildNodes(RecombSingleNodes &nodes);

  /**
   * Initializes the node pair values for the Baum backward algorithm.
   *
   * @param nodes the node pair values to be initialized
   */
  // void setInitialBackwardValues(RecombSingleNodes &nodes);

  /**
   * Sets the Baum backward algorithm values for this level of the HMM
   * and stores the parent node pair values in the specified
   * {@code nodes} parameter.  When the method call returns, the
   * {@code nodes} parameter will be reset to the parent
   * node pair values for this level of the HMM.
   *
   * @param nodes parent node pair values at the next level of HMM
   */
  // void setBackwardValues(RecombSingleNodes &nodes);

  /**
   * Returns the directed acyclic graph that determines the transition
   * probabilities.
   */
  const ImmutableDag &dag() const { return _dag; }

  /**
   * Returns the emission probabilities.
   */
  const SplicedGL &gl() const { return _gl; }

  /**
   * Return the level of the HMM.
   */
  int marker() const { return _marker; }

  /**
   * Return the number of possible genotypes at this level of the HMM.
   */
  // int nGenotypes() const { return _nGenotypes; }

  /**
   * Returns the current capacity of this level.
   */
  // int capacity() const { return _edges1.length(); }

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
    return _dag.parentNode(_marker, _edges1[state]);
  }

  /**
   * Returns the parent node of the second edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int parentNode2(int state) const {
    checkIndex(state);
    return _dag.parentNode(_marker, _edges2[state]);
  }

  /**
   * Returns the child node of the first edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int childNode1(int state) const {
    checkIndex(state);
    return _dag.childNode(_marker, _edges1[state]);
  }

  /**
   * Returns the child node of the second edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int childNode2(int state) const {
    checkIndex(state);
    return _dag.childNode(_marker, _edges2[state]);
  }

  /**
   * Returns the symbol for the first edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int symbol1(int state) const { return _dag.symbol(_marker, edge1(state)); }

  /**
   * Returns the symbol for the second edge of the specified HMM state
   * with nonzero forward probability.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  int symbol2(int state) const { return _dag.symbol(_marker, edge2(state)); }

  /**
   * Returns the normalized forward value for the specified HMM state
   * with nonzero forward probability.
   * The normalized forward value is obtained by dividing the
   * forward value by the sum of the forward values at this level
   * of the HMM.
   *
   * @param state an index of a HMM state with nonzero forward probability
   */
  float forwardValue(int state) const {
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
  float backwardValue(int state) const {
    checkIndex(state);
    return _bwdValues[state];
  }

  /**
   * Returns the sum of the forward values at this level of the HMM
   * when the forward values are computed using forward values
   * from the previous level that are normalized to sum to 1.
   */
  float forwardValuesSum() const {
    return _fwdValueSum;
  }

  /**
   * Returns the sum of the backward values at this level of the HMM
   * when the backward values are computed using backward
   * values from the next level that are normalized to sum to 1.
   */
  float backwardValuesSum() const {
    return _bwdValueSum;
  }

private:
  void setStates(const RecombSingleNodes &nodes,
                 SinglePermittedStates &permittedStates);
  float fwdValue(int edge1, int edge2, const RecombSingleNodes &nodes);
  void checkIndex(int state) const;

  const SamplerData &_samplerData;
  const ImmutableDag &_dag;
  const SplicedGL &_gl;

  int _marker;
  int _sample;
  int _size;

  QList<int> _edges1;
  QList<int> _edges2;
  QList<float> _fwdValues;
  QList<float> _bwdValues;
  float _fwdValueSum;
  float _bwdValueSum;
};

class RestrictedDag;

/**
 * Class {@code RecombSingleBaum} implements the Baum forward and
 * backward algorithms for a hidden Markov model (HMM) of an individual's
 * genotype data. The HMM transition probabilities model recent
 * genetic recombination by allowing jumps between states that are not
 * connected by a node.
 *
 * "Instances of class {@code RecombSingleBaum} are not thread-safe."
 */
class RecombSingleBaum
{
public:

  /**
   * Creates a new {@code RecombSingleBaum} instance from the specified data.
   *
   * @param samplerData the analysis data
   * @param seed the random seed
   * @param nSamplingsPerIndividual the number of haplotype pairs that
   * will be sampled for each individual
   * @param lowMem {@code true} if a low memory algorithm should be used, and
   * {@code false} otherwise
   */
  RecombSingleBaum(const SamplerData &samplerData /* , int seed */ , int nSamplingsPerIndividual,   /// %%%
                   bool lowMem);

  QList<HapPair> randomSample(int sample);

  // QList<HapPair> randomSample(int sample, QList<double> gProbs);

  const ImmutableDag &dag() const { return _dag; }

  const SplicedGL &gl() const { return _gl; }

  int nSamplingsPerIndividual() const {
    return _nSamplingsPerIndividual;
  }

  int seed() const { return _seed; }

private:
  void pruneLevels();
  /// int estMeanSize();
  QList<HapPair> hapList(int sample) const;
  void initSampleAlleles(const RecombSingleBaumLevel &level, int sample);
  int initialRandomState(const RecombSingleBaumLevel &level, int copy);
  void saveCurrentData(const RecombSingleBaumLevel &level, int sample,
                       int copy, int stateIndex);
  void sampleAlleles(const RecombSingleBaumLevel &level, int sample);
  int randomPreviousState(const RecombSingleBaumLevel &level, int copy);
  RecombSingleBaumLevel &nextLevel();

  RecombSingleBaumLevel &currentLevel() { return _levels[_arrayIndex]; }

  RecombSingleBaumLevel &previousLevel(int sample,
                                       SinglePermittedStates &permittedStates);
  void forwardAlgorithm(int sample, SinglePermittedStates &permittedStates);

  const SamplerData &_samplerData;
  const ImmutableDag &_dag;
  const RestrictedDag &_rdag;
  const SplicedGL &_gl;
  int _nMarkers;
  int _nSamplingsPerIndividual;
  long _seed;
  // Random _random;

  QList<int> _node1;
  QList<int> _node2;
  QList<float> _baseTrProb;
  QList<float> _maxSum;
  QList<double> _nodeValue;

  QList< QList<int> > _alleles1;
  QList< QList<int> > _alleles2;

  QList<RecombSingleBaumLevel> _levels;
  RecombSingleNodes _fwdNodes;

  int _windowIndex;
  int _arrayIndex;
};

#endif

