#ifndef BAUMHMM_H
#define BAUMHMM_H

#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"
#include "impute/haplotypepair.h"
#include "impute/dag.h"

/**
 * Class {@code SingleNodes} stores ordered node pairs and associated values.
 *
 * "Instances of class {@code SingleNodes} are not thread safe."
 */
class SingleNodes
{
public:

  /**
   * Creates a new instance of {@code SingleNodes} that has an
   * initial value of 0 for each ordered node pair. The first node
   * has index 0.
   */
  SingleNodes();

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
   * Returns the first node of the specified node pair in the list of
   * node pairs with non-zero value.  Repeated invocations of this
   * method with the same parameter will return the same value if
   * node values are not modified between invocations. If
   * {@code (index >= 0 && index < this.size())}, then the following
   * expression will always evaluate to {@code true}:<br>
   * {@code (this.value(this.enumNode1(index),
   * this.enumNode2(index)) == this.enumValue(index))}.
   *
   * @param index an index in a list of node pairs with non-zero
   * value
   * @return the first node of the specified node pair in a list of
   * node pairs with non-zero value
   */
  int enumNode1(int index) const {
    checkSize(index);
    return _node1[_index[index]];
  }

  /**
   * Returns the second node of the specified node pair in a list of
   * node pairs with non-zero value.  Repeated invocations of this
   * method with the same parameter will return the same value if
   * node values are not modified between invocations. If
   * {@code (index >= 0 && index < this.size())}, then the following
   * expression will always evaluate to {@code true}:<br>
   * {@code (this.value(this.enumNode1(index),
   * this.enumNode2(index)) == this.enumValue(index))}.
   *
   * @param index an index in a list of node pairs with non-zero value
   * @return the second node of the specified node pair in a list of
   * node pairs with non-zero value
   */
  int enumNode2(int index) const {
    checkSize(index);
    return _node2[_index[index]];
  }

  /**
   * Returns the value of the specified node pair in a list of
   * node pairs with non-zero value.  Repeated invocations of this
   * method with the same parameter will return the same value if
   * node values are not modified between invocations. If
   * {@code (index >= 0 && index < this.size())}, then the following
   * expression will always evaluate to {@code true}:<br>
   * {@code (this.value(this.enumNode1(index),
   * this.enumNode2(index)) == this.enumValue(index))}.
   *
   * @param index an index in a list of node pairs with non-zero value
   * @return the value of the specified ordered node pair in a list of
   * node pairs with non-zero value
   */
  float enumValue(int index) const {
    checkSize(index);
    return _value[_index[index]];
  }

  /**
   * Returns the value of the specified node pair.
   *
   * @param node1 the first node
   * @param node2 the second node
   * @return the value of the specified node pair
   */
  float value(int node1, int node2) const;

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

  int _size;
  int _capacity; // required to be a power of 2
  int _rehashThreshold;

  QVector<int> _index;
  QVector<int> _node1;
  QVector<int> _node2;
  QVector<float> _value;
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
  SingleBaumLevel(const Dag &dag, const SplicedGL &gl);

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
  const Dag &dag() const { return _dag; }

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
  void setStates(const SingleNodes &nodes);
  void checkIndex(int state) const;

  const Dag &_dag;
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
  SingleBaum(const Dag &dag, const SplicedGL &gl, int seed, int nSamplingsPerIndividual,
             bool lowMem);

  QList<HapPair> randomSample(int sample);

  // QList<HapPair> randomSample(int sample, QList<double> gtProbs);

  const Dag &dag() const {
    return _dag;
  }

  const SplicedGL &gl() const {
    return _gl;
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

  const Dag &_dag;
  const SplicedGL &_gl;
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
