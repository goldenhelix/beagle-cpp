/* Copyright notice. */
#ifndef SAMPLERDATA_H
#define SAMPLERDATA_H

#include "impute/ibshapsegments.h"

class ImmutableDag;

/**
 * Class {@code RestrictedDag} is a wrapper for a {@code ImmutableDag}
 * object that stores segments of identity by descent.
 *
 * "Instances of class {@code RestrictedDag} are immutable."
 */
class RestrictedDag
{
  friend class SinglePermittedStates;

public:

  /**
   * Constructs a {@code RestrictedDag} instance.
   * @param dag the {@code ImmutableDag} object being wrapped
   * @param haps the sample haplotypes
   * @param ibdLength the minimum length of an IBD segment
   * @param ibdExtend the length by which an IBD segment will be extended
   */
  RestrictedDag(const ImmutableDag &dag, const SampleHapPairs &haps,
		double ibdLength, double ibdExtend);

  RestrictedDag& operator=(const RestrictedDag& other);

  /**
   * Returns the haplotypes used to construct {@code this}.
   */
  const SampleHapPairs &sampleHaps() const { return _haps; }

  /**
   * Returns the DAG.
   */
  const ImmutableDag &dag() const { return _dag; }

protected:

  /**
   * Facilitates constructing the object (of class
   * SinglePermittedStates) containing the permitted states for the
   * specified sample.
   *
   * @param sample the sample index
   * @param output: number of markers
   * @param output: number of haps
   * @param output: hap segs for hap1 (2*sample)
   * @param output: hap segs for hap2 (2*sample + 1)
   */
  void singleStates(int sample, int &outNMarkers, int &outNHaps,
                    QList<HapSegment> &outHapSegs1,
                    QList<HapSegment> &outHapSegs2) const;

  /**
   * Other methods and variables needing access from the SinglePermittedStates object:
   */

  const QList<double>& pos() const { return _pos; }
  double ibdExtend() const { return _ibdExtend; }

  QList< QVector<int> > _hapStates;

private:
  void initPos();

  /**
   * Initializes QList< QVector<int> > _hapStates, whose (j,k)-th
   * element is the edge state at the j-th marker that is traversed
   * by the k-th haplotype.
   * @param dag the DAG constructed by the specified haplotypes
   * @param haps the haplotype pairs used to construct the specified DAG
   */
  void initHapStates();

  void ibsSegs(QList<HapSegment> &outHapSegs, int hap) const;

  /* filter if minimum requirements are not met */
  void containmentFilter(QList<HapSegment> &ibdSegments, int minEndDiff) const;

  const SampleHapPairs &_haps;
  const ImmutableDag &_dag;
  IbsHapSegments _hapSegments;
  QList<double> _pos;
  double _ibdExtend;
};

class Par;
class CurrentData;

/**
 * Class {@code SamplerData} contains immutable input data for the
 * current marker window.
 *
 * "Instances of class {@code SamplerData} are immutable."
 */
class SamplerData
{
 public:

  /**
   * Constructs a new {@code SamplerData} instance from the specified data.
   *
   * @param rdag the {@code RestrictedDag} object
   * @param par the analysis parameters
   * @param cd the input data for the current marker window
   * @param revMarkers {@code true} if and only if the order of markers is reversed
   * // @param runStats the object to which run-time statistics will be written
   */
  SamplerData(const RestrictedDag &rdag, const Par &par, const CurrentData &cd,
              bool revMarkers /* , RunStats runStats */ );
 
  SamplerData& operator=(const SamplerData& other);

  /**
   * Returns {@code true} if the order of markers is reversed, and
   * {@code false} otherwise
   */
  bool markersAreReversed() const { return _revMarkers; }

  /**
   * Returns the number of markers.
   */
  int nMarkers() const { return _gl.nMarkers(); }

  /**
   * Returns the number of samples.
   */
  int nSamples() const { return _gl.nSamples(); }

  /**
   * Returns the number of haplotypes.
   */
  int nHaps() const { return 2*_gl.nSamples(); }

  /**
   * returns the list of markers.
   */
  Markers markers() const { return _gl.markers(); }

  /**
   * Returns the analysis parameters.
   */
  const Par &par() const { return _par; }

  /**
   * Returns the DAG model.
   */
  const RestrictedDag &rdag() const { return _rdag; }

  /**
   * Returns the genotype likelihoods for the
   * target samples at the target data markers.
   */
  const FuzzyGL &gl() const { return _gl; }

  /**
   * Returns the allele error rate parameter
   */
  float err() const;

  /**
   * Returns the probability of recombination between {@code (marker - 1)}
   * and {@code marker}.
   * @param marker a marker index
   */
  float pRecomb(int marker) const { return _recombRate[marker]; }

private:
  void findDagRecombRate(const ImmutableDag &dag, float xdist);

  const Par &_par;
  bool _revMarkers;
  const RestrictedDag &_rdag;
  FuzzyGL _gl;
  QVector<float> _recombRate;
};

#endif
