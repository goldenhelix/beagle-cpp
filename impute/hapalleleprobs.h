/* Copyright notice.... */
#ifndef HAPALLELEPROBS_H
#define HAPALLELEPROBS_H

#include "impute/haplotypepair.h"

#include <QVector>

/**
 * Class {@code HapAlleleProbs} stores allele probabilities for
 * a haplotype.
 *
 * "All instances of {@code HapAlleleProbs} are required to be immutable."
 */
class HapAlleleProbs
{
public:

  /**
   * Constructs a new {@code HapAlleleProbs} instance.  The
   * {@code alleleProbs} array lists the probability of each allele for
   * each marker, sorted first by marker index, and then by allele index.
   * @param markers the markers
   * @param samples the samples
   * @param hap the haplotype index
   * @param alleleProbs the allele probabilities
   */
  HapAlleleProbs(const Markers &markers, const Samples &samples, int hap,
                 const QVector<float> &alleleProbs);

  /**
   * Returns the specified allele probability.
   *
   * @param marker a marker index
   * @param allele an allele index
   */
  float allele(int marker, int allele) const;

  /**
   * Returns the allele with maximum posterior probability.  If more than
   * one allele has maximum posterior probability, one of the
   * alleles with maximum posterior probability will be returned.
   * @param marker a marker index
   */
  int alleleWithMaxProb(int marker) const;

  /**
   * Returns the number of markers.
   */
  int nMarkers() const {
    return _markers.nMarkers();
  }

  /**
   * Returns the list of markers.
   */
  Markers markers() const {
    return _markers;
  }

  /**
   * Returns the specified marker.
   * @param marker a marker index
   */
  Marker marker(int marker) const {
    return _markers.marker(marker);
  }

  /**
   * Returns the haplotype index. The two haplotypes for sample {@code k}
   * have indices {@code 2*k} and {@code 2*k + 1}.
   */
  int hapIndex() const {
    return _hap;
  }

  /**
   * Returns the list of samples.
   */
  Samples samples() const {
    return _samples;
  }

private:
  float lastAlleleProb(int marker) const;

  Markers _markers;
  Samples _samples;
  int _hap;
  QVector<signed char> _alleleBin;
};

/**
 * Class {@code ConstrainedAlleleProbs} represents per-haplotype
 * allele probabilities for a list of samples. An internal index array
 * determines whether to use the object's list of HapAlleleProbs
 * objects or (instead) to use the object's SampleHapPairs instance to
 * obtain data for any given marker.
 *
 * "All instances of {@code [Constrained]AlleleProbs} are required to
 * be immutable."
 */
class ConstrainedAlleleProbs
{
public:

  /**
   * Construct a new {@code ConstrainedAlleleProbs} instance.  The alleles
   * in the specified {@code SampleHapPairs} object will
   * have probability 1 in the new {@code ConstrainedAlleleProbs} object.
   * All other allele probabilities in the constructed object
   * are determined by the list of {@code HapAlleleProbs} objects.
   * @param shp phased haplotype pairs for a subset of markers
   * @param alProbs the allele probabilities for each haplotype
   * @param indexMap an array of length {@code alProbs[h].nMarkers()}
   * whose {@code j}-th element is the index of marker
   * {@code alProbs[h].marker(j)} in {@code shp.markers()}, or is -1 if the
   * marker is not present in {@code shp.markers()}
   */
  ConstrainedAlleleProbs(SampleHapPairs shp, QList<HapAlleleProbs> alProbs,
                         QList<int> indexMap);

  /**
   * Constructs a new {@code ConstrainedAlleleProbs} instance that
   * merely wraps the specified {@code SampleHapPairs} object.  The
   * alleles in the specified {@code SampleHapPairs} instance will
   * have probability 1. The internal index array will be set so that
   * the {@code SampleHapPairs} object is always used to obtain the
   * data.
   *
   * @param sampleHapPairs the sample haplotype pairs that will
   * be wrapped by {@code this}
   */
  ConstrainedAlleleProbs(SampleHapPairs sampleHapPairs);

  /**
   * Returns {@code true} if the specified marker is not present in
   * the (target) input data and returns {@code false} otherwise.
   * @param marker a marker index
   */
  bool isImputed(int marker) const {
    return (_indexMap[marker] == -1);
  }

  /**
   * Returns the probability that the specified marker allele is
   * present on the first haplotype of the specified sample.
   *
   * @param marker a marker index
   * @param sample a sample index
   * @param allele an allele index
   */
  float alProb1(int marker, int sample, int allele) const;

  /**
   * Returns the probability that the specified marker allele is
   * present on the second haplotype of the specified sample.
   *
   * @param marker a marker index
   * @param sample a sample index
   * @param allele an allele index
   */
  float alProb2(int marker, int sample, int allele) const;

  /**
   * Returns the phased genotype probability, equal to
   * {@code (this.alProb1(marker, sample, allele1)
   *        * this.alProb2(marker, sample, allele2))}.
   *
   * @param marker a marker index
   * @param sample a sample index
   * @param allele1 allele index of the allele on the first haplotype
   * @param allele2 allele index of the allele on the second haplotype
   */
  float gtProb(int marker, int sample, int allele1, int allele2) const {
    return alProb1(marker, sample, allele1)*alProb2(marker, sample, allele2);
  }

  /**
   * Returns the marker allele with maximum probability for the
   * first haplotype of the specified sample. If more than one allele
   * has maximum probability, one of the alleles with maximum
   * probability will be returned.
   * @param marker a marker index
   * @param sample a sample index
   */
  int allele1(int marker, int sample) const;

  /**
   * Returns the marker allele with maximum probability for the
   * second haplotype of the specified sample. If more than one allele
   * has maximum probability, one of the alleles with maximum
   * probability will be returned.
   * @param marker a marker index
   * @param sample a sample index
   */
  int allele2(int marker, int sample) const;

  /**
   * Returns the number of markers.
   */
  int nMarkers() const { return _markers.nMarkers(); }

  /**
   * Returns the specified marker.
   * @param marker a marker index
   */
  Marker marker(int marker) const { return _markers.marker(marker); }

  /**
   * Returns the list of markers.
   */
  Markers markers() const { return _markers; }

  /**
   * Returns the number of samples.
   */
  int nSamples() const { return _samples.nSamples(); }

  /**
   * Returns the list of samples.
   */
  Samples samples() const { return _samples; }

private:
  SampleHapPairs _shp;     // Use _alleleProbs if _indexMap[j] == -1. Otherwise,
  QVector<int> _indexMap;  // use _shp on (target) marker #(_indexMap[j]).

  Markers _markers;
  Samples _samples;
  QList<HapAlleleProbs> _alleleProbs;
};

#endif
