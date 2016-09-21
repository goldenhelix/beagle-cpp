/* Copyright 2016 Golden Helix, Inc. */
#ifndef HAPLOTYPEPAIR_H
#define HAPLOTYPEPAIR_H

#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"

#include <QList>

/**
 * Class {@code HapPair} represents a pair of haplotypes for a sample.
 * The pair of haplotypes is guaranteed to have non-missing alleles
 * for each marker.
 */
class HapPair
{
public:
  /**
   * Constructs a new {@code HapPair} instance.
   * @param markers the sequence of markers
   * @param samples the list of samples
   * @param sampleIndex the sample index
   * @param alleles1 the sequence of allele indices for the first haplotype
   * @param alleles2 the sequence of alleles indices for the second haplotype
   */
  HapPair(const Markers &markers, const Samples &samples, int sampleIndex, QList<int> &alleles1,
          QList<int> &alleles2);

  /**
   * (A type of) copy constructor. Can copy in reversed order if specified.
   */
  HapPair(const HapPair &other, bool reverse);

  /**
   * Returns the first allele for the specified marker.
   * @param marker a marker index
   */
  int allele1(int marker) const { return allele(_alleles1Data, marker); }
  /**
   * Returns the second allele for the specified marker.
   * @param marker a marker index
   */
  int allele2(int marker) const { return allele(_alleles2Data, marker); }
  /**
   * Returns the Markers object.
   */
  Markers markers() const { return _markers; }
  /**
   * Returns the specified marker.
   */
  Marker marker(int marker) const { return _markers.marker(marker); }
  /**
   * Returns the number of markers.
   */
  int nMarkers() const { return _markers.nMarkers(); }
  /**
   * Returns the Samples object associated with this haplotype pair.
   */
  Samples samples() const { return _samples; }
  /**
   * Returns the index of the sample associated with this haplotype pair.
   */
  int sampleIndex() const { return _sampleIndex; }
private:
  QBitArray toBitArray(const Markers &markers, QList<int> &alleles);
  int allele(const QBitArray &bitset, int marker) const;

  Markers _markers;
  Samples _samples;
  int _sampleIndex;
  QBitArray _alleles1Data;
  QBitArray _alleles2Data;
};

/**
 * Class {@code HapPairs} represents a list of haplotype pairs.  Each
 * haplotype pair is guaranteed to have two non-missing alleles at
 * each marker.
 *
 * However, in this (base) class, there is no fixed correlation
 * between haplotype pairs and samples, and there might be more than
 * one set of samples represented here.
 */
class HapPairs
{
public:
  /**
   * Constructor with "reverse" flag.
   */
  HapPairs(const QList<HapPair> &hapPairList, bool reverse);

  /**
   * Copy constructor.
   */
  HapPairs(const HapPairs &other);

  HapPairs() : _isReversed(false), _numOfMarkersM1(-1) {}

  /**
  * Returns the allele for the specified marker and haplotype.
  * @param marker a marker index
  * @param haplotype a haplotype index
  */
  int allele(int marker, int haplotype) const;

  /**
   * Returns the first allele for the specified marker and haplotype pair.
   * @param marker a marker index
   * @param hapPair a haplotype pair index
   */
  int allele1(int marker, int hapPair) const;

  /**
   * Returns the second allele for the specified marker and haplotype pair.
   * @param marker a marker index
   * @param hapPair a haplotype pair index
   */
  int allele2(int marker, int hapPair) const;

  /**
   * Returns the number of markers.
   */
  int nMarkers() const { return _numOfMarkersM1 + 1; }
  /**
   * Returns the markers.
   */
  Markers markers() const { return _markers; }
  /**
   * Returns the specified marker.
   * @param marker a marker index
   */
  Marker marker(int mnum) const { return _markers.marker(mnum); }
  /**
   * Returns the number of marker alleles.
   * @param marker a marker index
   */
  int nAlleles(int mnum) const { return _markers.marker(mnum).nAlleles(); }
  /**
   * Returns the number of haplotypes.
   */
  int nHaps() const { return 2 * _hapPairs.length(); }
  /**
   * Returns the number of haplotype pairs.  The returned value is
   * equal to {@code this.nHaps()/2}.
   */
  int nHapPairs() const { return _hapPairs.length(); }
  /**
   * Returns the list of samples containing the sample associated with
   * the specified haplotype pair.
   * @param hapPair a haplotype pair index
   */
  Samples samples(int hapPair) const { return _hapPairs[hapPair].samples(); }
  /**
   * Returns the index of the sample associated with the specified
   * haplotype pair in the list of samples returned by {@code this.samples()}.
   * @param hapPair a haplotype pair index
   */
  int sampleIndex(int hapPair) const { return _hapPairs[hapPair].sampleIndex(); }

protected:
  void checkAndExtractMarkers(bool reverse);

  bool _isReversed;
  int _numOfMarkersM1;
  Markers _markers;
  QList<HapPair> _hapPairs;
};

/**
 * Class {@code SampleHapPairs} represents a list of samples and a
 * haplotype pair for each sample. Each haplotype pair is guaranteed
 * to have two non-missing alleles at each marker.
 *
 * Also, there is a one-to-one correspondence between haplotype pairs
 * and samples represented.
 */
class SampleHapPairs : public HapPairs
{
public:
  SampleHapPairs(const Samples &samples, const QList<HapPair> &hapPairList, bool reverse);

  SampleHapPairs(const SampleHapPairs &other) : HapPairs(other), _samples(other._samples) {}

  SampleHapPairs() : HapPairs() {}

  /**
   * Returns the list of samples containing the sample associated with
   * the specified haplotype pair. (Redundant calling interface needed
   * because of how subclassing works.)
   *
   * @param hapPair a haplotype pair index
   */
  Samples samples(int hapPair) const { return HapPairs::samples(hapPair); }
  /**
   * Returns the samples.  The {@code k}-th sample corresponds to
   * the {@code k}-th haplotype pair.
   */
  Samples samples() const { return _samples; }
  /**
   * Returns the number of samples.
   */
  int nSamples() const { return _samples.nSamples(); }
  /**
   * Returns the number of haplotypes.
   */
  int nHaps() const { return 2 * _samples.nSamples(); }
  /**
   * Returns the number of haplotype pairs.
   */
  int nHapPairs() const { return _samples.nSamples(); }
protected:
  SampleHapPairs(const Samples &samples) : HapPairs(), _samples(samples) {}
  void checkSamples();

  Samples _samples;
};

/**
 * Class {@code RefHapPairs} accesses a list of samples and a
 * haplotype pair for each sample, using BitSetRefGT records.
 */
class RefHapPairs : public SampleHapPairs
{
public:
  /**
   * Constructs a new {@code RefHapPairs} instance. A Markers object
   * is created from the data and is available for copying (through
   * the appropriate base class method).
   *
   * @param samples the sequence of samples
   * @param refVcfRecs the sequence of phased, non-missing per-marker
   * genotype data
   */
  RefHapPairs(const Samples &samples, const QList<BitSetRefGT> &refVcfRecs);

  RefHapPairs() : SampleHapPairs() {}

  int allele1(int marker, int hapPair) const { return _refVcfRecs[marker].allele1(hapPair); }
  int allele2(int marker, int hapPair) const { return _refVcfRecs[marker].allele2(hapPair); }
  int allele(int marker, int haplotype) const;

  Samples samples(int hapPair) const;
  Samples samples() const { return _samples; }
  int sampleIndex(int hapPair) const;

  int nAlleles(int marker) const { return _refVcfRecs[marker].nAlleles(); }
private:
  void createMarkers();

  QList<BitSetRefGT> _refVcfRecs;
};

/**
 * Class {@code GLSampleHapPairs} ("GL" means Genotype Likelihoods),
 * in addition to representing a list of samples and a haplotype pair
 * for each sample, can output a genotype likelihood for any given
 * sample. This (base-class) member of the "GL" family is always
 * categorized as "reference data".
 */
class GLSampleHapPairs : public SampleHapPairs
{
public:
  /**
   * Constructs a new {@code GLSampleHapPairs} instance from the
   * specified data. Will always be in forward direction.
   *
   * @param another instance of a "GL"-family class. In default mode,
   * must only contain phased, non-missing genotype data.
   * @param optional parameter to be used only by the SplicedGL
   * copy-and-maybe-reverse constructor.
   * @param optional parameter to be used only by the SplicedGL
   * copy-and-maybe-reverse constructor.
   */
  GLSampleHapPairs(const GLSampleHapPairs &otherGL, bool checkRef = true, bool reverse = false);

  GLSampleHapPairs() : SampleHapPairs() {}

  /**
   * Returns {@code true} if the observed data for each marker and
   * sample includes a phased genotype that has no missing alleles,
   * and returns {@code false} otherwise.
   */
  virtual bool isRefData() const { return true; }

  int allele(int marker, int haplotype) const
  {
    int sample = haplotype / 2;
    if ((haplotype & 1) == 0) {
      return allele1(marker, sample);
    } else {
      return allele2(marker, sample);
    }
  }

  int allele1(int marker, int hapPair) const;

  int allele2(int marker, int hapPair) const;

  Samples samples(int hapPair) const;
  Samples samples() const { return _samples; }
  int sampleIndex(int hapPair) const;

  /**
   * Returns the number of markers (overall).
   */
  int nMarkers() const { return _numOfGlMarkersM1 + 1; }
  /**
   * Returns the markers object (for all of the markers).
   */
  Markers markers() const { return _glMarkers; }
  /**
   * Returns the specified marker. (Can be any marker.)
   * @param marker a marker index
   */
  Marker marker(int mnum) const { return _glMarkers.marker(mnum); }
  /**
   * Returns the number of marker alleles. (Can be of any marker.)
   * @param marker a marker index
   */
  int nAlleles(int mnum) const { return _glMarkers.marker(mnum).nAlleles(); }
  /**
   * Returns the number of haplotypes in the marker. (Can be of any marker.)
   */

protected:
  // Partial constructor working on behalf of the SplicedGL(samples,
  // vma) constructor.
  GLSampleHapPairs(const Samples &samples);

  // Constructor working on behalf of the SplicedGL(haps, otherGL)
  // constructor.
  GLSampleHapPairs(const SampleHapPairs &haps, const GLSampleHapPairs &otherGL);

  int _overlap;
  bool _glIsReversed;
  QList<BitSetGT> _vcfRecs;
  int _numOfGlMarkersM1;
  Markers _glMarkers;
};

class SplicedGL : public GLSampleHapPairs
{
public:
  /**
   * Constructs a new {@code SplicedGL} instance from the specified
   * data. The resulting object will always have data in the forward
   * direction, and there will be no overlapping haplotypes embedded
   * in the resulting object. Also, the genotype data may or may not
   * be phased.
   *
   * This constructor is equivalent to a "BasicGL" constructor.
   *
   * @param A Samples object.
   * @param A QList of BitSetGT objects.
   */
  SplicedGL(const Samples &samples, const QList<BitSetGT> &vma);

  /**
   * Constructs a new {@code SplicedGL} instance from the specified
   * data. The resulting object will always have data in the forward
   * direction. Overlapping haplotypes will presumably be embedded
   * in the resulting object. Also, the genotype data may or may not
   * be phased.
   *
   * @param A SampleHapPairs object.
   * @param another instance of a "GL"-family class, for which it is
   * OK to contain unphased and/or missing genotype data.
   */
  SplicedGL(const SampleHapPairs &haps, const GLSampleHapPairs &otherGL);

  /**
   * Copy-and-possibly-reverse {@code SplicedGL} constructor.
   *
   * @param Another instance of a "GL"-family class.
   * @param Whether or not to reverse the data.
   */
  SplicedGL(const GLSampleHapPairs &otherGL, bool reverse);

  SplicedGL() : GLSampleHapPairs(), _isRefData(false) {}

  /**
   * Returns {@code true} if the observed data for each marker and
   * sample includes a phased genotype that has no missing alleles,
   * and returns {@code false} otherwise.
   */
  bool isRefData() const { return _isRefData; }
  /**
   * Returns the probability of the observed data for the specified marker
   * and sample if the specified pair of ordered alleles is the true
   * ordered genotype.
   * @param marker the marker index
   * @param sample the sample index
   * @param allele1 the first allele index
   * @param allele2 the second allele index
   */
  double gl(int marker, int sample, int allele1, int allele2) const;

  /**
   * Returns {@code true} if the observed data for the specified
   * marker and sample includes a phased genotype, and returns {@code false}
   * otherwise.
   * @param marker the marker index
   * @param sample the sample index
   */
  bool isPhased(int marker, int sample) const;

protected:
  void checkSamples();
  void createMarkers();
  void checkIfIsRefData();

  bool _isRefData;
};

/**
 * Class {@code FuzzyGL} is a {@code GL} class that
 * incorporates a fixed error rate for the
 * observed (emitted) allele to differ from the true allele.  Allele
 * errors are independent.
 */

class FuzzyGL : public SplicedGL
{
public:
  /**
   * Constructs a {@code FuzzyGL} instance.
   * @param gl the genotype likelihoods without error
   * @param err the allele error rate
   */
  FuzzyGL(const SplicedGL &gl, double err, bool reverse);

  FuzzyGL() : SplicedGL() {}

  double gl(int marker, int sample, int a1, int a2);

private:
  double phasedGL(int obs1, int obs2, int a1, int a2);

  double _ee;
  double _ef;
  double _ff;
};

#endif
