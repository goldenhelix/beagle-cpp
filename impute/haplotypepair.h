/* Copyright 2016 Golden Helix, Inc. */
#ifndef HAPLOTYPEPAIR_H
#define HAPLOTYPEPAIR_H

#include "impute/markers.h"
#include "impute/samples.h"

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
  HapPair(Markers markers, Samples samples, int sampleIndex, QList<int> &alleles1,
          QList<int> &alleles2);

  /**
   * (A type of) copy constructor. Can copy in reversed order if specified.
   */
  HapPair(HapPair &other, bool reverse);

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
  QBitArray toBitArray(Markers markers, QList<int> &alleles);
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
  HapPairs(QList<HapPair> hapPairList, bool reverse);

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
   * Returns the number of haplotypes.  The returned value is equal to
   * {@code 2*this.nHapPairs()}.
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
  SampleHapPairs(Samples samples, QList<HapPair> hapPairList, bool reverse);

  /**
   * Returns the samples.  The {@code k}-th sample corresponds to
   * the {@code k}-th haplotype pair.
   */
  Samples samples() const { return _samples; }
  /**
   * Returns the number of samples.
   */
  int nSamples() const { return _samples.nSamples(); }
private:
  void checkSamples();

  Samples _samples;
};

#endif
