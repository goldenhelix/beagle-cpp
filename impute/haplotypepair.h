/* Copyright 2016 Golden Helix, Inc. */
#ifndef HAPLOTYPEPAIR_H
#define HAPLOTYPEPAIR_H

#include "impute/samples.h"
#include "impute/markers.h"

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
    HapPair(Markers markers, Samples samples, int sampleIndex,
            QList<int> &alleles1, QList<int> &alleles2);

    /**
     * (A type of) copy constructor. Can copy in reversed order if specified.
     */
    HapPair(HapPair &other, bool reverse);

   /**
    * Returns the first allele for the specified marker.
    * @param marker a marker index
    */
    int allele1(int marker) { return allele(_alleles1Data, marker); }

   /**
    * Returns the second allele for the specified marker.
    * @param marker a marker index
    */
    int allele2(int marker) { return allele(_alleles2Data, marker); }

    /**
     * Returns the Markers object.
     */
    Markers markers() { return _markers; }

    /**
     * Returns the specified marker.
     */
    Marker marker(int marker) { return _markers.marker(marker); }

    /**
     * Returns the number of markers.
     */
    int nMarkers() { return _markers.nMarkers(); }

    /**
     * Returns the Samples object associated with this haplotype pair.
     */
    Samples samples() { return _samples; }

    /**
     * Returns the index of the sample associated with this haplotype pair.
     */
    int sampleIndex() { return _sampleIndex; }

private:
  QBitArray toBitArray(Markers markers, QList<int> &alleles);
  int allele(QBitArray &bitset, int marker);

  Markers _markers;
  Samples _samples;
  int _sampleIndex;
  QBitArray _alleles1Data;
  QBitArray _alleles2Data;
};
#endif
