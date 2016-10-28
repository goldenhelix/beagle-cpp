/* Copyright notice.... */
#ifndef IMPUTATIONDATA_H
#define IMPUTATIONDATA_H

#include "impute/hapalleleprobs.h"

#include <QVector>

#endif

/**
 * Class {@code PositionMap} represents a genetic map obtained by
 * multiplying chromosome position by a scale factor.
 * 
 * "Instances of class {@code PositionMap} are immutable."
 */
class PositionMap
{
public:

  /**
   * Constructs a new {@code PositionMap} instance.
   * @param scaleFactor the factor that is multiplied by
   * a base position to obtain the corresponding genetic map position
   */
 PositionMap(double scaleFactor) : _scaleFactor(scaleFactor) {}

  /**
   * Returns the base position corresponding to the specified genetic map
   * position. If the genetic position is not a map position then the base
   * position is estimated from the nearest genetic map positions using
   * linear interpolation.
   *
   * @param chrom the chromosome index
   * @param geneticPosition the genetic position on the chromosome
   */
  int basePos(int chrom, double geneticPosition) const { return (int) (0.5 + geneticPosition / _scaleFactor); }

  /**
   * Returns the genetic map position of the specified marker. The
   * genetic map position is estimated using linear interpolation.
   *
   * @param marker a genetic marker
   */
  double genPos(const Marker &marker) const { return _scaleFactor * marker.pos(); }

  /**
   * Returns the genetic map position of the specified genome coordinate.
   * The genetic map position is estimated using linear interpolation.
   *
   * @param chrom the chromosome index
   * @param basePosition the base coordinate on the chromosome
   */
  double genPos(int chrom, int basePosition) const { return _scaleFactor * basePosition; }

  /**
   * Returns the scale factor that is multiplied by the chromosome position
   * to obtain the corresponding genetic map position
   * @return the scale factor.
   */
  double scaleFactor() const { return _scaleFactor; }

private:
  double _scaleFactor;
};

/**
 * Class {@code RefHapSeg} represents a chromosome segment of
 * reference haplotypes.
 *
 * "Instances of class {@code RefHapSeg} are immutable."
 */
class RefHapSeg
{
public:

  /**
   * Constructs a new {@code RefHapSegs} instance from the specified data.
   * @param refHapPairs the reference haplotype pairs
   * @param start the starting marker index (inclusive)
   * @param end the ending marker index (exclusive)
   */
  RefHapSeg(const SampleHapPairs &refHapPairs, int start, int end);

  /**
   * Returns the reference haplotype pairs.
   */
  const SampleHapPairs& refHapPairs() const { return _refHapPairs; }

  /**
   * Return the number of reference allele sequences in this segment.
   */
  int nSeq() const { return _seqToHap.length(); }

  /**
   * Return the index of the reference allele sequence in this segment
   * for the specified reference haplotype.
   * @param hap a haplotype index
   */
  int seq(int hap) const { return _hapToSeq[hap]; }

  /**
   * Return the specified reference haplotype allele.
   * @param marker index of a marker in this segment
   * @param seq index of a reference allele sequence in this segment
   */
  int allele(int marker, int seq) const;

  /**
   * Returns the starting marker index (inclusive) of this segment.
   */
  int start() const { return _start; }

  /**
   * Returns the ending marker index (exclusive) of this segment.
   */
  int end() const { return _end; }

private:
  SampleHapPairs _refHapPairs;
  int _start;     // inclusive
  int _end;       // exclusive
  QVector<int> _hapToSeq;
  QVector<int> _seqToHap;
};


/**
 * Class {@code RefHapSegs} represents reference haplotypes that span
 * segments determined by non-overlapping clusters of markers.
 *
 * "Instances of class {@code RefHapSegs} are immutable."
 */
class RefHapSegs
{
public:

  /**
   * Constructs a new {@code RefHapSegs} instance from the specified data.
   * @param refHapPairs the reference haplotype pairs
   * @param segStart an array whose {@code j}-th element is the
   * starting reference marker index (inclusive) for the {@code j}-th
   * segment of markers
   * @param segEnd an array whose {@code j}-th element is the
   * ending reference marker index (exclusive) for the {@code j}-th
   * segment of markers
   */
  RefHapSegs(const SampleHapPairs &refHapPairs, const QVector<int> &segStart, const QVector<int> &segEnd);

  /**
   * Returns the reference haplotype pairs.
   */
  const SampleHapPairs& refHapPairs() const { return _refHapPairs; }

  /**
   * Return the number of distinct reference allele sequences in the
   * specified chromosome segment.
   * @param segment index of a chromosome segment determined by
   * the marker clusters
   */
  int nSeq(int segment) const { return _refHapSegs[segment].nSeq(); }

  /**
   * Return the number of markers in the specified chromosome segment.
   * @param segment index of a chromosome segment determined by
   * the marker clusters
   */
  int nMarkers(int segment) const { return _refHapSegs[segment].end() - _refHapSegs[segment].start(); }

  /**
   * Return the index of the allele sequence in the specified chromosome
   * segment for the specified reference haplotype.
   *
   * @param segment index of a chromosome segment determined by
   * the marker clusters
   * @param hap a haplotype index
   */
  int seq(int segment, int hap) const { return _refHapSegs[segment].seq(hap); }

  /**
   * Return the specified reference haplotype allele.
   *
   * @param segment index of a chromosome segment determined by
   * the marker clusters
   * @param marker index of a marker in the specified interval
   * @param seq index of a reference allele sequence in the specified
   * interval
   */
  int allele(int segment, int marker, int seq) const { return _refHapSegs[segment].allele(marker, seq); }

  /**
   * Returns the number of marker segments.
   */
  int nSegs() const { return _segStart.length(); }

  /**
   * Returns the index of the first marker (inclusive) in the specified
   * marker segment.
   * @param segment an index of a marker segment
   */
  int segStart(int segment) const { return _segStart[segment]; }

  /**
   * Returns the index of the last marker (exclusive) in the specified
   * marker segment.
   * @param segment an index of a marker segment
   */
  int segEnd(int segment) const { return _segEnd[segment]; }

private:
  void checkClusters(const QVector<int> &starts, const QVector<int> &ends, int nMarkers);

  QVector<int> _segStart;
  QVector<int> _segEnd;
  SampleHapPairs _refHapPairs;
  QList<RefHapSeg> _refHapSegs;
};

