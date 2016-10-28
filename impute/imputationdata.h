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
  int basePos(int chrom, double geneticPosition){ return (int) (0.5 + geneticPosition / _scaleFactor); }

  /**
   * Returns the genetic map position of the specified marker. The
   * genetic map position is estimated using linear interpolation.
   *
   * @param marker a genetic marker
   */
  double genPos(const Marker &marker) { return _scaleFactor * marker.pos(); }

  /**
   * Returns the genetic map position of the specified genome coordinate.
   * The genetic map position is estimated using linear interpolation.
   *
   * @param chrom the chromosome index
   * @param basePosition the base coordinate on the chromosome
   */
  double genPos(int chrom, int basePosition) { return _scaleFactor * basePosition; }

  /**
   * Returns the scale factor that is multiplied by the chromosome position
   * to obtain the corresponding genetic map position
   * @return the scale factor.
   */
  double scaleFactor() { return _scaleFactor; }

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
  RefHapSeg(SampleHapPairs* refHapPairs, int start, int end);

  /**
   * Returns the reference haplotype pairs.
   */
  SampleHapPairs* refHapPairs() const { return _refHapPairs; }

  /**
   * Return the number of reference allele sequences in this segment.
   */
  int nSeq() { return _seqToHap.length(); }

  /**
   * Return the index of the reference allele sequence in this segment
   * for the specified reference haplotype.
   * @param hap a haplotype index
   */
  int seq(int hap) { return _hapToSeq[hap]; }

  /**
   * Return the specified reference haplotype allele.
   * @param marker index of a marker in this segment
   * @param seq index of a reference allele sequence in this segment
   */
  int allele(int marker, int seq);

  /**
   * Returns the starting marker index (inclusive) of this segment.
   */
  int start() { return _start; }

  /**
   * Returns the ending marker index (exclusive) of this segment.
   */
  int end() { return _end; }

private:
  SampleHapPairs* _refHapPairs;
  int _start;     // inclusive
  int _end;       // exclusive
  QVector<int> _hapToSeq;
  QVector<int> _seqToHap;
};
