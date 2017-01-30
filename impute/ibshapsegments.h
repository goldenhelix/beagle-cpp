/* Copyright notice. */
#ifndef IBSHAPSEGMENTS_H
#define IBSHAPSEGMENTS_H

#include "impute/haplotypepair.h"

#include <QFile>

#include <QVector>

namespace IbsHapSegUtility
{
  int lowerBoundIndex(const QList<double> &pos, int first, double value);
}


/**
 * Class {@code HapSegment} represents a marker interval for a
 * haplotype.
 *
 * "Instances of class {@code HapSegment} are immutable."
 */
class HapSegment
{
public:

  /**
   * Constructs a new {@code HapSegment} instance.
   * @param hap the haplotype index
   * @param start the start marker index (inclusive)
   * @param end the end marker index (inclusive)
   */
  HapSegment(int hap, int start, int end);

  /**
   * Returns the first haplotype index.
   */
  int hap() const { return _hap; }

  /**
   * Returns the start marker index (inclusive).
   */
  int start() const { return _start; }

  /**
   * Returns the end marker index (inclusive).
   */
  int end() const { return _end; }

  /**
   * Compares this object with the specified object for order.
   * Returns true if this object is less than the specified object.
   *
   * {@code HapSegment} instances are ordered first by {@code
   * this.start()}, then by {@code this.end()}, and finally by
   * {@code this.hap()}.
   *
   * @param other the {@code HapSegment} to be compared
   */

  bool operator< (const HapSegment &other) const;

protected:
  int _hap;
  int _start;
  int _end;
};

/**
 * Class {@code HapSegmentES} is a {@code HapSegment} which orders
 * first by {@code this.end()}, then by {@code this.start()}, and
 * finally by {@code this.hap()}.
 */
class HapSegmentES : public HapSegment
{
public:

  /**
   * Constructs a new {@code HapSegmentES} instance from a {@code
   * HapSegment} instance.  @param hs the HapSegment instance
   */
  HapSegmentES(const HapSegment &hs);

  /**
   * Compares this object with the specified object for order.
   * Returns true if this object is less than the specified object.
   *
   * {@code HapSegmentES} instances are ordered first by {@code
   * this.end()}, then by {@code this.start()}, and finally by
   * {@code this.hap()}.
   *
   * @param other the {@code HapSegmentES} to be compared
   */

  bool operator< (const HapSegmentES &other) const;
};

class IndexMap;

/**
 * Class {@code IbsHapSegments} identifies IBS haplotype segments in
 * a list of sample halotype pairs.
 *
 * "Instances of {@code IbsHapSegments} are immutable."
 *
 *
 * Reference: Gusev A, Lowe JK, Stoffel M, Daly MJ, Altshuler D, Breslow JL,
 *            Friedman JM, Pe'er I (2008) Whole population, genomewide mapping
 *            of hidden relatedness.  Genome Research 2009;19(2):318-26.
 */
class IbsHapSegments
{
public:

  /**
   * Default constructor.
   */
  IbsHapSegments(){}

  /**
   * Initializes a (newly-created) {@code IbsHapSegments} object from the specified data.
   * @param haps the sample haplotype pairs
   * @param pos an array of non-decreasing marker positions whose {@code j}-th
   * element is the position of marker {@code haps.marker(j)}
   * @param minIbsLength the minimum length of a reported IBS segment
   */
  void initialize(const SampleHapPairs &haps, const QList<double> &pos, double minIbsLength);

  /**
   * Returns the sample haplotype pairs.
   */
  // const SampleHapPairs &haps() const { return _haps; }

  /**
   * Returns an array of non-decreasing marker positions whose {@code j}-th
   * element is the position of marker {@code this.haps().marker(j)}.
   */
  QList<double> pos() const { return _pos; }

  /**
   * Returns the minimum length of an IBS segment.
   */
  double minIbsLength() const { return _minIbsLength; }

  /**
   * Returns the list of haplotype segments for other haplotypes that
   * are IBS with the specified haplotype and have length greater
   * than or equal to {@code this.minIbsLength()}.
   *
   * (NOTE: This class is used by ibd code but not for Version 4.1 phasing.)
   *
   * @param outHapSegs receives the list of haplotype segments to be returned
   * @param hap the haplotype index
   */
  void find(QList<HapSegment> &outHapSegs, int hap) const;

  /**
   * Returns a list of haplotype segments for other haplotypes
   * that are IBS with the specified haplotype and that have length greater
   * than or equal to {@code this.minIbsLength()}. An IBS segment is
   * permitted (but not required) to be excluded from the returned
   * list if both end-points of the IBS segment are interior points of
   * another IBS segment.
   *
   * @param outHapSegs receives the list of haplotype segments to be returned
   * @param hap the haplotype index
   */
  void filteredFind(QList<HapSegment> &outHapSegs, int hap) const;

private:
  void checkArguments(const SampleHapPairs &haps, const QList<double> &pos, double minIbsLength) const;
  void findWindowStarts();

  void findIdSets();
  void fillOneIdSet(int ipFirst, int ipSecond);

  void matches(int hap, int window, IndexMap &map) const;

  /* Returns minimum start window index from extended segments */
  int extend(IndexMap &prev, IndexMap &next) const;

  /* (Used by ibd but not for 4.1 phasing....) */
  void save(int hap1, const IndexMap &prev,
            int prevExclEnd, QList<HapSegment> &segments) const;

  void filteredSave(int hap1, const IndexMap &prev,
                    int minExtendedStartWindow,
                    int prevExclEnd, QList<HapSegment> &segments) const;

  int findStart(int hap1, int hap2, int start) const;
  int findInclusiveEnd(int hap1, int hap2, int end) const;

  SampleHapPairs _haps;
  QList<double> _pos;
  double _minIbsLength;
  QList<int> _windowStarts; 
  QList< QVector< QList<int> > > _idSets;
};

class IndexSet;

class RSBDump
{
public:
  RSBDump() : _okToDump(false), _sampleNum(0) {}
  void setOkToDump(bool otd);
  void setNewSample(int sampNum);
  void dumpNext(int rv, int e1, int e2, int i1, int i2);
  void dumpFwd1(int m, int sy1, int sy2, float emp);
  void dumpFwd2(float pr, float r0, float r1, float r2, float ep1, float ep2);
  void dumpFwd3(int pn1, int pn2, float pnp1, float pnp2);
  void cr();
  void dumpSCD(int m, int e1, int e2, int s1, int s2, int node1, int node2, float p1, float p2, float maxSum, float fv, float fs, float gl);
  void dumpIndices(int m, const IndexSet &indices);
  void dumpHSegs(int hap, const QList<HapSegment> &ibsSegs);
  void dumpMSegs(int marker, const QList<HapSegment> &ibsSegs);
  void dumpWs(const QList<int> &ws);
  void dumpIdSets(const QList< QVector< QList<int> > > &idSets);

  bool okToDump() { return _okToDump; }
  void treeHeader(int hap, int start, int end, int size);
  void dumpTSegs(int center, char *type, int depth, const QList<HapSegment> &ibsSegs);
  void dumpModifier(char *svse, int maxStart, int minEnd, const QList<HapSegment> &ibsSegs);

private:
  QFile _out;
  bool _outHasBeenOpened;
  bool _okToDump;
  // int _levelNum;  // Merging also involves the previous level (in addition to later levels).
  int _sampleNum;
};

extern RSBDump globalRsb;

#endif
