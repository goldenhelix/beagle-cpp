#include "impute/imputationdata.h"

#include <math.h>

#define MIN_CM_DIST     1e-7

class HapSegData
{
public:

  HapSegData(const RefHapPairs &refHapPairs, int start, int end);

  QVector<int> hap2Seq() const { return _hap2seq; }

  QVector<int> seq2Hap() const;

private:
  void highMafUpdate(const RefHapPairs &refHapPairs, int marker);
  int indexOfAllele(const QList<int> &list, int allele) const;
  void mapAlleleToNewSeq(int seq, int allele);
  void updateHap2Seq(int h, int seq);

  QVector<int> _hap2seq;
  QList<int> _seq2Cnt;

  QList<int> _emptyList;
  QList< QList<int> > _seq2AlleleMap;
};

/**
 * Class {@code HaplotypeCoder} indexes the observed allele sequences
 * in reference and target haplotype pairs for a list of consecutive markers.
 *
 * "Instances of class {@code HaplotypeCoder} are immutable."
 */
class HaplotypeCoder
{
public:

  /**
   * Constructs a new {@code HaplotypeCoder} instance from the specified
   * data.
   * @param refHapPairs the reference haplotype pairs
   * @param targetHapPairs the target haplotype pairs
   */
  HaplotypeCoder(const SampleHapPairs &refHapPairs, const SampleHapPairs &targetHapPairs);

  /**
   * Returns the reference haplotype pairs used to construct this.
   */
  SampleHapPairs refHapPairs() { return _refHapPairs; }

  /**
   * Returns the target haplotype pairs used to construct this.
   */
  SampleHapPairs targetHapPairs() { return _targetHapPairs; }

  /**
   * Returns (sets) two lists--the first list maps each reference
   * haplotype index to the index of the allele sequence carried by that
   * reference haplotype, and the second list maps each target haplotype
   * index to the index of the allele sequence carried by that target
   * haplotype. The size of the first list is
   * {@code this.refHapPairs().nHaps()}, and the size of the second
   * list is {@code this.targetHapPairs().nHaps()}.
   *
   * @param refMap the first "returned" list containing the haplotype indices for the reference samples.
   * @param targMap the second "returned" list containing the haplotype indices for the target samples.
   * @param start the first marker index (inclusive)
   * @param end the last marker index (exclusive)
   */
  void findHaplotypes(QList<quint16> &refMap, QList<quint16> &targMap, int start, int end);

private:
  QList<int> partition(int marker, QList<int> lastEnds);
  void setAlleles(int marker);
  int partition(int start, int end, int splitAllele);
  void setAllelesToHapIndices(QList<int> ends);

  int _nRefHaps;
  int _nHaps;
  QVector<int> _alleles;
  QVector<int> _haps;
  SampleHapPairs _refHapPairs;
  SampleHapPairs _targetHapPairs;
};

HapSegData::HapSegData(const RefHapPairs &refHapPairs, int start, int end)
{
  _hap2seq.fill(0, refHapPairs.nHaps());
  _seq2Cnt.append(_hap2seq.length());
  _seq2AlleleMap.append(_emptyList);

  for (int m=start; m<end; ++m)
    highMafUpdate(refHapPairs, m);
}

QVector<int> HapSegData::seq2Hap() const
{
  QVector<int> seqToHap(_seq2AlleleMap.size(), -1);

  for (int h=0; h < _hap2seq.length(); ++h)
  {
    int seq = _hap2seq[h];

    if (seqToHap[seq] == -1)
      seqToHap[seq] = h;
  }

  return seqToHap;
}

void HapSegData::highMafUpdate(const RefHapPairs &refHapPairs, int marker)
{
  for (int j=0, n = _seq2AlleleMap.size(); j<n; ++j)
    _seq2AlleleMap[j].clear();

  for (int h=0; h < _hap2seq.length(); ++h)
  {
    int seq = _hap2seq[h];

    int allele = refHapPairs.allele(marker, h);

    if (_seq2AlleleMap[seq].isEmpty())
    {
      _seq2AlleleMap[seq].append(allele);
      _seq2AlleleMap[seq].append(seq);
    }
    else
    {
      int index = indexOfAllele(_seq2AlleleMap[seq], allele);

      if (index==_seq2AlleleMap[seq].size())
        mapAlleleToNewSeq(seq, allele);

      updateHap2Seq(h, _seq2AlleleMap[seq][index+1]);
    }
  }
}

int HapSegData::indexOfAllele(const QList<int> &list, int allele) const
{
  int index=0;

  while (index < list.size()  &&  list[index] != allele)
    index += 2;

  return index;
}

void HapSegData::mapAlleleToNewSeq(int seq, int allele)
{
  _seq2AlleleMap[seq].append(allele);
  _seq2AlleleMap[seq].append(_seq2AlleleMap.size());
  _seq2AlleleMap.append(_emptyList);
  _seq2Cnt.append(0);
}

void HapSegData::updateHap2Seq(int h, int seq)
{
  if(h >= _hap2seq.length())
    int a = h;
  if(_hap2seq[h] >= _seq2Cnt.length())
    int c = _hap2seq[h];
  if(seq >= _seq2Cnt.length())
    int e = seq;

  --_seq2Cnt[_hap2seq[h]];
  _hap2seq[h] = seq;
  ++_seq2Cnt[seq];
}

RefHapSeg::RefHapSeg(const RefHapPairs &refHapPairs, int start, int end)
{
  Q_ASSERT_X(start >= 0  &&  start < end  &&  end <= refHapPairs.nMarkers(),
             "RefHapSeg::RefHapSeg",
             "invalid start and/or end values");

  HapSegData hapSegData(refHapPairs, start, end);
  _refHapPairs = refHapPairs;
  _start = start;
  _end = end;
  _hapToSeq = hapSegData.hap2Seq();
  _seqToHap = hapSegData.seq2Hap();
}

int RefHapSeg::allele(int marker, int seq) const
{
  int refIndex = _start + marker;
  Q_ASSERT_X(marker >= 0  &&  refIndex < _end,
             "RefHapSeg::allele",
             "marker < 0 || refIndex >= _end");

  return _refHapPairs.allele(refIndex, _seqToHap[seq]);
}

RefHapSegs::RefHapSegs(const RefHapPairs &refHapPairs, const QVector<int> &segStart, const QVector<int> &segEnd)
{
  int nMarkers = refHapPairs.nMarkers();
  checkClusters(segStart, segEnd, nMarkers);
  _segStart = segStart;
  _segEnd = segEnd;
  _refHapPairs = refHapPairs;

  for(int j=0, n=segStart.length(); j<n; j++)
  {
    RefHapSeg rhs(refHapPairs, segStart[j], segEnd[j]);
    _refHapSegs.append(rhs);
  }
}

void RefHapSegs::checkClusters(const QVector<int> &starts, const QVector<int> &ends, int nMarkers)
{
  Q_ASSERT_X(starts.length() == ends.length(),
             "RefHapSegs::checkClusters",
             "starts.length() != ends.length()");

  for (int j=0; j < starts.length(); ++j)
  {
    Q_ASSERT_X(starts[j] >= 0  &&  starts[j] < ends[j]  &&  ends[j] <= nMarkers,
               "RefHapSegs::checkClusters",
               "starts[j] < 0 || starts[j] >= ends[j] || ends[j] > nMarkers");
  }
}

HaplotypeCoder::HaplotypeCoder(const SampleHapPairs &refHapPairs, const SampleHapPairs &targetHapPairs)
{
  Q_ASSERT_X(refHapPairs.markers() == targetHapPairs.markers(),
             "HaplotypeCoder::HaplotypeCoder",
             "inconsistent markers");

  _nRefHaps = refHapPairs.nHaps();
  _nHaps = _nRefHaps + targetHapPairs.nHaps();
  _refHapPairs = refHapPairs;
  _targetHapPairs = targetHapPairs;
}

void HaplotypeCoder::findHaplotypes(QList<quint16> &refMap, QList<quint16> &targMap,
                                    int start, int end)
{
  Q_ASSERT_X(start < end,
             "HaplotypeCoder::findHaplotypes",
             "start >= end");

  _alleles.fill(0, _nHaps);

  _haps.fill(0, _nHaps);
  for (int j=0; j < _nHaps; j++)
    _haps[j] = j;

  QList<int> lastEnds;
  lastEnds.append(_nHaps);

  for (int m=start; m<end; ++m)
    lastEnds = partition(m, lastEnds);

  setAllelesToHapIndices(lastEnds);

  for (int j = 0; j < _nRefHaps; j++)
    refMap.append(_alleles[j]);
  for (int j = _nRefHaps; j < _nHaps; j++)
    targMap.append(_alleles[j]);
}

QList<int> HaplotypeCoder::partition(int marker, QList<int> lastEnds)
{
  QList<int> nextEnds;
  nextEnds.reserve((4*lastEnds.size())/3 + 1 );
  setAlleles(marker);
  int nAlleles = _refHapPairs.marker(marker).nAlleles();
  int lastAllele = nAlleles - 1;

  int start = 0;
  for (int j=0, n=lastEnds.size(); j<n; ++j)
  {
    int end = lastEnds[j];
    for (int al=0; al<lastAllele; ++al)
    {
      int nextStart = partition(start, end, al);
      if (nextStart > start)
      {
        nextEnds.append(nextStart);
        start = nextStart;
      }
    }
    if (end > start)
      nextEnds.append(end);

    start = end;
  }

  return nextEnds;
}

void HaplotypeCoder::setAlleles(int marker)
{
  for (int j=0; j < _nRefHaps; ++j)
    _alleles[j] = _refHapPairs.allele(marker, j);

  for (int j = _nRefHaps; j < _nHaps; ++j)
    _alleles[j] = _targetHapPairs.allele(marker, j - _nRefHaps);
}

/* Returns the start index of the second partitioned set */
int HaplotypeCoder::partition(int start, int end, int splitAllele)
{
  int nextStart = end;
  while (start < nextStart)
  {
    int allele = _alleles[_haps[start]];
    if (allele == splitAllele)
      ++start;

    else
    {
      --nextStart;
      int tmp = _haps[nextStart];
      _haps[nextStart] = _haps[start];
      _haps[start] = tmp;
    }
  }

  return nextStart;
}

void HaplotypeCoder::setAllelesToHapIndices(QList<int> ends)
{
  int start = 0;
  for (int j=0, n=ends.size(); j<n; ++j)
  {
    int end = ends[j];
    for (int k=start; k<end; ++k)
      _alleles[_haps[k]] = j;

    start = end;
  }
}

// Helper functions for ImputationData....

static QList<int> findTargClustEnd(Markers targetMarkers, PositionMap genMap,
                                   float clusterDist)
{
  int nMarkers = targetMarkers.nMarkers();
  QList<int> ends;
  double startPos = genMap.genPos(targetMarkers.marker(0));
  int index = 0;
  for (int m=1; m<nMarkers; ++m)
  {
    double pos = genMap.genPos(targetMarkers.marker(m));
    if ((pos - startPos) > clusterDist)
    {
      ends.append(m);
      startPos = pos;
    }
  }
  ends.append(nMarkers);
  return ends;
}

static void setCodedAlleles(const SampleHapPairs &refHapPairs,
			    const SampleHapPairs &targetHapPairs, const QList<int> targEnd,
            QList< QList<quint16> > &refAlleles, QList< QList<quint16> > &targAlleles)
{
  HaplotypeCoder coder(refHapPairs, targetHapPairs);
  int start = 0;
  for (int j=0; j<targEnd.length(); ++j)
  {
    QList<quint16> refMap;
    QList<quint16> targMap;
    coder.findHaplotypes(refMap, targMap, start, targEnd[j]);
    refAlleles.append(refMap);
    targAlleles.append(targMap);
    start = targEnd[j];
  }
}

static QVector<float> err(float errRate, QList<int> targClustEnd)
{
  float maxErrProb = 0.5f;
  QVector<float> err(targClustEnd.length());
  int start = 0;
  for (int j=0; j<err.length(); ++j)
  {
    err[j] = errRate * (targClustEnd[j] - start);

    if (err[j] > maxErrProb)
      err[j] = maxErrProb;

    start = targClustEnd[j];
  }

  return err;
}

static double maxd(double a, double b)
{
  return (a > b) ? a : b;
}

static QVector<float> findPRcomb(int chrom, QVector<int> midPos, int nHaps,
                                 PositionMap map, float ne)
{
  QVector<float> rr(midPos.length());
  double c = -(0.04*ne/nHaps);    // 0.04 = 4/(100 cM/M)
  double lastGenPos = map.genPos(chrom, midPos[0]);
  rr[0] = 0.0;
  for (int j=1; j<rr.length(); ++j)
  {
    double genPos = map.genPos(chrom, midPos[j]);
    double genDist = maxd(abs(genPos - lastGenPos), MIN_CM_DIST);
    rr[j] = -expm1f(c*genDist);
    lastGenPos = genPos;
  }
  return rr;
}

static QVector<int> findMidPos(const Markers &refMarkers, const RefHapSegs &refHapSegs)
{
  QVector<int> midPos(refHapSegs.nSegs() - 1);
  for (int j=0; j<midPos.length(); ++j)
  {
    int startPos = refMarkers.marker(refHapSegs.segStart(j+1)).pos();
    int endPos = refMarkers.marker(refHapSegs.segEnd(j) - 1).pos();
    midPos[j] = (startPos + endPos) / 2;
  }
  return midPos;
}

static QVector<float> findPRecomb(const RefHapSegs &refHapSegs, const PositionMap &map, float ne)
{
  const SampleHapPairs &refHaps = refHapSegs.refHapPairs();
  const Markers &refMarkers = refHaps.markers();
  int nHaps = refHaps.nHaps();
  QVector<int> midPos = findMidPos(refMarkers, refHapSegs);
  int chrom = refMarkers.marker(0).chromIndex();
  return findPRcomb(chrom, midPos, nHaps, map, ne);
}

static QVector<double> findCumPos(const Markers &markers, const PositionMap &map)
{
  QVector<double> cumPos(markers.nMarkers());
  double lastGenPos = map.genPos(markers.marker(0));
  cumPos[0] = 0.0;
  for (int j=1; j<cumPos.length(); ++j)
  {
    double genPos = map.genPos(markers.marker(j));
    double genDist = maxd(abs(genPos - lastGenPos), MIN_CM_DIST);
    cumPos[j] = cumPos[j-1] + genDist;
    lastGenPos = genPos;
  }
  return cumPos;
}

static QVector<float> findWts(const RefHapSegs &refHapSegs, const PositionMap &map)
{
  const Markers &refMarkers = refHapSegs.refHapPairs().markers();
  QVector<double> cumPos = findCumPos(refMarkers, map);
  int nMarkers = refMarkers.nMarkers();
  int nClusters = refHapSegs.nSegs() - 1;
  QVector<float> wts(cumPos.length(), 0.0);

  for (int j=1; j<nClusters; ++j)
  {
    int start = refHapSegs.segStart(j);
    int end = refHapSegs.segEnd(j-1);
    int nextStart = refHapSegs.segStart(j+1);
    double nextStartPos = cumPos[nextStart];
    double totalLength = nextStartPos - cumPos[end - 1];
    for (int m=start; m<end; ++m)
      wts[m] = 1.0;
    for (int m=end; m<nextStart; ++m)
      wts[m] = (float) ( (nextStartPos - cumPos[m]) / totalLength );
  }

  return wts;
}

static RefHapSegs findRefHapSegs(RefHapPairs refHapPairs,
                             QList<int> targClustEnd, QList<int> targToRef)
{
  int n = targClustEnd.length();
  QVector<int> refClusterStart(n + 1, 0);
  QVector<int> refClusterEnd(n + 1, 0);

  refClusterStart[0] = 0;
  refClusterEnd[0] = targToRef[targClustEnd[0] - 1] + 1;

  refClusterStart[1] = targToRef[0];
  refClusterEnd[1] = targToRef[targClustEnd[1] - 1] + 1;

  for (int j=2; j<n; ++j)
  {
    refClusterStart[j] = targToRef[targClustEnd[j-2]];
    refClusterEnd[j] = targToRef[targClustEnd[j] - 1] + 1;
  }

  refClusterStart[n] = targToRef[targClustEnd[n-2]];
  refClusterEnd[n] = refHapPairs.nMarkers();

  return RefHapSegs(refHapPairs, refClusterStart, refClusterEnd);
}

ImputationData::ImputationData(const Par &par, const CurrentData &cd,
                               const SampleHapPairs &targetHapPairs,
                               const PositionMap &map)
{
  Q_ASSERT_X(cd.targetMarkers() == targetHapPairs.markers(),
             "ImputationData::ImputationData",
             "inconsistent markers");

  Q_ASSERT_X(cd.targetSamples() == targetHapPairs.samples(),
             "ImputationData::ImputationData",
             "inconsistent samples");

  QList<int> targClustEnd = findTargClustEnd(targetHapPairs.markers(), map, par.cluster());

  setCodedAlleles(cd.restrictedRefSampleHapPairs(), targetHapPairs,
                  targClustEnd, _refAlleles, _targAlleles);

  _nClusters = targClustEnd.length();
  _refHapPairs = cd.refSampleHapPairs();
  _refHapSegs = findRefHapSegs(_refHapPairs, targClustEnd, cd.markerIndices());
  _targHapPairs = targetHapPairs;
  _errProb = err(par.err(), targClustEnd);
  _pRecomb = findPRecomb(_refHapSegs, map, par.ne());
  _weight = findWts(_refHapSegs, map);
}

// Utilities for LSHapBaum....

static float sum(const QVector<float> &fa)
{
  float sum = 0.0;
  foreach (float f, fa)
    sum += f;

  return sum;
}

static void scale(QVector<float> &fa, float divisor)
{
  for (int j=0; j<fa.length(); ++j)
    fa[j] /= divisor;
}

static float findThreshold(int nSeq)
{
  float tentThresh = 1.0 / nSeq;
  return (tentThresh < 0.005) ? tentThresh : 0.005;
}

LSHapBaum::LSHapBaum(const ImputationData &impData, bool lowMem)
{
  _windowIndex = -9999;
  _arrayIndex = -9999;
  _impData = impData;
  _lowMem = lowMem;
  _n = impData.refHapPairs().nHaps();
  _refMarkers = impData.refHapPairs().markers();
  _alleleProbs.fill(0.0, _refMarkers.sumAlleles());

  int nClusters = impData.nClusters();
  int size = lowMem ? (int) ceil(sqrt(1 + 8*nClusters)/2.0) + 1
                    : nClusters;
  _fwdValueIndex2Marker.fill(0, size);

  QVector<float> zeros(_n, 0.0);
  for (int s=0; s < size; s++)
    _fwdVal.append(zeros);

  _bwdVal.fill(0.0, _n);

  _refHapSegs = impData.refHapSegs();

  for (int j=0; j < nClusters; ++j)
  {
    QVector<float> fwdp(_refHapSegs.nSeq(j+1), 0.0);
    QVector<float> bwdp(_refHapSegs.nSeq(j));
    _fwdHapProbs.append(fwdp);
    _bwdHapProbs.append(bwdp);
  }
}

HapAlleleProbs LSHapBaum::randomHapSample(int hap)
{
  _alleleProbs.fill(0.0);
  int nMarkers = _impData.nClusters();
  _windowIndex = 0;
  _arrayIndex = -1;
  setForwardValues(0, nMarkers, hap);
  _bwdVal.fill(1.0/_n);
  setStateProbs(nMarkers-1, currentIndex());
  for (int m=nMarkers-2; m>=0; --m)
  {
    setBwdValue(m, hap);
    setStateProbs(m, previousIndex(hap));
  }
  setAlleleProbs();
  return HapAlleleProbs(_refMarkers, _impData.targetSamples(),
                                  hap, _alleleProbs);
}

void LSHapBaum::setForwardValues(int start, int end, int hap)
{
  float lastSum = 1.0;
  for (int m=start; m<end; ++m)
  {
    float probRec = _impData.pRecomb(m);
    float probNoRec = 1.0f - probRec;
    float noErrProb = _impData.noErrProb(m);
    float errProb = _impData.errProb(m);
    float shift = probRec/_n;
    float scale = probNoRec/lastSum;
    int prev = currentIndex();
    int next = nextIndex();
    float sum = 0.0;
    _fwdValueIndex2Marker[next] = m;
    int a = _impData.targetAllele(m, hap);
    for (int h=0; h < _n; ++h)
    {
      float em = (a == _impData.refAllele(m, h)) ? noErrProb : errProb;
      
	  _fwdVal[next][h] = (m==0) ? em : em*(scale*_fwdVal[prev][h] + shift);
      sum += _fwdVal[next][h];
    }
    lastSum = sum;
  }
}

void LSHapBaum::setBwdValue(int m, int hap)
{
  int mP1 = m + 1;
  float probRec = _impData.pRecomb(mP1);
  float probNoRec = 1.0f - probRec;
  float noErrProb = _impData.noErrProb(mP1);
  float errProb = _impData.errProb(mP1);
  float sum = 0.0;
  int al = _impData.targetAllele(mP1, hap);
  for (int h=0; h < _n; ++h)
  {
    float em = (al == _impData.refAllele(mP1, h)) ? noErrProb : errProb;
    _bwdVal[h] *= em;
    sum += _bwdVal[h];
  }
  float scale = probNoRec/sum;
  float shift = probRec/_n;
  for (int h=0; h < _n; ++h)
    _bwdVal[h] = scale*_bwdVal[h] + shift;
}

void LSHapBaum::setStateProbs(int m, int fwdIndex)
{
  _fwdHapProbs[m].fill(0.0);
  _bwdHapProbs[m].fill(0.0);
  for (int h=0; h < _n; ++h)
  {
    float stateProbs = _fwdVal[fwdIndex][h]*_bwdVal[h];
    _fwdHapProbs[m][_refHapSegs.seq(m+1, h)] += stateProbs;
    _bwdHapProbs[m][_refHapSegs.seq(m, h)] += stateProbs;
  }
  float sumhp = sum(_fwdHapProbs[m]);
  scale(_fwdHapProbs[m], sumhp);
  scale(_bwdHapProbs[m], sumhp);
}

void LSHapBaum::setAlleleProbs()
{
  setFirstAlleleProbs();
  int nSegsM1 = _refHapSegs.nSegs() - 1;
  for (int j=1; j<nSegsM1; ++j)
    setAlleleProbs(j);
  setLastAlleleProbs();
}

void LSHapBaum::setFirstAlleleProbs()
{
  int segment = 0;
  int nSeq = _refHapSegs.nSeq(segment);
  int endRefMarker = _refHapSegs.segStart(segment + 1);
  float threshold = findThreshold(nSeq);
  for (int seq=0; seq<nSeq; ++seq)
  {
    if (_bwdHapProbs[segment][seq] >= threshold)
    {
      for (int m=0; m<endRefMarker; ++m)
      {
        int start = _refMarkers.sumAlleles(m);
        int allele = _refHapSegs.allele(segment, m, seq);
        _alleleProbs[start + allele] += _bwdHapProbs[segment][seq];
      }
    }
  }
}

void LSHapBaum::setAlleleProbs(int segment)
{
  Q_ASSERT_X(segment > 0,
             "LSHapBaum::setAlleleProbs",
             "segment <= 0");

  int clustStart = _refHapSegs.segStart(segment);
  int clustEnd = _refHapSegs.segEnd(segment - 1);
  int nextClustStart = _refHapSegs.segStart(segment + 1);
  int nSeq = _refHapSegs.nSeq(segment);
  float threshold = findThreshold(nSeq);
  for (int seq=0; seq<nSeq; ++seq)
  {
    bool useFwd = _fwdHapProbs[segment-1][seq] >= threshold;
    bool useBwd = _bwdHapProbs[segment][seq] >= threshold;

    if (useFwd)
    {
      for (int m=clustStart; m<clustEnd; ++m)
      {
	int start = _refMarkers.sumAlleles(m);
	int allele = _refHapSegs.allele(segment, m - clustStart, seq);
	_alleleProbs[start + allele] += _fwdHapProbs[segment-1][seq];
      }
    }

    if (useFwd || useBwd)
    {
      for (int m=clustEnd; m<nextClustStart; ++m)
      {
	int start = _refMarkers.sumAlleles(m);
	int allele = _refHapSegs.allele(segment, m - clustStart, seq);
	double wt = _impData.weight(m);
	_alleleProbs[start + allele] += wt*_fwdHapProbs[segment-1][seq];
	_alleleProbs[start + allele] += (1-wt)*_bwdHapProbs[segment][seq];
      }
    }
  }
}

void LSHapBaum::setLastAlleleProbs()
{
  int segment = _refHapSegs.nSegs() - 1;
  int cluster = segment - 1;
  int refMarkerStart = _refHapSegs.segStart(segment);
  int refMarkerEnd = _refHapSegs.segEnd(segment);
  int nSeq = _refHapSegs.nSeq(segment);
  float threshold = findThreshold(nSeq);
  for (int seq=0; seq<nSeq; ++seq)
  {
    if (_fwdHapProbs[cluster][seq] >= threshold)
    {
      for (int m=refMarkerStart; m<refMarkerEnd; ++m)
      {
	int start = _refMarkers.sumAlleles(m);
	int allele = _refHapSegs.allele(segment, m - refMarkerStart, seq);
	_alleleProbs[start + allele] += _fwdHapProbs[cluster][seq];
      }
    }
  }
}

int LSHapBaum::nextIndex()
{
  ++_arrayIndex;
  if (_arrayIndex == _fwdVal.length())
  {
    ++_windowIndex;
    _arrayIndex = _windowIndex;
  }
  return _arrayIndex;
}

int LSHapBaum::previousIndex(int hap)
{
  if (_arrayIndex == _windowIndex)
  {
    --_windowIndex;
    _arrayIndex = _windowIndex;
    int start = _fwdValueIndex2Marker[_arrayIndex] + 1;
    int end = start + ( _fwdVal.length() - (_arrayIndex + 1) );
    setForwardValues(start, end, hap);
    return _arrayIndex;
  }
  else
    return --_arrayIndex;
}

