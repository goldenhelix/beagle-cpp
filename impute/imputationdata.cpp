#include "impute/imputationdata.h"

class HapSegData
{
public:

  HapSegData(const SampleHapPairs &refHapPairs, int start, int end);

  QVector<int> hap2Seq() const { return _hap2seq; }

  QVector<int> seq2Hap() const;

private:
  void highMafUpdate(const SampleHapPairs &refHapPairs, int marker);
  int indexOfAllele(const QList<int> &list, int allele) const;
  void mapAlleleToNewSeq(QList<int> &list, int allele);
  void updateHap2Seq(int h, int seq);

  QVector<int> _hap2seq;
  QList<int> _seq2Cnt;

  QList<int> _emptyList;
  QList< QList<int> > _seq2AlleleMap;
};

HapSegData::HapSegData(const SampleHapPairs &refHapPairs, int start, int end)
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

void HapSegData::highMafUpdate(const SampleHapPairs &refHapPairs, int marker)
{
  for (int j=0, n = _seq2AlleleMap.size(); j<n; ++j)
    _seq2AlleleMap[j].clear();

  for (int h=0; h < _hap2seq.length(); ++h)
  {
    int seq = _hap2seq[h];
    int allele = refHapPairs.allele(marker, h);
    QList<int> &list = _seq2AlleleMap[seq];

    if (list.isEmpty())
    {
      list.append(allele);
      list.append(seq);
    }
    else
    {
      int index = indexOfAllele(list, allele);

      if (index==list.size())
        mapAlleleToNewSeq(list, allele);

      updateHap2Seq(h, list[index+1]);
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

void HapSegData::mapAlleleToNewSeq(QList<int> &list, int allele)
{
  list.append(allele);
  list.append(_seq2AlleleMap.size());
  _seq2AlleleMap.append(_emptyList);
  _seq2Cnt.append(0);
}

void HapSegData::updateHap2Seq(int h, int seq)
{
  _seq2Cnt[_hap2seq[h]]--;
  _hap2seq[h] = seq;
  _seq2Cnt[_hap2seq[h]]++;
}

RefHapSeg::RefHapSeg(const SampleHapPairs &refHapPairs, int start, int end)
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

RefHapSegs::RefHapSegs(const SampleHapPairs &refHapPairs, const QVector<int> &segStart, const QVector<int> &segEnd)
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
