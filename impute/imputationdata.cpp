#include "impute/imputationdata.h"


class HapSegData
{
public:

  HapSegData(SampleHapPairs* refHapPairs, int start, int end);

  QVector<int> hap2Seq() { return _hap2seq; }

  QVector<int> seq2Hap();

private:
  void highMafUpdate(int marker);
  int indexOfAllele(QList<int> list, int allele);
  void mapAlleleToNewSeq(QList<int> list, int allele);
  void updateHap2Seq(int h, int seq);

  SampleHapPairs* _refHapPairs;
  QVector<int> _hap2seq;
  QList<int> _seq2Cnt;

  QList<int> _emptyList;
  QList< QList<int> > _seq2AlleleMap;
};

HapSegData::HapSegData(SampleHapPairs* refHapPairs, int start, int end)
{
  _refHapPairs = refHapPairs;
  _hap2seq.fill(0, refHapPairs->nHaps());
  _seq2Cnt.append(_hap2seq.length());
  _seq2AlleleMap.append(_emptyList);

  for (int m=start; m<end; ++m)
    highMafUpdate(m);
}

QVector<int> HapSegData::seq2Hap()
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

void HapSegData::highMafUpdate(int marker)
{
  for (int j=0, n = _seq2AlleleMap.size(); j<n; ++j)
    _seq2AlleleMap[j].clear();

  for (int h=0; h < _hap2seq.length(); ++h)
  {
    int seq = _hap2seq[h];
    int allele = _refHapPairs->allele(marker, h);
    QList<int> list = _seq2AlleleMap[seq];

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

int HapSegData::indexOfAllele(QList<int> list, int allele)
{
  int index=0;

  while (index < list.size()  &&  list[index] != allele)
    index += 2;

  return index;
}

void HapSegData::mapAlleleToNewSeq(QList<int> list, int allele)
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

RefHapSeg::RefHapSeg(SampleHapPairs* refHapPairs, int start, int end)
{
  Q_ASSERT_X(start >= 0  &&  start < end  &&  end <= refHapPairs->nMarkers(),
             "RefHapSeg::RefHapSeg",
             "invalid start and/or end values");

  HapSegData hapSegData(refHapPairs, start, end);
  _refHapPairs = refHapPairs;
  _start = start;
  _end = end;
  _hapToSeq = hapSegData.hap2Seq();
  _seqToHap = hapSegData.seq2Hap();
}

int RefHapSeg::allele(int marker, int seq)
{
  int refIndex = _start + marker;
  Q_ASSERT_X(marker >= 0  &&  refIndex < _end,
             "RefHapSeg::allele",
             "marker < 0 || refIndex >= end");

  return _refHapPairs->allele(refIndex, _seqToHap[seq]);
}
