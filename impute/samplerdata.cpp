#include "impute/samplerdata.h"

#include "impute/dag.h"
#include "impute/iointerface.h"
#include "impute/haplotypepair.h"

#include <QLinkedList>
#include <QMutableLinkedListIterator>

#include <float.h>
#include <math.h>

#define MIN_CM_DIST        1e-7
#define END_FILTER         1

// Utility for RestrictedDag:

// For a given start field, the end field is sorted in reverse order.
static bool modStartComparator(const HapSegment &hs1, const HapSegment &hs2)
{
  if (hs1.start() != hs2.start())
    return (hs1.start() < hs2.start());
  else if (hs1.end() != hs2.end())
    return (hs1.end() > hs2.end());
  else
    return (hs1.hap() < hs2.hap());
}

RestrictedDag::RestrictedDag(const ImmutableDag &dag, const SampleHapPairs &haps,
                             double ibdLength, double ibdExtend)
  : _dag(dag), _haps(haps), _ibdExtend(ibdExtend)
{
  Q_ASSERT_X(ibdLength > 0.0,
             "RestrictedDag::RestrictedDag",
             "ibdLength <= 0d");

  Q_ASSERT_X(ibdExtend > 0.0,
             "RestrictedDag::RestrictedDag",
             "ibdExtend <= 0d");

  initPos();
  initHapStates();
  _hapSegments.initialize(_haps, _pos, ibdLength);
}

void RestrictedDag::initPos()
{
  _pos = _dag.posArray();
  double scaleFactor = 0.2;
  for (int j=0; j < _pos.length(); ++j)
    _pos[j] *= scaleFactor;
}

void RestrictedDag::initHapStates()
{
  Q_ASSERT_X(_dag.nLevels()==_haps.nMarkers(),
             "RestrictedDag::initHapStates",
             "_dag.nLevels()!=_haps.nMarkers()");

  int nHaps = _haps.nHaps();
  int nLevels = _dag.nLevels();

  QVector<int> zeroes(nHaps, 0);  // "int[][] _hapStates = new int[nLevels][nHaps];"
  for (int j=0; j < nLevels; j++)
    _hapStates.append(zeroes);

  for (int h=0; h<nHaps; ++h)
  {
    int node = 0;
    int symbol = _haps.allele(0, h);
    _hapStates[0][h] = _dag.outEdgeBySymbol(0, node, symbol);

    for (int j=1; j < nLevels; ++j)
    {
      node = _dag.childNode(j-1, _hapStates[j-1][h]);
      symbol = _haps.allele(j, h);
      _hapStates[j][h] = _dag.outEdgeBySymbol(j, node, symbol);

      Q_ASSERT_X(_hapStates[j][h] != -1,
                 "RestrictedDag::initHapStates",
                 "_hapStates[j][h] == -1");
    }
  }
}

void RestrictedDag::singleStates(int sample, int &outNMarkers, int &outNHaps,
                                 QList<HapSegment> &outHapSegs1, QList<HapSegment> &outHapSegs2) const
{
  Q_ASSERT_X(sample >= 0  &&  sample < _haps.nSamples(),
             "RestrictedDag::singleStates",
             "sample < 0 || sample >= _haps.nSamples()");

  outNMarkers = _haps.nMarkers();
  outNHaps = _haps.nHaps();

  int hap1 = 2*sample;
  int hap2 = 2*sample + 1;
  ibsSegs(outHapSegs1, hap1);
  ibsSegs(outHapSegs2, hap2);
}

void RestrictedDag::ibsSegs(QList<HapSegment> &outHapSegs, int hap) const
{
  _hapSegments.filteredFind(outHapSegs, hap);
  containmentFilter(outHapSegs, END_FILTER);
}

/* filter if minimum requirements are not met */
void RestrictedDag::containmentFilter(QList<HapSegment> &outIbsSegments,
                                      int minEndDiff) const
{
  Q_ASSERT_X(minEndDiff >= 0, "RestrictedDag::containmentFilter", "minEndDiff < 0");

  if (!outIbsSegments.isEmpty())
  {
    qStableSort(outIbsSegments.begin(), outIbsSegments.end(), modStartComparator);

    QLinkedList<HapSegment> list;
    QList<HapSegment> filteredSegments;
    for (int k=0, m=outIbsSegments.size(); k<m; ++k) {
      const HapSegment &hs = outIbsSegments[k];
      bool exclude = false;
      QMutableLinkedListIterator<HapSegment> it(list);
      while (it.hasNext() && exclude==false) {
        const HapSegment &cover = it.next();
        int cStart = cover.start();
        int cEnd = cover.end();
        if (cEnd <= hs.start() ) {
          it.remove();
        }
        else {
          if ( (hs.start() - cStart) >= minEndDiff
               && (cEnd - hs.end()) >= minEndDiff) {
            exclude = true;
          }
        }
      }
      if (exclude==false) {
        list.append(hs);
        filteredSegments.append(hs);
      }
    }
    outIbsSegments = filteredSegments;
  }
}

int RestrictedDag::modifyStart(const HapSegment &targetHS, CenteredIntIntervalTree<HapSegment> &tree) const
{
  int targetHsStart = targetHS.start();
  int maxStart = IbsHapSegUtility::lowerBoundIndex(_pos, targetHsStart, _pos[targetHsStart] - _ibdExtend);

  int minEnd = targetHS.end();

  QList<HapSegment> list;
  tree.intersectAll(maxStart, minEnd, list);
  return list.isEmpty() ? maxStart : targetHS.start();
}

int RestrictedDag::modifyEnd(const HapSegment &targetHS, CenteredIntIntervalTree<HapSegment> &tree) const
{
  int maxStart = targetHS.start();

  int targetHsEnd = targetHS.end();
  double targetValue = _pos[targetHsEnd] + _ibdExtend;
  int minEnd = IbsHapSegUtility::lowerBoundIndex(_pos, targetHsEnd, targetValue);

  if(minEnd == _pos.length()  ||  _pos[minEnd] > targetValue)  // end is inclusive
    --minEnd;

  QList<HapSegment> list;
  tree.intersectAll(maxStart, minEnd, list);
  return list.isEmpty() ? minEnd : targetHS.end();
}

SamplerData::SamplerData(const RestrictedDag &rdag, const Par &par, const CurrentData &cd,
                         bool revMarkers /* , RunStats runStats */ )
  : _rdag(rdag), _par(par), _revMarkers(revMarkers)
{
  _gl = FuzzyGL(cd.targetGL(), par.err(), revMarkers);
  findDagRecombRate(rdag.dag(), par.mapscale());
}

void SamplerData::findDagRecombRate(const ImmutableDag &dag, float xdist)
{
  QList<double> bglDist = dag.posArray();
  for (int j=0; j < bglDist.length(); ++j)
    bglDist[j] *= 0.2;

  double c = -2.0*xdist;

  _recombRate.fill(0.0f, dag.nLevels());

  double lastGenPos = bglDist[0];
  for (int j=1; j<_recombRate.length(); ++j)
  {
    double genPos = bglDist[j];
    double genDist = fmax(abs(genPos - lastGenPos), MIN_CM_DIST);
    _recombRate[j] = -expm1f(c*genDist);
    lastGenPos = genPos;
  }
}

float SamplerData::err() const
{
  return _par.err();
}
