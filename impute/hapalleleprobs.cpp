#include "impute/hapalleleprobs.h"

static signed char convertToSignedChar(float f)
{
  if (f >= 1.0) {
    f = (float) 0.99999;
  }
  int bin = ((int) (f*256)) - 128;
  return (signed char) bin;
}

static float INCREMENT = (1.0 / 256);

static float convertToFloat(signed char b)
{
  return (b + 128.5f)*INCREMENT;
}

HapAlleleProbs::HapAlleleProbs(const Markers &markers, const Samples &samples, int hap,
                               const QVector<float> &alleleProbs)
{
  Q_ASSERT_X(alleleProbs.length() == markers.sumAlleles(),
             "HapAlleleProbs::HapAlleleProbs",
             "inconsistent data");
  Q_ASSERT_X(hap >= 0  &&  hap < 2*samples.nSamples(),
             "HapAlleleProbs::HapAlleleProbs",
             "index out of bounds");

  _alleleBin.fill(alleleProbs.length() - markers.nMarkers());
  int index = 0;
  for (int m=0, n=markers.nMarkers(); m<n; ++m)
  {
    float sum = 0.0;
    int start = markers.sumAlleles(m);
    int end = markers.sumAlleles(m+1);
    for (int j=start; j<end; ++j)
    {
      float p = alleleProbs[j];
      Q_ASSERT_X(p >= 0  &&  p <= 1.01,
                 "HapAlleleProbs::HapAlleleProbs",
                 "illegal value in alleleProbs array");

      sum += alleleProbs[j];
    }

    Q_ASSERT_X(sum <= 1.01,  
               "HapAlleleProbs::HapAlleleProbs",
               "sum of allele probs. > 1.01");

    for (int j=start; j<end - 1; ++j)
      _alleleBin[index++] = convertToSignedChar(alleleProbs[j] / sum);
  }

  Q_ASSERT_X(index == _alleleBin.length(),
             "HapAlleleProbs::HapAlleleProbs",
             "index != _alleleBin.length()");

  _markers = markers;
  _samples = samples;
  _hap = hap;
}

float HapAlleleProbs::allele(int marker, int allele) const
{
  int nAlleles = _markers.marker(marker).nAlleles();

  Q_ASSERT_X(allele >= 0  &&  allele < nAlleles,
             "HapAlleleProbs::allele",
             "allele number out of range");

  int start = _markers.sumAlleles(marker) - marker;
  if (nAlleles == 2)
  {
    float f = convertToFloat(_alleleBin[start]);
    return (allele == 0) ? f : (1.0 - f);
  }
  else if (allele == nAlleles - 1)
    return lastAlleleProb(marker);

  else
    return convertToFloat(_alleleBin[start + allele]);
}

float HapAlleleProbs::lastAlleleProb(int marker) const
{
  int nAlleles = _markers.marker(marker).nAlleles();
  int start = _markers.sumAlleles(marker) - marker;
  int end = start + nAlleles - 1;
  float sum = 1.0;
  for (int j = start; j < end; ++j) {
    sum -= convertToFloat(_alleleBin[j]);
  }
  return (sum < 0.0) ? 0.0 : sum;
}

int HapAlleleProbs::alleleWithMaxProb(int marker) const
{
  int nAlleles = _markers.marker(marker).nAlleles();
  int start = _markers.sumAlleles(marker) - marker;

  if (nAlleles == 2)
    return (_alleleBin[start] >= 0) ? 0 : 1;
  else
  {
    int bestIndex = start;
    int end = start + nAlleles - 1;
    float sumProb = 0.0;
    for (int j = start; j<end; ++j) {
      sumProb += convertToFloat(_alleleBin[j]);
      if (_alleleBin[j] > _alleleBin[bestIndex])
	bestIndex = j;
    }

    if ( sumProb < 0.5)
      return nAlleles - 1;
    else
      return bestIndex - start;
  }
}


static bool hapAPComparator(const HapAlleleProbs &p1, const HapAlleleProbs &p2)
{
  return (p1.hapIndex() < p2.hapIndex());
}

// "Normal" constructor:
ConstrainedAlleleProbs::ConstrainedAlleleProbs(SampleHapPairs shp, QList<HapAlleleProbs> alProbs,
                                               QList<int> indexMap)
{
  Q_ASSERT_X(alProbs.length() > 0,
             "ConstrainedAlleleProbs::ConstrainedAlleleProbs",
             "alProbs.length()==0");

  _alleleProbs = alProbs;
  qStableSort(_alleleProbs.begin(), _alleleProbs.end(), hapAPComparator);

  _markers = _alleleProbs[0].markers();
  _samples = _alleleProbs[0].samples();

  for (int j=1; j<_alleleProbs.length(); ++j)
  {
    Q_ASSERT_X(_alleleProbs[j].markers() == _markers,
             "ConstrainedAlleleProbs::ConstrainedAlleleProbs",
             "inconsistent markers between allele probs");

    Q_ASSERT_X(_alleleProbs[j].samples() == _samples,
             "ConstrainedAlleleProbs::ConstrainedAlleleProbs",
             "inconsistent samples between allele probs");
  }

  Q_ASSERT_X(shp.markers() == _markers,
             "ConstrainedAlleleProbs::ConstrainedAlleleProbs",
             "markers inconsistent between allele probs and shp");

  Q_ASSERT_X(shp.samples() == _samples,
             "ConstrainedAlleleProbs::ConstrainedAlleleProbs",
             "samples inconsistent between allele probs and shp");

  Q_ASSERT_X(indexMap.length() == _markers.nMarkers(),
             "ConstrainedAlleleProbs::ConstrainedAlleleProbs",
             "indexMap.length() != _markers.nMarkers()");

  for (int j=0; j<indexMap.length(); ++j)
  {
    if (indexMap[j] != -1)
    {
      Q_ASSERT_X(_markers.marker(j) == shp.marker(indexMap[j]),
             "ConstrainedAlleleProbs::ConstrainedAlleleProbs",
             "markers inconsistent with index map");
    }
  }

  _shp = shp;
  _indexMap = indexMap;
}

// SampleHapPairs-based constructor:
ConstrainedAlleleProbs::ConstrainedAlleleProbs(SampleHapPairs shp)
{
  _shp = shp;
  _markers = shp.markers();
  _samples = shp.samples();

  // Create trivial index map.
  for(int j=0, n=_markers.nMarkers(); j<n; j++)
    _indexMap.append(j);
}

float ConstrainedAlleleProbs::alProb1(int marker, int sample, int allele) const
{
  int targetMarker = _indexMap[marker];
  if (targetMarker == -1)
  {
    Q_ASSERT_X(_alleleProbs[2*sample].hapIndex()/2 == sample,
               "ConstrainedAlleleProbs::alProb1",
               "inconsistent hapIndex");
    return _alleleProbs[2*sample].allele(marker, allele);
  }
  else
    return (_shp.allele1(targetMarker, sample) == allele) ? 1.0 : 0.0;
}

float ConstrainedAlleleProbs::alProb2(int marker, int sample, int allele) const
{
  int targetMarker = _indexMap[marker];
  if (targetMarker == -1)
  {
    Q_ASSERT_X(_alleleProbs[2*sample + 1].hapIndex()/2 == sample,
               "ConstrainedAlleleProbs::alProb2",
               "inconsistent hapIndex");
    return _alleleProbs[2*sample + 1].allele(marker, allele);
  }
  else
    return (_shp.allele2(targetMarker, sample) == allele) ? 1.0 : 0.0;
}

int ConstrainedAlleleProbs::allele1(int marker, int sample) const
{
  int targetMarker = _indexMap[marker];
  if (targetMarker == -1)
  {
    return _alleleProbs[2*sample].alleleWithMaxProb(marker);
  }
  else
    return _shp.allele1(targetMarker, sample);
}

int ConstrainedAlleleProbs::allele2(int marker, int sample) const
{
  int targetMarker = _indexMap[marker];
  if (targetMarker == -1)
  {
    return _alleleProbs[2*sample + 1].alleleWithMaxProb(marker);
  }
  else
    return _shp.allele2(targetMarker, sample);
}

