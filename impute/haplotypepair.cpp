#include "haplotypepair.h"

HapPair::HapPair(const Markers &markers, const Samples &samples, int sampleIndex,
                 const QList<int> &alleles1, const QList<int> &alleles2)
{
  Q_ASSERT_X(alleles1.length() == markers.nMarkers() && alleles2.length() == markers.nMarkers(),
             "HapPair constructor", "inconsistent markers");
  Q_ASSERT_X(sampleIndex >= 0 && sampleIndex < samples.nSamples(), "HapPair constructor",
             "sample index out of bounds");

  _markers = markers;
  _samples = samples;
  _sampleIndex = sampleIndex;
  _alleles1Data = toBitArray(markers, alleles1);
  _alleles2Data = toBitArray(markers, alleles2);
}

HapPair::HapPair(const HapPair &other, bool reverse)
{
  _markers = Markers(other.markers(), reverse);
  _samples = other.samples();
  _sampleIndex = other.sampleIndex();

  if (reverse) {
    int nmrks = _markers.nMarkers();
    QList<int> newAlList1;
    QList<int> newAlList2;
    for (int mrk = nmrks - 1; mrk >= 0; mrk--) {
      newAlList1.append(other.allele1(mrk));
      newAlList2.append(other.allele2(mrk));
    }
    _alleles1Data = toBitArray(_markers, newAlList1);
    _alleles2Data = toBitArray(_markers, newAlList2);
  } else {
    _alleles1Data = other._alleles1Data;  // Can reach inside the other object's private
    _alleles2Data = other._alleles2Data;  // data because this object is of the same class.
  }
}

QBitArray HapPair::toBitArray(const Markers &markers, const QList<int> &alleles)
{
  int index = 0;
  QBitArray bs(_markers.sumHaplotypeBits(), false);
  for (int k = 0; k < alleles.length(); ++k) {
    int allele = alleles[k];
    Q_ASSERT_X(allele >= 0 && allele < _markers.marker(k).nAlleles(), "HapPair::toBitArray",
               "allele out of bounds for a marker");
    int mask = 1;
    int nBits = markers.sumHaplotypeBits(k + 1) - markers.sumHaplotypeBits(k);
    for (int l = 0; l < nBits; ++l) {
      bs.setBit(index++, (allele & mask) == mask);
      mask <<= 1;
    }
  }
  return bs;
}

int HapPair::allele(const QBitArray &bitset, int marker) const
{
  int start = _markers.sumHaplotypeBits(marker);
  int end = _markers.sumHaplotypeBits(marker + 1);
  if (end == (start + 1)) {
    return bitset.testBit(start) ? 1 : 0;
  }
  int allele = 0;
  int mask = 1;
  for (int j = start; j < end; ++j) {
    if (bitset.testBit(j)) {
      allele += mask;
    }
    mask <<= 1;
  }
  return allele;
}

HapPairs::HapPairs(const QList<HapPair> &hapPairList, bool reverse)
{
  Q_ASSERT_X(!hapPairList.isEmpty(), "HapPairs::HapPairs", "Empty HapPair list");

  _hapPairs = hapPairList;
  checkAndExtractMarkers(reverse);
}

HapPairs::HapPairs(const HapPairs &other)
  : _isReversed(other._isReversed),
    _numOfMarkersM1(other._numOfMarkersM1),
    _markers(other._markers),
    _hapPairs(other._hapPairs)
{
}

/**
 * Checks that all haplotype pairs have alleles for one and
 * the same marker list, and sets the list of markers into the
 * current object's data.
 *
 * @param reverse whether to have the markers in reverse order
 */
void HapPairs::checkAndExtractMarkers(bool reverse)
{
  _markers = _hapPairs[0].markers();

  // Check for self-consistent Markers objects in the HapPair list
  // before reversing the Markers object of this HapPairs object.

  for (int j = 1, n = _hapPairs.length(); j < n; ++j)
    Q_ASSERT_X(_hapPairs[j].markers() == _markers,
               "HapPairs::checkAndExtractMarkers",
               "inconsistent markers");

  if (reverse)
    _markers = Markers(_markers, true);

  _isReversed = reverse;
  _numOfMarkersM1 = _markers.nMarkers() - 1;
}

int HapPairs::allele(int marker, int haplotype) const
{
  int hapPair = haplotype / 2;
  if ((haplotype & 1) == 0) {
    return allele1(marker, hapPair);
  } else {
    return allele2(marker, hapPair);
  }
}

int HapPairs::allele1(int marker, int hapPair) const
{
  if (_isReversed)
    return _hapPairs[hapPair].allele1(_numOfMarkersM1 - marker);
  else
    return _hapPairs[hapPair].allele1(marker);
}

int HapPairs::allele2(int marker, int hapPair) const
{
  if (_isReversed)
    return _hapPairs[hapPair].allele2(_numOfMarkersM1 - marker);
  else
    return _hapPairs[hapPair].allele2(marker);
}

/**
 * A sort compare class whose operator() method returns whether or not
 * the firstHapPair is "less than" the second HapPair. "Less than"
 * means {@code samples.index(hp1.idIndex(hp1.sampleIndex()))} is less
 * than {@code samples.index(hp2.idIndex(hp2.sampleIndex()))}. An
 * "Assert" statement checks if both HapPair objects are
 * associated with a sample contained within the Samples object passed
 * to this class' object.
 *
 * @param samples the list of samples used to compare {@code HapsPair}
 * objects
 */

class CompareSamplesUsed
{
public:
  CompareSamplesUsed(const Samples &samples) : _samples(samples) {}
  bool operator()(const HapPair &hp1, const HapPair &hp2)
  {
    int id1 = hp1.samples().idIndex(hp1.sampleIndex());
    int id2 = hp2.samples().idIndex(hp2.sampleIndex());
    int i1 = _samples.findLocalIndex(id1);
    int i2 = _samples.findLocalIndex(id2);
    Q_ASSERT_X(i1 != -1 && i2 != -1,
               "CompareSamplesUsed::operator()",
               "sample not listed in the Samples object.");
    return i1 < i2;
  }

private:
  Samples _samples;
};

SampleHapPairs::SampleHapPairs(const Samples &samples, const QList<HapPair> &hapPairList,
                               bool reverse)
  : HapPairs(hapPairList, reverse), _samples(samples)
{
  CompareSamplesUsed csu(_samples);  // Sort comparison object.
  qStableSort(_hapPairs.begin(), _hapPairs.end(), csu);

  checkSamples();
}

void SampleHapPairs::checkSamples()
{
  Q_ASSERT_X(_samples.nSamples() == _hapPairs.length(), "SampleHapPairs::checkSamples",
             "inconsistent numbers of samples");

  for (int j = 0, n = _hapPairs.length(); j < n; ++j)
  {
    int i1 = _samples.idIndex(j);
    HapPair &hp = _hapPairs[j];

    if (_samples == hp.samples())
    {
      Q_ASSERT_X(hp.sampleIndex() == j, "SampleHapPairs::checkSamples",
                 "inconsistent samples");
    }
    else
    {
      int i2 = hp.samples().idIndex(hp.sampleIndex());
      Q_ASSERT_X(i1 == i2, "SampleHapPairs::checkSamples", "inconsistent samples");
    }
  }
}

Samples SampleHapPairs::samples(int hapPair) const
{
  Q_ASSERT_X(hapPair >= 0 && hapPair < _samples.nSamples(),
             "SamplesHapPairs::samples(int)",
             "out of bounds HapPair index");
  return _samples;
}

int SampleHapPairs::sampleIndex(int hapPair) const
{
  Q_ASSERT_X(hapPair >= 0 && hapPair < _samples.nSamples(),
             "SampleHapPairs::samples(int)",
             "out of bounds HapPair index");
  return hapPair;
}

RefHapPairs::RefHapPairs(const Samples &samples, const QList<BitSetRefGT> &refVcfRecs)
  : SampleHapPairs(samples), _refVcfRecs(refVcfRecs)
{
  createMarkers();

  for (int j = 0; j < _refVcfRecs.length(); ++j)
    Q_ASSERT_X(_refVcfRecs[j].samples() == _samples, "RefHapPairs::RefHapPairs",
               "sample inconsistency");
}

void RefHapPairs::createMarkers()
{
  _numOfMarkersM1 = _refVcfRecs.length() - 1;

  QList<Marker> markerList;
  for (int j = 0; j <= _numOfMarkersM1; ++j)
    markerList.append(_refVcfRecs[j].marker());
  _markers = Markers(markerList);
}

int RefHapPairs::allele(int marker, int haplotype) const
{
  int hapPair = haplotype / 2;
  if ((haplotype & 1) == 0)
    return allele1(marker, hapPair);
  else
    return allele2(marker, hapPair);
}

GLSampleHapPairs::GLSampleHapPairs(const GLSampleHapPairs &otherGL, bool checkRef, bool reverse)
  : SampleHapPairs(otherGL),
    _overlap(otherGL._overlap),
    _vcfRecs(otherGL._vcfRecs),
    _numOfGlMarkersM1(otherGL._numOfGlMarkersM1)
{
  // Copy constructor for GLSampleHapPairs, which also acts as a
  // "utility" constructor for the SplicedGL copy-and-maybe-reverse
  // constructor.

  if (checkRef)
    Q_ASSERT_X(otherGL.isRefData(),
               "GLSampleHapPairs::GLSampleHapPairs",
               "other GL is not reference data");

  _glIsReversed = (reverse) ? (!otherGL._glIsReversed) : otherGL._glIsReversed;

  _glMarkers = Markers(otherGL.markers(), reverse);
}

GLSampleHapPairs::GLSampleHapPairs(const Samples &samples)
  : SampleHapPairs(samples), _overlap(0), _glIsReversed(false), _numOfGlMarkersM1(-1)
{
  // Partial samples-object-only otherwise-empty constructor.
}

GLSampleHapPairs::GLSampleHapPairs(const SampleHapPairs &haps, const GLSampleHapPairs &otherGL)
  : SampleHapPairs(haps),
    _vcfRecs(otherGL._vcfRecs),
    _numOfGlMarkersM1(otherGL._numOfGlMarkersM1),
    _glIsReversed(otherGL._glIsReversed),
    _glMarkers(otherGL.markers())
{
  // Does the real work of the "SplicedGL constructor".

  // Not necessarily needed, but let's check if everything stays clean, here.
  Q_ASSERT_X(otherGL._numOfMarkersM1 == -1,
             "SplicedGL::SplicedGL(haps, otherGL)",
             "OtherGL already has overlap information in it.");

  // Make sure everything that is supposed to match between the two
  // parts of this class' data matches.

  for (int j = 0, n = SampleHapPairs::nMarkers(); j < n; ++j)
    Q_ASSERT_X(SampleHapPairs::marker(j) == marker(j),
               "SplicedGL::SplicedGL(haps, otherGL)",
               "inconsistent markers");

  if (samples().nSamples()) {
    // Only test if we had info. in the haps object.
    Q_ASSERT_X(samples() == otherGL.samples(),
               "SplicedGL::SplicedGL(haps, otherGL)",
               "inconsistent samples");
  } else
    // Else use otherGL's Samples object.
    _samples = otherGL._samples;

  _overlap = SampleHapPairs::nMarkers();
}

int GLSampleHapPairs::allele1(int marker, int sampNum) const
{
  if (_glIsReversed) {
    if ((_numOfGlMarkersM1 - marker) < _overlap)
      return HapPairs::allele1(_numOfGlMarkersM1 - marker, sampNum);
    else
      return _vcfRecs[_numOfGlMarkersM1 - marker].allele1(sampNum);
  } else {
    if (marker < _overlap)
      return HapPairs::allele1(marker, sampNum);
    else
      return _vcfRecs[marker].allele1(sampNum);
  }
}

int GLSampleHapPairs::allele2(int marker, int sampNum) const
{
  if (_glIsReversed) {
    if ((_numOfGlMarkersM1 - marker) < _overlap)
      return HapPairs::allele2(_numOfGlMarkersM1 - marker, sampNum);
    else
      return _vcfRecs[_numOfGlMarkersM1 - marker].allele2(sampNum);
  } else {
    if (marker < _overlap)
      return HapPairs::allele2(marker, sampNum);
    else
      return _vcfRecs[marker].allele2(sampNum);
  }
}

int GLSampleHapPairs::allele(int marker, int haplotype) const
{
  int sample = haplotype / 2;
  if ((haplotype & 1) == 0) {
    return allele1(marker, sample);
  } else {
    return allele2(marker, sample);
  }
}

SplicedGL::SplicedGL(const Samples &samples, const QList<BitSetGT> &vma) : GLSampleHapPairs(samples)
{
  // "BasicGL constructor"....
  _glIsReversed = false;
  _vcfRecs = vma;
  checkSamples();
  createMarkers();
  checkIfIsRefData();
}

SplicedGL::SplicedGL(const SampleHapPairs &haps, const GLSampleHapPairs &otherGL)
  : GLSampleHapPairs(haps, otherGL)
{
  // "SplicedGL constructor". But also see GLSampleHapPairs(haps,
  // otherGL) above.

  _isRefData = otherGL.isRefData();
}

SplicedGL::SplicedGL(const GLSampleHapPairs &otherGL, bool reverse)
  : GLSampleHapPairs(otherGL, false, reverse), _isRefData(otherGL.isRefData())
{}

void SplicedGL::checkSamples()
{
  for (int j = 0, n = _vcfRecs.length(); j < n; ++j)
    Q_ASSERT_X(_vcfRecs[j].samples() == _samples, "SplicedGL::checkSamples",
               "inconsistent samples");
}

void SplicedGL::createMarkers()
{
  _numOfGlMarkersM1 = _vcfRecs.length() - 1;

  QList<Marker> markerList;
  for (int j = 0; j <= _numOfGlMarkersM1; ++j)
    markerList.append(_vcfRecs[j].marker());
  _glMarkers = Markers(markerList);
}

void SplicedGL::checkIfIsRefData()
{
  _isRefData = true;
  for (int j = 0; j < _vcfRecs.length() && _isRefData; ++j) {
    if (!_vcfRecs[j].isRefData()) {
      _isRefData = false;
      break;
    }
  }
}

float SplicedGL::gl(int marker, int sample, int allele1, int allele2) const
{
  if (_glIsReversed) {
    if ((_numOfGlMarkersM1 - marker) < _overlap) {
      int a1 = HapPairs::allele1(_numOfGlMarkersM1 - marker, sample);
      int a2 = HapPairs::allele2(_numOfGlMarkersM1 - marker, sample);
      return (allele1 == a1 && allele2 == a2) ? 1.0f : 0.0f;
    } else
      return _vcfRecs[_numOfGlMarkersM1 - marker].gl(sample, allele1, allele2);
  } else {
    if (marker < _overlap) {
      int a1 = HapPairs::allele1(marker, sample);
      int a2 = HapPairs::allele2(marker, sample);
      return (allele1 == a1 && allele2 == a2) ? 1.0f : 0.0f;
    } else
      return _vcfRecs[marker].gl(sample, allele1, allele2);
  }
}

bool SplicedGL::isPhased(int marker, int sample) const
{
  if (_glIsReversed) {
    if ((_numOfGlMarkersM1 - marker) < _overlap)
      return true;
    else
      return _vcfRecs[_numOfGlMarkersM1 - marker].isPhased(sample);
  } else {
    if (marker < _overlap)
      return true;
    else
      return _vcfRecs[marker].isPhased(sample);
  }
}

FuzzyGL::FuzzyGL(const SplicedGL &gl, float err, bool reverse) : SplicedGL(gl, reverse)
{
  Q_ASSERT_X(err >= 0.0 && err < 1.0, "FuzzyGL::FuzzyGL", "err out of range.");

  float e = err;
  float f = 1.0 - err;
  _ee = e * e;
  _ef = e * f;
  _ff = f * f;
}

float FuzzyGL::gl(int marker, int sample, int a1, int a2) const
{
  // following algorithm is for both diallelic and multi-allelic markers
  int obs1 = allele1(marker, sample);
  int obs2 = allele2(marker, sample);

  if (obs1 >= 0 && obs2 >= 0) {
    if (obs1 == obs2 || isPhased(marker, sample)) {
      return phasedGL(obs1, obs2, a1, a2);
    } else {
      return phasedGL(obs1, obs2, a1, a2) + phasedGL(obs2, obs1, a1, a2);
    }
  } else
    return SplicedGL::gl(marker, sample, a1, a2);
}

float FuzzyGL::phasedGL(int obs1, int obs2, int a1, int a2) const
{
  if (obs1 == a1) {
    return obs2 == a2 ? _ff : _ef;
  } else {
    return obs2 == a2 ? _ef : _ee;
  }
}

MaskedEndsGL::MaskedEndsGL(const SplicedGL &gl, int start, int end)
    : _gl(gl), _start(start), _end(end)
{
  Q_ASSERT_X(start >= 0  &&  start <= end  &&  end <= gl.nMarkers(),
             "MaskedEndsGL::MaskedEndsGL",
             "start<0 || start>end || end>gl.nMarkers()");
}

void MaskedEndsGL::checkMarkerAndSample(int marker, int sample) const
{
  Q_ASSERT_X(marker >= 0  &&  marker < _gl.nMarkers(),
             "MaskedEndsGL::checkMarkerAndSample",
             "marker<0 || marker>=_gl.nMarkers()");

  Q_ASSERT_X(sample >= 0  &&  sample < _gl.nSamples(),
             "MaskedEndsGL::checkMarkerAndSample",
             "sample<0 || sample>=_gl.nSamples()");
}

void MaskedEndsGL::checkAllele(int marker, int allele) const
{
  Q_ASSERT_X(allele<0 || allele>=_gl.marker(marker).nAlleles(),
             "MaskedEndsGL::checkAllele",
             "allele<0 || allele>=_gl.marker(marker).nAlleles()");
}

float MaskedEndsGL::gl(int marker, int sample, int allele1, int allele2) const
{
  if (marker < _start  ||  marker >= _end)
  {
    checkMarkerAndSample(marker, sample);
    checkAllele(marker, allele1);
    checkAllele(marker, allele2);
    return 1.0f;
  }
  else
    return _gl.gl(marker, sample, allele1, allele2);
}

bool MaskedEndsGL::isPhased(int marker, int sample) const
{
  if (marker < _start  ||  marker >= _end)
  {
    checkMarkerAndSample(marker, sample);
    return false;
  }
  else
    return _gl.isPhased(marker, sample);
}

int MaskedEndsGL::allele1(int marker, int sample) const
{
  if (marker < _start  ||  marker >= _end)
  {
    checkMarkerAndSample(marker, sample);
    return -1;
  }
  else
    return _gl.allele1(marker, sample);
}

int MaskedEndsGL::allele2(int marker, int sample) const
{
  if (marker < _start  ||  marker >= _end)
  {
    checkMarkerAndSample(marker, sample);
    return -1;
  }
  else
    return _gl.allele2(marker, sample);
}

QList<HapPair> HapUtility::createHapPairList(const Markers &markers,
                                             const SampleHapPairs &targetHapPairs,
                                             const QList<int> &mapping)
{
  Samples samples = targetHapPairs.samples();
  int nSamples = samples.nSamples();
  int nMarkers = markers.nMarkers();
  QList<HapPair> list;
  QList<int> a1;
  QList<int> a2;
  for (int m = 0; m < nMarkers; ++m) {
    a1.append(0);
    a2.append(0);
  }
  for (int s = 0; s < nSamples; ++s) {
    for (int m = 0; m < nMarkers; ++m) {
      a1[m] = targetHapPairs.allele1(mapping[m], s);
      a2[m] = targetHapPairs.allele2(mapping[m], s);
    }
    list.append(HapPair(markers, samples, s, a1, a2));
  }
  return list;
}

