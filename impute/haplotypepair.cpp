#include "haplotypepair.h"

    /**
     * Constructs a new {@code HapPair} instance.
     * @param markers the sequence of markers
     * @param samples the list of samples
     * @param sampleIndex the sample index
     * @param alleles1 the sequence of allele indices for the first haplotype
     * @param alleles2 the sequence of alleles indices for the second haplotype
     */
HapPair::HapPair(Markers markers, Samples samples, int sampleIndex,
            QList<int> &alleles1, QList<int> &alleles2)
{
  Q_ASSERT_X(alleles1.length() == markers.nMarkers()
             &&  alleles2.length() == markers.nMarkers(),
	     "HapPair constructor", "inconsistent markers");
  Q_ASSERT_X(sampleIndex >= 0  &&  sampleIndex < samples.nSamples(),
	     "HapPair constructor", "sample index out of bounds");

        _markers = markers;
        _samples = samples;
        _sampleIndex = sampleIndex;
        _alleles1Data = toBitArray(markers, alleles1);
        _alleles2Data = toBitArray(markers, alleles2);
}

HapPair::HapPair(HapPair &other, bool reverse)
{
  _markers = Markers(other.markers(), reverse);
  _samples = other.samples();
  _sampleIndex = other.sampleIndex();

  if(reverse)
  {
    int nmrks = _markers.nMarkers();
    QList<int> newAlList1;
    QList<int> newAlList2;
    for (int mrk = nmrks - 1; mrk >= 0; mrk--)
    {
      newAlList1.append(other.allele1(mrk));
      newAlList2.append(other.allele2(mrk));
    }
    _alleles1Data = toBitArray(_markers, newAlList1);
    _alleles2Data = toBitArray(_markers, newAlList2);
  }
  else
  {
    _alleles1Data = other._alleles1Data;    // Can reach inside the other object's private
    _alleles2Data = other._alleles2Data;    // data because this object is of the same class.
  }
}

QBitArray HapPair::toBitArray(Markers markers, QList<int> &alleles)
{
        int index = 0;
        QBitArray bs(_markers.sumHaplotypeBits(), false);
        for (int k=0; k<alleles.length(); ++k) {
            int allele = alleles[k];
            Q_ASSERT_X(allele >= 0  &&  allele < _markers.marker(k).nAlleles(),
		       "HapPair::toBitArray", "allele out of bounds for a marker");
            int mask = 1;
            int nBits = markers.sumHaplotypeBits(k+1) - markers.sumHaplotypeBits(k);
            for (int l=0; l<nBits; ++l) {
	      bs.setBit(index++, (allele & mask)==mask );
	      mask <<= 1;
            }
        }
        return bs;
}

int HapPair::allele(const QBitArray &bitset, int marker) const
{
        int start = _markers.sumHaplotypeBits(marker);
        int end = _markers.sumHaplotypeBits(marker+1);
        if (end==(start+1)) {
            return bitset.testBit(start) ? 1 : 0;
        }
        int allele = 0;
        int mask = 1;
        for (int j=start; j<end; ++j) {
            if (bitset.testBit(j)) {
                allele += mask;
            }
            mask <<= 1;
        }
        return allele;
}


HapPairs::HapPairs(QList<HapPair> hapPairList, bool reverse)
{
  Q_ASSERT_X(!hapPairList.isEmpty(), "HapPairs::HapPairs", "Empty HapPair list");

  _hapPairs = hapPairList;
  checkAndExtractMarkers(reverse);
}

    /**
     * Checks that all haplotype pairs have alleles for lists of
     * markers of the same lengths (and thus presumably for one and
     * the same marker list), and sets the list of markers into the
     * current object's data.
     *
     * @param reverse whether to have the markers in reverse order
     */
void HapPairs::checkAndExtractMarkers(bool reverse)
{
  _markers = _hapPairs[0].markers();

  if(reverse)
    _markers = Markers(_markers, true);

  for (int j=1, n=_hapPairs.length(); j<n; ++j)
    Q_ASSERT_X(_hapPairs[j].markers() == _markers,
	       "HapPairs::checkAndExtractMarkers",
	       "inconsistent markers");

  _isReversed = reverse;
  _numOfMarkersM1 = _markers.nMarkers() - 1;
}

int HapPairs::allele(int marker, int haplotype) const
{
  int hapPair = haplotype / 2;
  if ((haplotype & 1) == 0) {
    return _hapPairs[hapPair].allele1(marker);
  } else {
    return _hapPairs[hapPair].allele2(marker);
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
 * "Assert" statement throws if both HapPair objects are not
 * associated with a sample contained within the Samples object passed
 * to this class' object.
 *
 * @param samples the list of samples used to compare {@code HapsPair}
 * objects
 */

class CompareSamplesUsed
{
public:
  CompareSamplesUsed( Samples samples ) 
    : _samples(samples) {}

  bool operator()( const HapPair& hp1, const HapPair& hp2 )
  {
    int id1 = hp1.samples().idIndex(hp1.sampleIndex());
    int id2 = hp2.samples().idIndex(hp2.sampleIndex());
    int i1 = _samples.findLocalIndex(id1);
    int i2 = _samples.findLocalIndex(id2);
    Q_ASSERT_X(i1 != -1  &&  i2 != -1,
	       "CompareSamplesUsed::operator()",
	       "sample not listed in the Samples object.");
    return i1 < i2;
  }

private:
  Samples _samples;
};

SampleHapPairs::SampleHapPairs(Samples samples, QList<HapPair> hapPairList, bool reverse)
  : HapPairs(hapPairList, reverse)
{
  _samples = samples;
  checkSamples();
  CompareSamplesUsed csu(_samples);
  qStableSort(_hapPairs.begin(), _hapPairs.end(), csu);
}

void SampleHapPairs::checkSamples()
{
  Q_ASSERT_X(_samples.nSamples() == _hapPairs.size(),
	     "SampleHapPairs::checkSamples", "inconsistent numbers of samples");

  for (int j=0, n=_hapPairs.length(); j<n; ++j) {
    if (!(_samples == _hapPairs[j].samples())) {
      HapPair hp = _hapPairs[j];
      int i1 = _samples.idIndex(j);
      int i2 = hp.samples().idIndex(hp.sampleIndex());
      Q_ASSERT_X(i1 == i2, "SampleHapPairs::checkSamples", "inconsistent samples");
    }
  }
}
