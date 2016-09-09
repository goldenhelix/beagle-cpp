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

int HapPair::allele(QBitArray &bitset, int marker) {
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

