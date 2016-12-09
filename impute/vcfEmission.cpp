#include "impute/vcfemission.h"

void BitSetRefGTSharedData::storeAllele(QBitArray &alleles, int sample, int bitsPerAllele,
                                        int allele)
{
  int index = sample * bitsPerAllele;
  int mask = 1;
  for (int k = 0; k < bitsPerAllele; ++k) {
    if ((allele & mask) == mask) {
      alleles.setBit(index);
    }
    ++index;
    mask <<= 1;
  }
}

int BitSetRefGTSharedData::alleleFromBits(const QBitArray &bits, int sample) const
{
  int start = bitsPerAllele * sample;
  int end = start + bitsPerAllele;
  int allele = 0;
  int mask = 1;
  for (int j = start; j < end; ++j) {
    if (bits.testBit(j)) {
      allele += mask;
    }
    mask <<= 1;
  }
  return allele;
}

int BitSetRefGTSharedData::allele1(int sample) const
{
  return bitsPerAllele == 1 ? allele1Data.testBit(sample) : alleleFromBits(allele1Data, sample);
}

int BitSetRefGTSharedData::allele2(int sample) const
{
  return bitsPerAllele == 1 ? allele2Data.testBit(sample) : alleleFromBits(allele2Data, sample);
}

int BitSetGTSharedData::allele1(int sample) const
{
  return isMissing1Data.testBit(sample) ? -1 : (bitsPerAllele == 1
                                                    ? allele1Data.testBit(sample)
                                                    : alleleFromBits(allele1Data, sample));
}

int BitSetGTSharedData::allele2(int sample) const
{
  return isMissing2Data.testBit(sample) ? -1 : (bitsPerAllele == 1
                                                    ? allele2Data.testBit(sample)
                                                    : alleleFromBits(allele2Data, sample));
}

float BitSetGTSharedData::gl(int sample, int a1, int a2) const
{
  bool isMiss1 = isMissing1Data.testBit(sample);
  bool isMiss2 = isMissing2Data.testBit(sample);
  int obsA1 = allele1(sample);
  int obsA2 = allele2(sample);
  bool isPh = isPhasedData.testBit(sample);

  if (isMiss1 && isMiss2) {
    return 1.0;
  } else if (isMiss1 || isMiss2) {
    bool consistent = (obsA1 < 0 || obsA1 == a1) && (obsA2 < 0 || obsA2 == a2);
    if (!isPh && !consistent) {
      consistent = (obsA1 < 0 || obsA1 == a2) && (obsA2 < 0 || obsA2 == a1);
    }
    return consistent ? 1.0 : 0.0;
  } else {
    if (isPh) {
      return (obsA1 == a1 && obsA2 == a2) ? 1.0 : 0.0;
    } else {
      bool isConsistent = (obsA1 == a1 && obsA2 == a2) || (obsA1 == a2 && obsA2 == a1);
      return isConsistent ? 1.0 : 0.0;
    }
  }
}

void BitSetRefGT::setIdInfo(int chromIndex, int pos, const CString &id)
{
  _d->chromIndex = chromIndex;
  _d->pos = pos;
  _d->id = id;
  _d->nGenotypes = 0;
}

void BitSetRefGT::addAllele(const CString &allele)
{
  _d->allelesInMarker.append(allele);

  int l = _d->allelesInMarker.length();
  _d->nGenotypes = (l * (1 + l)) / 2;
}

void BitSetRefGT::storePhasedAlleles(const QVector<int> &als1, const QVector<int> &als2)
{
  int nAllelesM1 = _d->allelesInMarker.length() - 1;
  int nStorageBits = 0;
  while (nAllelesM1 > 0) {
    nStorageBits++;
    nAllelesM1 >>= 1;
  }
  _d->bitsPerAllele = nStorageBits;

  int nSamp = nSamples();
  _d->allele1Data.fill(false, nSamp * nStorageBits);
  _d->allele2Data.fill(false, nSamp * nStorageBits);

  for (int samp = 0; samp < nSamp; ++samp) {
    int a1 = als1[samp];
    int a2 = als2[samp];

    Q_ASSERT_X(a1 != -1 && a2 != -1, "BitSetRefGT::storePhasedAlleles", "Missing reference genotype");

    _d->storeAllele(_d->allele1Data, samp, nStorageBits, a1);
    _d->storeAllele(_d->allele2Data, samp, nStorageBits, a2);
  }
}

float BitSetRefGT::gl(int sample, int al1, int al2) const
{
  bool matches = (al1 == allele1(sample) && al2 == allele2(sample));
  return matches ? 1.0 : 0.0;
}

int BitSetRefGT::allele(int hap) const
{
  int sample = hap / 2;
  return (hap & 1) == 0 ? allele1(sample) : allele2(sample);
}

void BitSetGT::setIdInfo(int chromIndex, int pos, const CString &id)
{
  _d->chromIndex = chromIndex;
  _d->pos = pos;
  _d->id = id;
  _d->nGenotypes = 0;
}

void BitSetGT::addAllele(const CString &allele)
{
  _d->allelesInMarker.append(allele);

  int l = _d->allelesInMarker.length();
  _d->nGenotypes = (l * (1 + l)) / 2;
}

void BitSetGT::storeAlleles(const QVector<int> &als1, const QVector<int> &als2, const QVector<bool> &arePhased)
{
  int nAllelesM1 = _d->allelesInMarker.length() - 1;
  int nStorageBits = 0;
  while (nAllelesM1 > 0) {
    nStorageBits++;
    nAllelesM1 >>= 1;
  }
  _d->bitsPerAllele = nStorageBits;

  bool referenceDataSoFar = true;

  int nSamp = nSamples();
  _d->allele1Data.fill(false, nSamp * nStorageBits);
  _d->allele2Data.fill(false, nSamp * nStorageBits);
  _d->isMissing1Data.fill(false, nSamp);
  _d->isMissing2Data.fill(false, nSamp);
  _d->isPhasedData.fill(false, nSamp);

  for (int samp = 0; samp < nSamp; ++samp) {
    int a1 = als1[samp];
    int a2 = als2[samp];

    if (arePhased[samp])
      _d->isPhasedData.setBit(samp);
    else
      referenceDataSoFar = false;

    if (a1 == -1) {
      _d->isMissing1Data.setBit(samp);
      referenceDataSoFar = false;
    } else
      _d->storeAllele(_d->allele1Data, samp, nStorageBits, a1);

    if (a2 == -1) {
      _d->isMissing2Data.setBit(samp);
      referenceDataSoFar = false;
    } else
      _d->storeAllele(_d->allele2Data, samp, nStorageBits, a2);
  }

  _d->isRefData = referenceDataSoFar;
}

int BitSetGT::allele(int hap) const
{
  int sample = hap / 2;
  return (hap & 1) == 0 ? allele1(sample) : allele2(sample);
}
