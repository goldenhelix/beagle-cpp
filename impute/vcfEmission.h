/* Copyright 2016 Golden Helix, Inc. */
#ifndef VCFEMISSION_H
#define VCFEMISSION_H

#include "impute/markers.h"
#include "impute/samples.h"

class BitSetRefGTSharedData : public MarkerSharedData
{
public:
  BitSetRefGTSharedData() : MarkerSharedData() {}
  BitSetRefGTSharedData(const Samples &samples) : MarkerSharedData(), samples(samples) {}
  BitSetRefGTSharedData(const BitSetRefGTSharedData &other) : MarkerSharedData(other)
  {
    Q_ASSERT_X(false, "BitSetRefGTSharedData", "Resetting a BitSetRefGTSharedData instance....");
  }

  ~BitSetRefGTSharedData() {}
  void storeAllele(QBitArray &alleles, int sample, int bitsPerAllele, int allele);
  int alleleFromBits(const QBitArray &bits, int sample) const;

  int allele1(int sample) const;
  int allele2(int sample) const;

  Samples samples;

  int bitsPerAllele;
  QBitArray allele1Data;
  QBitArray allele2Data;
};

class BitSetGTSharedData : public BitSetRefGTSharedData
{
public:
  BitSetGTSharedData() : BitSetRefGTSharedData() {}
  BitSetGTSharedData(const Samples &samples) : BitSetRefGTSharedData(samples) {}
  BitSetGTSharedData(const BitSetGTSharedData &other) : BitSetRefGTSharedData(other)
  {
    Q_ASSERT_X(false, "BitSetGTSharedData", "Resetting a BitSetGTSharedData instance....");
  }

  ~BitSetGTSharedData() {}
  int allele1(int sample) const;
  int allele2(int sample) const;
  double gl(int sample, int a1, int a2) const;

  bool isRefData;
  QBitArray isMissing1Data;
  QBitArray isMissing2Data;
  QBitArray isPhasedData;
};

class BitSetRefGT
{
public:
  BitSetRefGT() { _d = new BitSetRefGTSharedData; }
  BitSetRefGT(const Samples &samples) { _d = new BitSetRefGTSharedData(samples); }
  BitSetRefGT(const BitSetRefGT &other) : _d(other._d) {}
  void setIdInfo(int chromIndex, int pos, CString id);
  void setAllele(CString allele);

  void storePhasedAlleles(QList<int> &als1, QList<int> &als2);

  /**
   * Returns the number of samples.
   */
  int nSamples() const { return _d->samples.nSamples(); }
  /**
   * Returns the list of samples.
   */
  Samples samples() const { return _d->samples; }
  /**
  * Returns the allele on the specified haplotype.
  * @param haplotype a haplotype index
  */
  int allele(int hap) const;

  /**
   * Returns the first allele for the specified haplotype pair.
   * @param hapPair a haplotype pair index
   */
  int allele1(int sample) const { return _d->allele1(sample); }
  /**
   * Returns the second allele for the specified haplotype pair.
   * @param hapPair a haplotype pair index
   */
  int allele2(int sample) const { return _d->allele2(sample); }
  /**
   * Returns the number of haplotypes.  The returned value is equal to
   * {@code 2*this.nHapPairs()}.
   */
  int nHaps() const { return 2 * _d->samples.nSamples(); }
  /**
   * Returns the number of haplotype pairs.  The returned value is
   * equal to {@code this.nHaps()/2}.
   */
  int nHapPairs() const { return _d->samples.nSamples(); }
  /**
   * Returns the probability of the observed data if the specified pair
   * of ordered alleles is the true genotype in the specified sample.
   * @param sample the sample index
   * @param allele1 the first allele index
   * @param allele2 the second allele index
   */
  double gl(int sample, int allele1, int allele2) const;

  Marker marker() const { return Marker((MarkerSharedData *)&(*_d)); }
  /**
   * Returns the number of alleles for the marker, including the REF
   * allele.
   */
  int nAlleles() const { return _d->allelesInMarker.length(); }
private:
  QSharedDataPointer<BitSetRefGTSharedData> _d;
};

Q_DECLARE_TYPEINFO(BitSetRefGT, Q_MOVABLE_TYPE);

class BitSetGT
{
public:
  BitSetGT() { _d = new BitSetGTSharedData; }
  BitSetGT(const Samples &samples) { _d = new BitSetGTSharedData(samples); }
  BitSetGT(const BitSetGT &other) : _d(other._d) {}
  void setIdInfo(int chromIndex, int pos, CString id);
  void setAllele(CString allele);

  void storeAlleles(QList<int> &als1, QList<int> &als2, QList<bool> &arePhased);

  /**
   * Returns the number of samples.
   */
  int nSamples() const { return _d->samples.nSamples(); }
  /**
   * Returns the list of samples.
   */
  Samples samples() const { return _d->samples; }
  /**
  * Returns the allele on the specified haplotype.
  * @param haplotype a haplotype index
  */
  int allele(int hap) const;

  /**
   * Returns {@code true} if the genotype emission probabilities
   * for each sample are determined by a phased called genotype
   * that has no missing alleles, and returns {@code false} otherwise.
   */
  bool isRefData() const { return _d->isRefData; }
  /**
   * Returns {@code true} if the genotype emission probabilities for
   * the specified sample are determined by a phased, nonmissing genotype,
   * and returns {@code false} otherwise.
   * @param sample the sample index
   */
  bool isPhased(int samp) const { return _d->isPhasedData.testBit(samp); }
  /**
   * Returns the first allele for the specified haplotype pair.
   * @param hapPair a haplotype pair index
   */
  int allele1(int sample) const { return _d->allele1(sample); }
  /**
   * Returns the second allele for the specified haplotype pair.
   * @param hapPair a haplotype pair index
   */
  int allele2(int sample) const { return _d->allele2(sample); }
  /**
   * Returns the number of haplotypes.  The returned value is equal to
   * {@code 2*this.nHapPairs()}.
   */
  int nHaps() const { return 2 * _d->samples.nSamples(); }
  /**
   * Returns the number of haplotype pairs.  The returned value is
   * equal to {@code this.nHaps()/2}.
   */
  int nHapPairs() const { return _d->samples.nSamples(); }
  /**
   * Returns the probability of the observed data if the specified pair
   * of ordered alleles is the true genotype in the specified sample.
   * @param sample the sample index
   * @param allele1 the first allele index
   * @param allele2 the second allele index
   */
  double gl(int samp, int al1, int al2) const { return _d->gl(samp, al1, al2); }
  /**
   * Returns the Marker object associated with this record.
   */
  Marker marker() const { return Marker((MarkerSharedData *)&(*_d)); }
  /**
   * Returns the number of alleles for the marker, including the REF
   * allele.
   */
  int nAlleles() const { return _d->allelesInMarker.length(); }
private:
  QSharedDataPointer<BitSetGTSharedData> _d;
};

Q_DECLARE_TYPEINFO(BitSetGT, Q_MOVABLE_TYPE);

#endif
