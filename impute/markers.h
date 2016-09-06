/* Copyright 2016 Golden Helix, Inc. */
#ifndef MARKERS_H
#define MARKERS_H

#include "ghicore/cstring.h"
#include "impute/samples.h"

#include <QList>
#include <QSharedData>
#include <QBitArray>


namespace ChromeIds
{
  int nIds();

  CString chromeId(int index);

  int getIndex(CString name);
  int getIndexIfIndexed(CString name);
};


class MarkerSharedData : public QSharedData
{
public:
  MarkerSharedData() : QSharedData(), chromIndex(-1) {}
  MarkerSharedData(const MarkerSharedData &other)
    : QSharedData(other), chromIndex(-1)
  {
    throw("Resetting a MarkerSharedData instance....");
  }
  ~MarkerSharedData() { }

  int chromIndex;
  int pos;
  CString id;
  QList<CString> alleles;
  int nGenotypes;
};


class Marker
{
public:
  Marker() { _d = new MarkerSharedData; }
  Marker(MarkerSharedData *dpt) {_d = dpt;}               // For initializing from a "derived" class.
  Marker(const Marker &other) : _d(other._d) { }

  ~Marker() { }

  void setIdInfo(int chromIndex, int pos, CString id);
  void setAllele(CString allele);

  CString chrom() const;

  int chromIndex() const { return _d->chromIndex; }
  int pos() const { return _d->pos; }

  /**
   * Returns the marker identifier if there is
   * one, else returns chrome + ":" + pos.
   */
  CString id() const;

  /**
   * Returns the number of alleles for the marker, including the REF
   * allele.
   */
  int nAlleles() const { return _d->alleles.length(); }

  /**
   * Returns the specified allele.  The reference allele has index 0.
   */
  CString allele(int index) const { return _d->alleles[index]; }

  /**
   * Returns the alleles.
   */
  QList<CString> alleles() const { return _d->alleles; }

  /**
   * Returns the number of distinct genotypes, which equals
   * nAlleles()*(1 + nAlleles())/2.
   */
  int nGenotypes() const { return _d->nGenotypes; }

  /**
   * Returns true if the other Marker has the same chromosome,
   * position, and allele lists, and returns false otherwise.
   * Equality does not depend on the value of the ID field.
   */
  bool operator==(Marker otherMarker) const;

private:
  QSharedDataPointer<MarkerSharedData> _d;
};

Q_DECLARE_TYPEINFO(Marker, Q_MOVABLE_TYPE);


class MarkersPluralSharedData : public QSharedData
{
public:
  MarkersPluralSharedData() {}
  MarkersPluralSharedData(const MarkersPluralSharedData &other)
    : QSharedData(other)
  {
    throw("Resetting a MarkersPluralSharedData instance....");
  }
  ~MarkersPluralSharedData() { }

  QList<Marker> fwdMarkerArray;
  QList<int> fwdSumAlleles;
  QList<int> fwdSumGenotypes;
  QList<int> fwdSumHaplotypeBits;
};


class Markers
{
public:
  Markers() { initSharedDataPointers(); }
  Markers(const Markers &other) : _d(other._d), _drev(other._drev) { }

  Markers(QList<Marker> individualMarkers);

  ~Markers() { }

    /**
     * Constructs and returns a new {@code Markers} instance that is
     * equivalent to reversing the order of markers in {@code this}.
     */
  Markers reverse() const;

    /**
     * Returns the number of markers.
     */
    int nMarkers() const {
      return _d->fwdMarkerArray.length();
    }

    /**
     * Returns the specified marker.
     */
    Marker marker(int marker) const {
        return _d->fwdMarkerArray[marker];
    }

    /**
     * Returns the list of markers.
     */
    QList<Marker> markers() const {
        return _d->fwdMarkerArray;
    }

    /**
     * Returns a {@code Markers} instance that represents
     * the specified range of marker indices.
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     */
    Markers restrict(int start, int end) const;

    /**
     * Returns the sum of the number of alleles for
     * the markers with index less than the specified index.
     * @param marker a marker index
     */
    int sumAlleles(int marker) const {
        return _d->fwdSumAlleles[marker];
    }

    /**
     * Returns {@code this.sumAlleles(this.nMarkers())}.
     */
    int sumAlleles() const {
      return _d->fwdSumAlleles[_d->fwdMarkerArray.length()];
    }

    /**
     * Returns the sum of the number of possible genotypes for the markers
     * with index less than the specified index.
     * @param marker a marker index
     */
    int sumGenotypes(int marker) const {
        return _d->fwdSumGenotypes[marker];
    }

    /**
     * Returns {@code this.sumGenotypes(this.nMarkers())}.
     */
    int sumGenotypes() const {
        return _d->fwdSumGenotypes[_d->fwdMarkerArray.length()];
    }

    /**
     * Returns the number of bits requires to store a haplotype for the
     * markers with index less than the specified index.
     * @param marker a marker index
     */
    int sumHaplotypeBits(int marker) const {
        return _d->fwdSumHaplotypeBits[marker];
    }

    /**
     * Returns {@code this.sumHaplotypeBits(this.nMarkers())}.
     */
    int sumHaplotypeBits() const {
      return _d->fwdSumHaplotypeBits[_d->fwdMarkerArray.length()];
    }

private:
  void initSharedDataPointers();
  void checkMarkerPosOrder(QList<Marker> markers);
  void setReverseMarkers(QList<Marker> fwdList);
  void cumSumAlleles(QList<int> &sumAlleles, QList<Marker> markers);
  void cumSumGenotypes(QList<int> &sumGenotypes, QList<Marker> markers);
  void cumSumHaplotypeBits(QList<int> &sumHaplotypeBits, QList<Marker> markers);

  QSharedDataPointer<MarkersPluralSharedData> _d;
  QSharedDataPointer<MarkersPluralSharedData> _drev;
};

Q_DECLARE_TYPEINFO(Markers, Q_MOVABLE_TYPE);


#endif
