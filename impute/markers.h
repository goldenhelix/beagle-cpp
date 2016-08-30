/* Copyright 2016 Golden Helix, Inc. */
#ifndef MARKERS_H
#define MARKERS_H

#include "ghicore/cstring.h"

#include <QList>
#include <QSharedData>

class ChromeIds
{
public:
  int nIds();

  static CString chromeId(int index);

  int getIndex(CString name);
  int getIndexIfIndexed(CString name);
};


class MarkerSharedData : public QSharedData
{
public:
  MarkerSharedData() : chromIndex(-1) {}
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
  Marker() { d = new MarkerSharedData; }

  // Let's see if the copy constructor can be automatically compiled.
  // Let's see if the assignment operator can be automatically compiled.

  void setIdInfo(int chromIndex, int pos, CString id);
  void setAllele(CString allele);

  CString chrom() const;

  int chromIndex() const { return d->chromIndex; }
  int pos() const { return d->pos; }

  /**
   * Returns the marker identifier if there is
   * one, else returns chrome + ":" + pos.
   */
  CString id() const;

  /**
   * Returns the number of alleles for the marker, including the REF
   * allele.
   */
  int nAlleles() const { return d->alleles.length(); }

  /**
   * Returns the specified allele.  The reference allele has index 0.
   */
  CString allele(int index) const { return d->alleles[index]; }

  /**
   * Returns the alleles.
   */
  QList<CString> alleles() const { return d->alleles; }

  /**
   * Returns the number of distinct genotypes, which equals
   * nAlleles()*(1 + nAlleles())/2.
   */
  int nGenotypes() const { return d->nGenotypes; }

  /**
   * Returns true if the other Marker has the same chromosome,
   * position, and allele lists, and returns false otherwise.
   * Equality does not depend on the value of the ID field.
   */
  bool operator==(Marker otherMarker) const;

private:
  QSharedDataPointer<MarkerSharedData> d;
};

Q_DECLARE_TYPEINFO(Marker, Q_MOVABLE_TYPE);


/*
class Markers
{
public:
  int nMarkers();
  void setSamp(int sampleIndex);
  int index(int localIndex) { return _indexToSample[localIndex]; }
  int findIndex(int sampleIndex);
  CString name(int localIndex);

private:
  SampleNames _sNamesObject;
  QList<int> _indexToSample;
  QMap<int, int> _indexFromSample;
};
*/

#endif
