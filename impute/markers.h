/* Copyright 2016 Golden Helix, Inc. */
#ifndef MARKERS_H
#define MARKERS_H

#include "ghicore/cstring.h"

#include <QList>

class ChromeIds
{
public:
  int nIds();

  static CString chromeId(int index);

  int getIndex(CString name);
  int getIndexIfIndexed(CString name);
};

class Marker
{
public:
  Marker() : _chromIndex(-1) {}

  void setIdInfo(int chromIndex, int pos, CString id);
  void setAllele(CString allele);

  CString chrom();

  int chromIndex() { return _chromIndex; }
  int pos() { return _pos; }
  /**
   * Returns the marker identifier if there is
   * one, else returns chrome + ":" + pos.
   */
  CString id();

  /**
   * Returns the number of alleles for the marker, including the REF
   * allele.
   */
  int nAlleles() { return _alleles.length(); }
  /**
   * Returns the specified allele.  The reference allele has index 0.
   */
  CString allele(int index) { return _alleles[index]; }
  /**
   * Returns the alleles.
   */
  QList<CString> alleles() { return _alleles; }
  /**
   * Returns the number of distinct genotypes, which equals
   * nAlleles()*(1 + nAlleles())/2.
   */
  int nGenotypes();

  /**
   * Returns true if the other Marker has the same chromosome,
   * position, and allele lists, and returns false otherwise.
   * Equality does not depend on the value of the ID field.
   */
  bool operator==(Marker otherMarker);

private:
  int _chromIndex;
  int _pos;
  CString _id;
  QList<CString> _alleles;
  int _nGenotypes;
};

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
