/* Copyright 2016 Golden Helix, Inc. */
#ifndef SAMPLES_H
#define SAMPLES_H

#include "ghicore/cstring.h"

#include <QSharedData>
#include <QList>
#include <QMap>


namespace SampleNames
{
  int nNames();
  CString name(int index);

  /**
   * Gets the sample name index. Adds the name to the sample names if
   * the name does not already exist.
   */
  int getIndex(CString name);

  /**
   * Gets the sample name index if it exists. Returns -1 otherwise.
   */
  int getIndexIfIndexed(CString name);
};


class SamplesSharedData : public QSharedData
{
public:
  SamplesSharedData() : QSharedData() {}
  SamplesSharedData(const SamplesSharedData &other)
    : QSharedData(other)
  {
    throw("Resetting a SamplesSharedData instance....");
  }
  ~SamplesSharedData() { }

  QList<int> indexToSample;
  QMap<int, int> indexFromSample;
};


class Samples
{
public:
  Samples() { _d = new SamplesSharedData; }
  Samples(const Samples &other) : _d(other._d) { }

  ~Samples() { }

  void setSamp(int sampleIndex);

  int nSamples() const;
  int index(int localIndex) const ;
  int findIndex(int sampleIndex) const;
  CString name(int localIndex) const;

private:
  QSharedDataPointer<SamplesSharedData> _d;
};

Q_DECLARE_TYPEINFO(Samples, Q_MOVABLE_TYPE);

#endif
