/* Copyright 2016 Golden Helix, Inc. */
#ifndef SAMPLES_H
#define SAMPLES_H

#include "ghicore/cstring.h"

#include <QList>
#include <QMap>

class SampleNames
{
public:
  int nNames();
  CString name(int index);

  int getIndex(CString name);
  int getIndexIfIndexed(CString name);
};

class Samples
{
public:
  int nSamples();
  void setSamp(int sampleIndex);
  int index(int localIndex) { return _indexToSample[localIndex]; }
  int findIndex(int sampleIndex);
  CString name(int localIndex);

private:
  SampleNames _sNamesObject;
  QList<int> _indexToSample;
  QMap<int, int> _indexFromSample;
};

#endif
