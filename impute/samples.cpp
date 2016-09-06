#include "samples.h"
#include <QString>


static QList<CString> _names;
static QMap<CString, int> _nameMap;

int SampleNames::nNames()
{
  return _names.length();
}

CString SampleNames::name(int index)
{
  return _names[index];
}

int SampleNames::getIndex(CString name)
{
  if (_nameMap.contains(name))
    return _nameMap[name];
  else {
    int newIndex = _names.length();
    _nameMap.insert(name, newIndex);
    _names.append(name);
    return newIndex;
  }
}

int SampleNames::getIndexIfIndexed(CString name)
{
  if (_nameMap.contains(name))
    return _nameMap[name];
  else
    return -1;
}


void Samples::setSamp(int sampleIndex)
{
  if (_d->indexFromSample.contains(sampleIndex))
    throw(QString("Duplicate sample index %1.").arg(sampleIndex));

  _d->indexFromSample.insert(sampleIndex, _d->indexToSample.length());
  _d->indexToSample.append(sampleIndex);
}

int Samples::nSamples() const
{
  return _d->indexToSample.length();
}

int Samples::index(int localIndex) const
{
  return _d->indexToSample[localIndex];
}

int Samples::findIndex(int sampleIndex) const
{
  if (_d->indexFromSample.contains(sampleIndex))
    return _d->indexFromSample[sampleIndex];
  else
    return -1;
}

CString Samples::name(int localIndex) const
{
  return SampleNames::name(index(localIndex));
}
