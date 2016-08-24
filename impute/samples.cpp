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

/**
 * Gets the sample name index. Adds the name to the sample names if
 * the name does not already exist.
 */
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

/**
 * Gets the sample name index if it exists. Returns -1 otherwise.
 */
int SampleNames::getIndexIfIndexed(CString name)
{
  if (_nameMap.contains(name))
    return _nameMap[name];
  else
    return -1;
}

int Samples::nSamples()
{
  return _indexToSample.length();
}

void Samples::setSamp(int sampleIndex)
{
  if (_indexFromSample.contains(sampleIndex))
    throw(QString("Duplicate sample index %1.").arg(sampleIndex));

  _indexFromSample.insert(sampleIndex, _indexToSample.length());
  _indexToSample.append(sampleIndex);
}

int Samples::findIndex(int sampleIndex)
{
  if (_indexFromSample.contains(sampleIndex))
    return _indexFromSample[sampleIndex];
  else
    return -1;
}

CString Samples::name(int localIndex)
{
  return _sNamesObject.name(_indexToSample[localIndex]);
}
