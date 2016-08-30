#include "markers.h"
#include <QString>

#include <QList>
#include <QMap>

static QList<CString> _chrIds;
static QMap<CString, int> _chrIdMap;

int ChromeIds::nIds()
{
  return _chrIds.length();
}

CString ChromeIds::chromeId(int index)
{
  return _chrIds[index];
}

/**
 * Gets the chromosome ID index. Adds the name to the chromosome IDs if
 * the name does not already exist.
 */
int ChromeIds::getIndex(CString name)
{
  if (_chrIdMap.contains(name))
    return _chrIdMap[name];
  else {
    int newIndex = _chrIds.length();
    _chrIdMap.insert(name, newIndex);
    _chrIds.append(name);
    return newIndex;
  }
}

/**
 * Gets the chromosome ID index if it exists. Returns -1 otherwise.
 */
int ChromeIds::getIndexIfIndexed(CString name)
{
  if (_chrIdMap.contains(name))
    return _chrIdMap[name];
  else
    return -1;
}


void Marker::setIdInfo(int chromIndex, int pos, CString id)
{
  d->chromIndex = chromIndex;
  d->pos = pos;
  d->id = id;
  d->nGenotypes = 0;
}

void Marker::setAllele(CString allele)
{
  d->alleles.append(allele);

  int l = d->alleles.length();
  d->nGenotypes = (l * (1 + l)) / 2;
}

CString Marker::chrom() const
{
  return ChromeIds::chromeId(d->chromIndex);
}

CString Marker::id() const
{
  if (d->id.length())
    return d->id;
  else
    return CString(QString("%1:%2").arg((QString) chrom()).arg(pos()));
}

// Do not be concerned about equality of the marker ID's.
bool Marker::operator==(Marker otherMarker) const
{
  if (d->chromIndex != otherMarker.chromIndex()) {
    return false;
  }
  if (d->pos != otherMarker.pos()) {
    return false;
  }
  if (d->alleles != otherMarker.alleles()) {
    return false;
  }
  return true;
}

/*
int Markers::nMarkers()
{
  return _indexToSample.length();
}

void Markers::setSamp(int sampleIndex)
{
  if (_indexFromSample.contains(sampleIndex))
    throw(QString("Duplicate sample index %1.").arg(sampleIndex));

  _indexFromSample.insert(sampleIndex, _indexToSample.length());
  _indexToSample.append(sampleIndex);
}

int Markers::findIndex(int sampleIndex)
{
  if (_indexFromSample.contains(sampleIndex))
    return _indexFromSample[sampleIndex];
  else
    return -1;
}

CString Markers::name(int localIndex)
{
  return _sNamesObject.name(_indexToSample[localIndex]);
}
*/
