#include "markers.h"
#include <QString>

#include <QList>
#include <QMap>
#include <QSet>

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
  _d->chromIndex = chromIndex;
  _d->pos = pos;
  _d->id = id;
  _d->nGenotypes = 0;
}

void Marker::setAllele(CString allele)
{
  _d->alleles.append(allele);

  int l = _d->alleles.length();
  _d->nGenotypes = (l * (1 + l)) / 2;
}

CString Marker::chrom() const
{
  return ChromeIds::chromeId(_d->chromIndex);
}

CString Marker::id() const
{
  if (_d->id.length())
    return _d->id;
  else
    return CString(QString("%1:%2").arg((QString) chrom()).arg(pos()));
}

// Do not be concerned about equality of the marker ID's.
bool Marker::operator==(Marker otherMarker) const
{
  if (_d->chromIndex != otherMarker.chromIndex()) {
    return false;
  }
  if (_d->pos != otherMarker.pos()) {
    return false;
  }
  if (_d->alleles != otherMarker.alleles()) {
    return false;
  }
  return true;
}


Markers::Markers(QList<Marker> individualMarkers)
{
  initSharedDataPointers();
        checkMarkerPosOrder(individualMarkers);

        _d->fwdMarkerArray = individualMarkers;
        setReverseMarkers(individualMarkers);

        cumSumAlleles(_d->fwdSumAlleles, individualMarkers);
        cumSumGenotypes(_d->fwdSumGenotypes, individualMarkers);
        cumSumHaplotypeBits(_d->fwdSumHaplotypeBits, individualMarkers);

        cumSumAlleles(_drev->fwdSumAlleles, _drev->fwdMarkerArray);
        cumSumGenotypes(_drev->fwdSumGenotypes, _drev->fwdMarkerArray);
        cumSumHaplotypeBits(_drev->fwdSumHaplotypeBits, _drev->fwdMarkerArray);
}

void Markers::initSharedDataPointers()
{
  _d = new MarkersPluralSharedData;
  _drev = new MarkersPluralSharedData;
}

 void Markers::checkMarkerPosOrder(QList<Marker> markers) {
      if (markers.length() < 2) {
            return;
        }
        QSet<int> chromIndices;
        chromIndices.insert(markers[0].chromIndex());
        chromIndices.insert(markers[1].chromIndex());
        for (int j=2; j<markers.length(); ++j) {
            int chr0 = markers[j-2].chromIndex();
            int chr1 = markers[j-1].chromIndex();
            int chr2 = markers[j].chromIndex();
            if (chr0 == chr1 && chr1==chr2) {
                int pos0 = markers[j-2].pos();
                int pos1 = markers[j-1].pos();
                int pos2 = markers[j].pos();
                if ( (pos1<pos0 && pos1<pos2) || (pos1>pos0 && pos1>pos2) ) {
		  throw (QString("markers not in chromosomal order: \n%1\n%2\n%3")
                 .arg(markers[j - 2].id().asQString())
                 .arg(markers[j - 1].id().asQString())
                 .arg(markers[j].id().asQString()));
                }
            }
            else if (chr1!=chr2) {
                if (chromIndices.contains(chr2)) {
		  throw (QString("markers on chromosome are not contiguous: %1")
                 .arg(ChromeIds::chromeId(chr2).asQString()));
                }
                chromIndices.insert(chr2);
            }
        }
    }

 void Markers::setReverseMarkers(QList<Marker> fwdList)
 {
   for (int jrev=fwdList.length() - 1; jrev >= 0; --jrev)
     _drev->fwdMarkerArray.append( fwdList[jrev] );
 }

 void Markers::cumSumAlleles(QList<int> &sumAlleles, QList<Marker> markers) {
   int cumSum = 0;
        for (int j=0; j < markers.length(); ++j) {
	  sumAlleles.append(cumSum);
          cumSum += markers[j].nAlleles();
        }
	sumAlleles.append(cumSum);
    }

 void Markers::cumSumGenotypes(QList<int> &sumGenotypes, QList<Marker> markers) {
   int cumSum = 0;
        for (int j=0; j < markers.length(); ++j) {
	  sumGenotypes.append(cumSum);
          cumSum += markers[j].nGenotypes();
        }
	sumGenotypes.append(cumSum);
    }

 void Markers::cumSumHaplotypeBits(QList<int> &sumHaplotypeBits, QList<Marker> markers) {
   int cumSum = 0;
        for (int j=0; j < markers.length(); ++j) {
	  sumHaplotypeBits.append(cumSum);
	  int nAllelesM1 = markers[j].nAlleles() - 1;
	  int nStorageBits = 0;
	  while(nAllelesM1 > 0){
	    nStorageBits++;
	    nAllelesM1 >>= 1;
	  }
          cumSum += nStorageBits;
        }
	sumHaplotypeBits.append(cumSum);
    }

 Markers Markers::reverse() const {
   Markers reversedObject = (*this);
   reversedObject._d.swap(reversedObject._drev);
   return reversedObject;
 }

    Markers Markers::restrict(int start, int end) const {
      if (end > _d->fwdMarkerArray.length()) {
	throw (QString("For Markers::restrict: end (%1) > # markers (%2)").arg(end).arg(nMarkers()));
        }
	  return Markers(_d->fwdMarkerArray.mid(start, end - start));
    }
