#include "impute/singlepermittedstates.h"

#include "impute/samplerdata.h"

void IndexSet::initialize(int max)
{
  Q_ASSERT_X(max >= 0,
             "IndexSet::initialize",
             "max < 0");

  _inSet.fill(false, max+1);
  _indices.fill(0, max+1);
}

bool IndexSet::add(int element)
{
  if (_inSet[element]==false) {
    _indices[_size++] = element;
    _inSet[element]=true;
    return true;
  }
  else {
    return false;
  }
}

void IndexSet::clear()
{
  for (int j=0, n=_size; j<n; ++j) {
    _inSet[_indices[j]] = false;
  }
  _size = 0;
}

int IndexSet::enumeratedValue(int enumIndex) const
{
  Q_ASSERT_X(enumIndex < _size,
             "IndexSet::enumeratedValue",
             "enumIndex>=_size");

  return _indices[enumIndex];
}

SinglePermittedStates::SinglePermittedStates(const RestrictedDag &rdag, int sample)
  : _rdag(rdag), _marker(-1), _i1(0), _i2(0), _edge1(-1), _edge2(-1), _rev(false)
{
  int nHaps;
  QList<HapSegment> list1;
  QList<HapSegment> list2;
  rdag.singleStates(sample, _nMarkers, nHaps, list1, list2);

  int hap1 = 2*sample;
  int hap2 = 2*sample + 1;

  _indices1.initialize(nHaps);
  _indices2.initialize(nHaps);

  QList<HapSegment> extList1;
  QList<HapSegment> extList2;
  extendSegment(extList1, hap1, list1);
  extendSegment(extList2, hap2, list2);

  _tree1.initialize(_nMarkers, extList1);
  _tree2.initialize(_nMarkers, extList2);
}

void SinglePermittedStates::extendSegment(QList<HapSegment> &extendedSegs,
                                          int hap, const QList<HapSegment> &ibsSegs) const
{
  CenteredIntIntervalTree<HapSegment> tree(_nMarkers, ibsSegs);

  // permit states traversed by hap
  int lastMarker = _nMarkers-1;
  extendedSegs.append( HapSegment(hap, 0, lastMarker) );

  // permit states traversed by IBS haps
  for (int k=0, n=ibsSegs.size(); k<n; ++k) {
    HapSegment targetHS = ibsSegs[k];
    int start = _rdag.modifyStart(targetHS, tree);
    int end  = _rdag.modifyEnd(targetHS, tree);
    extendedSegs.append( HapSegment(targetHS.hap(), start, end) );
  }
}

void SinglePermittedStates::convertToIndices(int marker,
                                             const CenteredIntIntervalTree<HapSegment> &tree,
                                             IndexSet &set) const
{
  set.clear();

  QList<HapSegment> hsegs;
  tree.intersect(marker, hsegs);

  foreach(HapSegment hs, hsegs)
    set.add(_rdag._hapStates[marker][hs.hap()]);
}

void SinglePermittedStates::setMarker(int marker)
{
  _marker = marker;
  _i1 = 0;
  _i2 = 0;
  _edge1 = -1;
  _edge2 = -1;
  _rev = false;
  convertToIndices(marker, _tree1, _indices1);
  convertToIndices(marker, _tree2, _indices2);
}

void SinglePermittedStates::next()
{
  Q_ASSERT_X(hasNext(),
             "SinglePermittedStates::next",
             "hasNext()==false");

  if (_rev) {
    int tmp = _edge1;
    _edge1 = _edge2;
    _edge2 = tmp;
    ++_i2;
    if (_i2==_indices2.size()) {
      ++_i1;
      _i2 = 0;
    }
    _rev = false;
  }
  else {
    _edge1 = _indices1.enumeratedValue(_i1);
    _edge2 = _indices2.enumeratedValue(_i2);
    if (!_indices1.contains(_edge2)
        || !_indices2.contains(_edge1)) {
      _rev = true;
    }
    else {
      ++_i2;
      if (_i2==_indices2.size()) {
        ++_i1;
        _i2 = 0;
      }
    }
  }
}

