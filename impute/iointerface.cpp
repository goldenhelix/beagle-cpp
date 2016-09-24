#include "impute/iointerface.h"

#include "impute/haplotypepair.h"

void RefDataReader::makeNewWindow(int overlap)
{
  _oldVcfRefRecs = _vcfRefRecs.mid(_vcfRefRecs.length() - overlap);
  _vcfRefRecs.clear();
  _vcfRefRecs = _oldVcfRefRecs;
}

void RefDataReader::addNewDataToNewWindow(int desiredWindowSize)
{
  int cci = currentChromIndex();
  while (_vcfRefRecs.length() < desiredWindowSize && hasNextRec() &&
         nextRec().marker().chromIndex() == cci) {
    _vcfRefRecs.append(nextRec());
    advanceRec();
  }

  // add all markers at the same marker position
  BitSetRefGT last = _vcfRefRecs[_vcfRefRecs.length() - 1];
  while (hasNextRec() && samePosition(last, nextRec())) {
    _vcfRefRecs.append(nextRec());
    advanceRec();
  }

  QList<Marker> markerlist;
  for (int i = 0; i < _vcfRefRecs.length(); i++)
    markerlist.append(_vcfRefRecs[i].marker());

  _markers = Markers(markerlist);
}

bool RefDataReader::lastWindowOnChrom() const
{
  return (!hasNextRec() || !(sameChrom(nextRec(), _vcfRefRecs[0])));
}

bool RefDataReader::sameChrom(BitSetRefGT a, BitSetRefGT b) const
{
  return a.marker().chromIndex() == b.marker().chromIndex();
}

bool RefDataReader::samePosition(BitSetRefGT a, BitSetRefGT b) const
{
  return sameChrom(a, b) && a.marker().pos() == b.marker().pos();
}

int RefDataReader::currentChromIndex() const
{
  if (_vcfRefRecs.length()) {
    return _vcfRefRecs[0].marker().chromIndex();
  } else if (hasNextRec()) {
    return nextRec().marker().chromIndex();
  } else
    return -1;
}

void TargDataReader::makeNewWindow(int overlap)
{
  _oldVcfEmissions = _vcfEmissions.mid(_vcfEmissions.length() - overlap);
  _vcfEmissions.clear();
  _vcfEmissions = _oldVcfEmissions;
}

void TargDataReader::addNewDataToNewWindow(int desiredWindowSize)
{
  int cci = currentChromIndex();
  while (_vcfEmissions.length() < desiredWindowSize && hasNextRec() &&
         nextRec().marker().chromIndex() == cci) {
    _vcfEmissions.append(nextRec());
    advanceRec();
  }

  // add all markers at the same marker position
  BitSetGT last = _vcfEmissions[_vcfEmissions.length() - 1];
  while (hasNextRec() && samePosition(last, nextRec())) {
    _vcfEmissions.append(nextRec());
    advanceRec();
  }

  QList<Marker> markerlist;
  for (int i = 0; i < _vcfEmissions.length(); i++)
    markerlist.append(_vcfEmissions[i].marker());

  _markers = Markers(markerlist);
}

bool TargDataReader::lastWindowOnChrom() const
{
  return (!hasNextRec() || !(sameChrom(nextRec(), _vcfEmissions[0])));
}

bool TargDataReader::sameChrom(BitSetGT a, BitSetGT b) const
{
  return a.marker().chromIndex() == b.marker().chromIndex();
}

bool TargDataReader::samePosition(BitSetGT a, BitSetGT b) const
{
  return sameChrom(a, b) && a.marker().pos() == b.marker().pos();
}

int TargDataReader::currentChromIndex() const
{
  if (_vcfEmissions.length()) {
    return _vcfEmissions[0].marker().chromIndex();
  } else if (hasNextRec()) {
    return nextRec().marker().chromIndex();
  } else
    return -1;
}

void VcfWindow::advanceWindow(int overlap, int desiredWindowSize, GenericDataReader &dr)
{
  checkParameters(overlap, desiredWindowSize, dr);

  dr.makeNewWindow(overlap);
  dr.addNewDataToNewWindow(desiredWindowSize);
  _overlap = overlap;
  _cumMarkerCnt += (dr.windowSize() - overlap);
}

void VcfWindow::checkParameters(int overlap, int desiredWindowSize, GenericDataReader &dr)
{
  Q_ASSERT_X(overlap >= 0 && overlap < desiredWindowSize,
             "VcfWindow::checkParameters",
             "overlap < 0  or  overlap >= desiredWindowSize");

  Q_ASSERT_X(overlap <= dr.maxWindowSize(), "VcfWindow::checkParameters",
             "overlap > max. window size");

  Q_ASSERT_X(overlap <= 0 || !dr.lastWindowOnChrom(),
             "VcfWindow::checkParameters",
             "overlap > 0 after last window on chromosome");
}

bool TargetData::canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const
{
  return tr.canAdvanceWindow();
}

void TargetData::advanceWindow(int overlap,
                               int desiredWindowSize,
                               TargDataReader &tr,
                               RefDataReader &rr)
{
  _vcfWindow.advanceWindow(overlap, desiredWindowSize, tr);
  _gl = SplicedGL(tr.samples(), tr.vcfRecs());
  // _markers = tr.markers();
  _window++;
}

void TargetData::setCdData(CurrentData &cd, const Par &par, const SampleHapPairs &overlapHaps,
                           const TargDataReader &tr, const RefDataReader &rr)
{
  Q_ASSERT_X(overlapHaps.nHaps() == 0 || tr.samples() == overlapHaps.samples(),
             "TargetData::setCdData",
             "inconsistent samples");

  cd._window = _window;
  cd._initHaps = overlapHaps;
  cd._prevSpliceStart = _vcfWindow.overlap() / 2;
  cd._nextOverlapStart = nextOverlapStart(par.overlap(), tr);
  cd._nextSpliceStart = (tr.windowSize() + cd._nextOverlapStart) / 2;
  cd._nextTargetOverlapStart = cd._nextOverlapStart;
  cd._nextTargetSpliceStart = cd._nextSpliceStart;

  // cd._families = families;
  // cd._weights = new Weights(families);

  cd._targetGL = SplicedGL(overlapHaps, _gl);

  cd._nRefSamples = 0;
  // cd._refSamples is left alone.
  cd._targetSamples = tr.samples();
  cd._allSamples = tr.samples();
  cd._markers = tr.markers();
  cd._targetMarkers = tr.markers();

  for (int i = 0; i < cd._markers.nMarkers(); i++)  // Make trivial 1-to-1 mapping.
  {
    cd._targetMarkerIndex.append(i);
    cd._markerIndex.append(i);
  }

  // cd._restRefHapPairs is left alone.
  // cd._refSampleHapPairs is left alone.
  // cd._restrictedRefSampleHapPairs is left alone.

  // cd._recombRate = recombRate(targetMarkers, genMap, par.mapscale());
}

/* Returns the index of the first marker in the overlap */
int TargetData::nextOverlapStart(int targetOverlap, const TargDataReader &tr)
{
  Q_ASSERT_X(targetOverlap >= 0, "TargetData::nextOverlapStart", "targetOverlap < 0");

  if (targetOverlap == 0 || tr.lastWindowOnChrom())
    return tr.windowSize();

  Markers markers = tr.markers();
  int nextOverlap = tr.windowSize() - targetOverlap;
  if (nextOverlap < 0)
    nextOverlap = 0;
  while (nextOverlap > 0 &&
         markers.marker(nextOverlap).pos() == markers.marker(nextOverlap - 1).pos()) {
    --nextOverlap;
  }

  return nextOverlap;
}

bool AllData::canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const
{
  return rr.canAdvanceWindow();
}

void AllData::advanceWindow(int overlap,
                            int desiredWindowSize,
                            TargDataReader &tr,
                            RefDataReader &rr)
{
}

void AllData::setCdData(CurrentData &cd, const Par &par, const SampleHapPairs &overlapHaps,
                        const TargDataReader &tr, const RefDataReader &rr)
{}

