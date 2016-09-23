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

  Q_ASSERT_X(
      overlap <= dr.maxWindowSize(), "VcfWindow::checkParameters", "overlap > max. window size");

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
  // _markers = extractMarkers(markerData);
  _window++;
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
