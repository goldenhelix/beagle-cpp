#include "impute/iointerface.h"

#include "impute/hapalleleprobs.h"

#define MIN_R2_DEN      1e-8

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

bool RefDataReader::sameChrom(const BitSetRefGT &a, const BitSetRefGT &b) const
{
  return a.marker().chromIndex() == b.marker().chromIndex();
}

bool RefDataReader::samePosition(const BitSetRefGT &a, const BitSetRefGT &b) const
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

bool TargDataReader::sameChrom(const BitSetGT &a, const BitSetGT &b) const
{
  return a.marker().chromIndex() == b.marker().chromIndex();
}

bool TargDataReader::samePosition(const BitSetGT &a, const BitSetGT &b) const
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

QList<int> TargDataReader::restrictedMakeNewWindow(const QList<int> &oldRefIndices,
                                                   int overlapStart)
{
  int oldLen = oldRefIndices.length();
  int restrictedOverlapStart = oldLen - 1;
  while (restrictedOverlapStart >= 0 && oldRefIndices[restrictedOverlapStart] >= overlapStart)
    restrictedOverlapStart--;
  restrictedOverlapStart++;

  QList<int> overlapRefIndices;
  for (int oldIndex = restrictedOverlapStart; oldIndex < oldLen; oldIndex++)
    overlapRefIndices.append(oldRefIndices[oldIndex] - overlapStart);

  _oldVcfEmissions = _vcfEmissions.mid(restrictedOverlapStart);
  _vcfEmissions.clear();
  _vcfEmissions = _oldVcfEmissions;

  return overlapRefIndices;
}

void TargDataReader::restrictedAdvanceWindow(QList<int> &refIndices, int refOverlap,
                                             const Markers &nextMarkers)
{
  for (int j = refOverlap, n = nextMarkers.nMarkers(); j < n; ++j) {
    Marker m = nextMarkers.marker(j);
    if (hasNextRec() && nextRec().marker().chromIndex() == m.chromIndex()) {
      while (hasNextRec() && nextRec().marker().pos() < m.pos())
        advanceRec();

      while (hasNextRec() && nextRec().marker().pos() == m.pos() &&
             (nextRec().marker() == m) == false)
        advanceRec();
    }
    if (hasNextRec() && nextRec().marker() == m) {
      _restrictedCumMarkerCnt++;
      _vcfEmissions.append(nextRec());
      refIndices.append(j);
      advanceRec();
    }
  }

  Q_ASSERT_X(_vcfEmissions.length(),
             "TargDataReader::restrictedAdvanceWindow",
             "No target markers found that overlap the reference markers.");

  QList<Marker> markerlist;
  for (int i = 0; i < _vcfEmissions.length(); i++)
    markerlist.append(_vcfEmissions[i].marker());

  _markers = Markers(markerlist);
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

  Q_ASSERT_X(overlap <= 0 || !dr.lastWindowOnChrom(),
             "VcfWindow::checkParameters",
             "overlap > 0 after last window on chromosome");
}

int InputData::nextOverlapStart(int targetOverlap, const GenericDataReader &gr)
{
  Q_ASSERT_X(targetOverlap >= 0, "InputData::nextOverlapStart", "targetOverlap < 0");

  if (targetOverlap == 0 || gr.lastWindowOnChrom())
    return gr.windowSize();

  Markers markers = gr.markers();
  int nextOverlap = gr.windowSize() - targetOverlap;
  if (nextOverlap < 0)
    nextOverlap = 0;
  while (nextOverlap > 0 &&
         markers.marker(nextOverlap).pos() == markers.marker(nextOverlap - 1).pos()) {
    --nextOverlap;
  }

  return nextOverlap;
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
    cd._markerIndex.append(i);

  // cd._restRefHapPairs is left alone.
  // cd._refSampleHapPairs is left alone.
  // cd._restrictedRefSampleHapPairs is left alone.

  // cd._recombRate = recombRate(targetMarkers, genMap, par.mapscale());
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
  int oldWindowSize = rr.windowSize();

  _vcfWindow.advanceWindow(overlap, desiredWindowSize, rr);
  _refSampleHapPairs = RefHapPairs(rr.samples(), rr.refRecs());

  _refIndices = tr.restrictedMakeNewWindow(_refIndices, oldWindowSize - overlap);
  tr.restrictedAdvanceWindow(_refIndices, overlap, rr.markers());

  QList<int> oneToOneMapping;
  for (int i = 0; i < rr.windowSize(); i++)
    oneToOneMapping.append(i);
  _refHapPairs = HapUtility::createHapPairList(rr.markers(), _refSampleHapPairs, oneToOneMapping);
  _targetRefHapPairs =
      HapUtility::createHapPairList(tr.markers(), _refSampleHapPairs, _refIndices);

  _gl = SplicedGL(tr.samples(), tr.vcfRecs());
  _window++;
}

void AllData::setCdData(CurrentData &cd, const Par &par, const SampleHapPairs &overlapHaps,
                        const TargDataReader &tr, const RefDataReader &rr)
{
  Q_ASSERT_X(overlapHaps.nHaps() == 0 || tr.samples() == overlapHaps.samples(),
             "AllData::setCdData",
             "inconsistent samples");

  cd._window = _window;
  cd._initHaps = overlapHaps;
  cd._prevSpliceStart = _vcfWindow.overlap() / 2;
  cd._nextOverlapStart = nextOverlapStart(par.overlap(), rr);
  cd._nextSpliceStart = (rr.windowSize() + cd._nextOverlapStart) / 2;

  int ws = tr.windowSize();
  int i = 0;
  while (i < ws && _refIndices[i] < cd._nextOverlapStart)
    i++;
  cd._nextTargetOverlapStart = i;

  i = 0;
  while (i < ws && _refIndices[i] < cd._nextSpliceStart)
    i++;
  cd._nextTargetSpliceStart = i;

  // cd._families = families;
  // cd._weights = new Weights(families);

  cd._targetGL = SplicedGL(overlapHaps, _gl);

  cd._refSamples = rr.samples();
  cd._nRefSamples = cd._refSamples.nSamples();
  cd._targetSamples = tr.samples();
  cd._allSamples = allSamples(tr, rr);
  cd._markers = rr.markers();
  cd._targetMarkers = tr.markers();

  cd._markerIndex = _refIndices;

  cd._restRefHapPairs = _targetRefHapPairs;
  cd._refSampleHapPairs = _refSampleHapPairs;
  cd._restrictedRefSampleHapPairs = SampleHapPairs(cd._refSamples, _targetRefHapPairs, false);

  // cd._recombRate = recombRate(targetMarkers, genMap, par.mapscale());
}

Samples AllData::allSamples(const TargDataReader &tr, const RefDataReader &rr)
{
  if (!_allSamples.nSamples()) {
    checkSampleOverlap(rr.samples(), tr.samples());

    /*
      Create the all-samples Samples object. Target samples are listed
      first so that sample indices agree with sample indices in target
      data genotype likelihoods.
    */
    Samples refSamps = rr.samples();
    Samples targSamps = tr.samples();
    int nRef = refSamps.nSamples();
    int nTarget = targSamps.nSamples();
    for (int j = 0; j < nTarget; ++j) {
      _allSamples.setSamp(targSamps.idIndex(j));
    }
    for (int j = 0; j < nRef; ++j) {
      _allSamples.setSamp(refSamps.idIndex(j));
    }
  }

  return _allSamples;
}

void AllData::checkSampleOverlap(Samples ref, Samples nonRef)
{
  int nRef = ref.nSamples();
  int nNonRef = nonRef.nSamples();
  int n = nRef + nNonRef;
  QList<int> idIndices;

  for (int j = 0; j < nRef; ++j)
    idIndices.append(ref.idIndex(j));

  for (int j = 0; j < nNonRef; ++j)
    idIndices.append(nonRef.idIndex(j));

  qStableSort(idIndices.begin(), idIndices.end());
  for (int j = 1; j < idIndices.length(); ++j) {
    Q_ASSERT_X(idIndices[j - 1] != idIndices[j],
               "AllData::checkSampleOverlap",
               "Overlap between reference and non-reference samples.");
  }
}

static int checkDataAndReturnMaxIndex(QVector<double> probs)
{
  int maxIndex = 0;
  double sum = 0.0;
  for (int j=0; j<probs.length(); ++j)
  {
    Q_ASSERT_X(probs[j] >= 0,
               "checkDataAndReturnMaxIndex (iointerface.cpp)",
               "probs[j] < 0");

    sum += probs[j];

    if (probs[j] > probs[maxIndex])
      maxIndex = j;
  }

  Q_ASSERT_X(abs(sum - 1.0) <= 1e-5,
             "checkDataAndReturnMaxIndex (iointerface.cpp)",
             "abs(sum - 1.0) > 1e-5");

  return maxIndex;
}

void R2Estimator::clear()
{
  _nGenotypes = 0;
  _sumCall = 0.0;
  _sumSquareCall = 0.0;
  _sumExpected = 0.0;
  _sumExpectedSquare = 0.0;
  _sumSquareExpected = 0.0;
  _sumCallExpected = 0.0;
}

void R2Estimator::addSampleData(QVector<double> doseProbs)
{
  Q_ASSERT_X(doseProbs.length() == 3,
             "R2Estimator::addSampleData",
             "doseProbs.length() != 3");

  ++_nGenotypes;
  int call = checkDataAndReturnMaxIndex(doseProbs);
  double exp = (doseProbs[1] + 2*doseProbs[2]);
  double expSquare = (doseProbs[1] + 4*doseProbs[2]);
  _sumCall += call;
  _sumSquareCall += call*call;
  _sumExpected += exp;
  _sumExpectedSquare += expSquare;
  _sumSquareExpected += (exp*exp);
  _sumCallExpected += (call*exp);
}

double R2Estimator::allelicR2()
{
  double f = 1.0/_nGenotypes;
  double cov = _sumCallExpected - (_sumCall * _sumExpected * f);
  double varBest = _sumSquareCall - (_sumCall * _sumCall * f);
  double varExp = _sumExpectedSquare - (_sumExpected * _sumExpected * f);
  double den = varBest*varExp;
  return (den < MIN_R2_DEN) ? 0.0 : (cov*cov/den);
}

double R2Estimator::doseR2()
{
  double f = 1.0/_nGenotypes;
  double num = _sumSquareExpected - (_sumExpected * _sumExpected * f);
  double den = _sumExpectedSquare - (_sumExpected * _sumExpected * f);

  if (num < 0.0)
    num = 0.0;

  return (den < MIN_R2_DEN) ? 0.0 : (num / den);
}


void ImputeDataWriter::printWindowOutput(const CurrentData &cd,
                                         const SampleHapPairs &targetHapPairs,
                                         const ConstrainedAlleleProbs &alProbs,
                                         const Par &par)
{
  bool markersAreImputed = cd.nTargetMarkers() < cd.nMarkers();
  setIsImputed(cd);

  _start = cd.prevSpliceStart();
  _end = cd.nextSpliceStart();

  Q_ASSERT_X(_start <= _end,
             "ImputeDataWriter::printWindowOutput",
             "start > end");

  _printDS = markersAreImputed;     // Dosage
  _printGP = markersAreImputed  &&  par.gprobs();  // Genotype probs

  printWindowData(alProbs);

  /// if (par.ibd())
  ///   printIbd(cd, ibd);
}

void ImputeDataWriter::setIsImputed(const CurrentData &cd)
{
  _isImputed.fill(false, cd.nMarkers());

  if (cd.nTargetMarkers() < _isImputed.length())
  {
    for (int j=0, n=cd.nTargetMarkers(); j<n; j++)
      _isImputed[cd.markerIndex(j)] = true;
  }
}

void ImputeDataWriter::printWindowData(const ConstrainedAlleleProbs &alProbs)
{
  Q_ASSERT_X(_isImputed.length() == alProbs.nMarkers(),
             "ImputeDataWriter::printWindowData",
             "inconsistent data");

  initializeForWindow(4*alProbs.nSamples());

  for (int _mNum=_start; _mNum < _end; ++_mNum)
  {
    const Marker &marker = alProbs.marker(_mNum);
    resetRec(marker);
    _alProbs1.fill(0.0, marker.nAlleles());
    _alProbs2.fill(0.0, marker.nAlleles());
    for (int sample=0, n=alProbs.nSamples(); sample<n; ++sample)
    {
      for (int j=0; j < _alProbs1.length(); ++j)
      {
        _alProbs1[j] = alProbs.alProb1(_mNum, sample, j);
        _alProbs2[j] = alProbs.alProb2(_mNum, sample, j);
      }
      constructSampleDataForMarker();
    }
    finishAndWriteRec();
  }
}

void ImputeDataWriter::initializeForWindow(int initSize)
{
  _gt3Probs.fill(0.0, 3);

  initializeWindowBuffering(initSize);
}

void ImputeDataWriter::resetRec(const Marker &marker)
{
  _r2Est.clear();

  _marker = marker;
  _allele1 = -1;
  _allele2 = -1;
  _nAlleles = marker.nAlleles();
  _nGenotypes = marker.nGenotypes();
  _gtProbs.fill(0.0, marker.nGenotypes());
  _cumAlleleProbs.fill(0.0, marker.nAlleles());
  _dose.fill(0.0, marker.nAlleles());
}

void ImputeDataWriter::constructSampleDataForMarker()
{
  _gt3Probs.fill(0.0, 3);
  _allele1 = maxIndex(_alProbs1, _nAlleles);
  _allele2 = maxIndex(_alProbs2, _nAlleles);
  _dose[0] = _alProbs1[0] + _alProbs2[0];
  _gtProbs[0] = _alProbs1[0] * _alProbs2[0];
  _gt3Probs[0] = _gtProbs[0];
  int gt = 1;
  for (int a2=1; a2 < _alProbs1.length(); ++a2)
  {
    _dose[a2] = _alProbs1[a2] + _alProbs2[a2];

    for (int a1=0; a1 <= a2; ++a1)
    {
      double gtProb = _alProbs1[a1]*_alProbs2[a2];
      if (a1 != a2)
        gtProb += _alProbs1[a2]*_alProbs2[a1];

      _gtProbs[gt++] = gtProb;
      _gt3Probs[(a1==0) ? 1 : 2] += gtProb;
    }
  }

  for (int j=0; j < _dose.length(); ++j)
    _cumAlleleProbs[j] += _dose[j];

  _r2Est.addSampleData(_gt3Probs);

  appendPhasedVariantData();
}


int ImputeDataWriter::maxIndex(QVector<double> &da, int expLength)
{
  Q_ASSERT_X(da.length() == expLength,
             "ImputeDataWriter::maxIndex",
             "da.length() != expLength");

  int maxIndex = 0;
  double sum = 0;

  for (int j=0; j < da.length(); ++j)
  {
    Q_ASSERT_X(da[j] >= 0.0,
               "ImputeDataWriter::maxIndex",
               "da[j] < 0");

    sum += da[j];

    if (da[j] > da[maxIndex])
      maxIndex = j;
  }

  if (sum != 1.0)
  {
    for (int j=0; j < da.length(); ++j)
      da[j] /= sum;
  }

  return maxIndex;
}
