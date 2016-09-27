// Test-dedicated utility include file.

void loadTestDataForRefData4x4(Samples &samplesR, QList<BitSetRefGT> &refEmissions)
{
  // Set the data for "samplesR" before there are any other
  // references to the object.

  samplesR.setSamp(SampleNames::getIndex("SAMP001"));  // The first three names and global sample
  samplesR.setSamp(SampleNames::getIndex("SAMP003"));  // indexes already exist (globally).
  samplesR.setSamp(SampleNames::getIndex("SAMP005"));  // The "SAMP007" name will be new the first
  samplesR.setSamp(SampleNames::getIndex("SAMP007"));  // time we run this utility.

  // Build some VcfEmissions objects as "reference data" (containing
  // phased genotypes without missing data). Be careful to load the
  // data before doing any "copy" or "assignment" operations (which
  // operations are actually reference copy operations).

  BitSetRefGT r0(samplesR);

  r0.setIdInfo(ChromeIds::getIndex("17"), 12345, "RS12345");
  r0.setAllele("A");
  r0.setAllele("C");

  QList<int> r01; r01.append(1); r01.append(0); r01.append(0); r01.append(1);
  QList<int> r02; r02.append(1); r02.append(1); r02.append(0); r02.append(1);
  r0.storePhasedAlleles(r01, r02);

  refEmissions.append(r0);

  BitSetRefGT r1(samplesR);

  r1.setIdInfo(ChromeIds::getIndex("17"), 22345, "RS22345");
  r1.setAllele("G");
  r1.setAllele("T");

  QList<int> r11; r11.append(0); r11.append(0); r11.append(1); r11.append(0);
  QList<int> r12; r12.append(1); r12.append(1); r12.append(1); r12.append(1);
  r1.storePhasedAlleles(r11, r12);

  refEmissions.append(r1);

  BitSetRefGT r2(samplesR);  // This marker will be inbetween two markers that overlap the target data.

  r2.setIdInfo(ChromeIds::getIndex("17"), 22678, "RS22678");
  r2.setAllele("C");
  r2.setAllele("G");

  QList<int> r21; r21.append(0); r21.append(1); r21.append(1); r21.append(0);
  QList<int> r22; r22.append(1); r22.append(0); r22.append(1); r22.append(1);
  r2.storePhasedAlleles(r21, r22);

  refEmissions.append(r2);

  BitSetRefGT r3(samplesR);

  r3.setIdInfo(ChromeIds::getIndex("17"), 33345, "RS33345");
  r3.setAllele("A");
  r3.setAllele("T");

  QList<int> r31; r31.append(0); r31.append(1); r31.append(0); r31.append(0);
  QList<int> r32; r32.append(1); r32.append(0); r32.append(1); r32.append(1);
  r3.storePhasedAlleles(r31, r32);

  refEmissions.append(r3);
}

void loadTestDataForTargetData3x3(Samples &samplesT, QList<BitSetGT> &targetEmissions, int missingVal=-1, bool defaultPhasing=false, bool read6x3=false)
{
  // Set the data for "samplesT" before there are any other
  // references to the object.

  samplesT.setSamp(SampleNames::getIndex("SAMP073"));
  samplesT.setSamp(SampleNames::getIndex("SAMP087"));
  samplesT.setSamp(SampleNames::getIndex("SAMP095"));

  // Build some VcfEmissions objects as "target data" (containing
  // normal unphased genotypes). Be careful to load the data before
  // doing any "copy" or "assignment" operations (which are actually
  // reference copy operations).

  BitSetGT t0(samplesT);

  t0.setIdInfo(ChromeIds::getIndex((defaultPhasing != read6x3) ? "1" : "17"), 12345, "RS12345");
  t0.setAllele("A");
  t0.setAllele("C");

  QList<int> t01; t01.append(1); t01.append(0); t01.append(missingVal);
  QList<int> t02; t02.append(1); t02.append(1); t02.append(0);
  // Have arePhased[1] be "true" on purpose....
  QList<bool> arePhased; arePhased.append(defaultPhasing); arePhased.append(true); arePhased.append(defaultPhasing);
  t0.storeAlleles(t01, t02, arePhased);

  targetEmissions.append(t0);

  BitSetGT t1(samplesT);

  t1.setIdInfo(ChromeIds::getIndex("17"), 22345, "RS22345");
  t1.setAllele("G");
  t1.setAllele("T");

  QList<int> t11; t11.append(0); t11.append(0); t11.append(1);
  QList<int> t12; t12.append(1); t12.append(missingVal); t12.append(1);
  arePhased[0] = true; arePhased[2] = true;
  t1.storeAlleles(t11, t12, arePhased);

  targetEmissions.append(t1);

  BitSetGT t2(samplesT);

  t2.setIdInfo(ChromeIds::getIndex("17"), 33345, "RS33345");
  t2.setAllele("A");
  t2.setAllele("T");

  QList<int> t21; t21.append(0); t21.append(1); t21.append(1);
  QList<int> t22; t22.append(missingVal); t22.append(0); t22.append(0);
  arePhased[0] = defaultPhasing; arePhased[2] = defaultPhasing;
  t2.storeAlleles(t21, t22, arePhased);

  targetEmissions.append(t2);

  if(read6x3)
  {
    BitSetGT t3(samplesT);

    t3.setIdInfo(ChromeIds::getIndex("17"), 42345, "RS42345");
    t3.setAllele("A");
    t3.setAllele("C");

    QList<int> t31; t31.append(1); t31.append(0); t31.append(missingVal);
    QList<int> t32; t32.append(1); t32.append(1); t32.append(0);
    // Have arePhased[1] be "true" on purpose....
    QList<bool> arePhased; arePhased.append(defaultPhasing); arePhased.append(true); arePhased.append(defaultPhasing);
    t3.storeAlleles(t31, t32, arePhased);

    targetEmissions.append(t3);

    BitSetGT t4(samplesT);

    t4.setIdInfo(ChromeIds::getIndex("17"), 52345, "RS52345");
    t4.setAllele("G");
    t4.setAllele("T");

    QList<int> t41; t41.append(0); t41.append(0); t41.append(1);
    QList<int> t42; t42.append(1); t42.append(missingVal); t42.append(1);
    arePhased[0] = true; arePhased[2] = true;
    t4.storeAlleles(t41, t42, arePhased);

    targetEmissions.append(t4);

    BitSetGT t5(samplesT);

    t5.setIdInfo(ChromeIds::getIndex("17"), 63345, "RS63345");
    t5.setAllele("A");
    t5.setAllele("T");

    QList<int> t51; t51.append(0); t51.append(1); t51.append(1);
    QList<int> t52; t52.append(missingVal); t52.append(0); t52.append(0);
    arePhased[0] = defaultPhasing; arePhased[2] = defaultPhasing;
    t5.storeAlleles(t51, t52, arePhased);

    targetEmissions.append(t5);
  }
}

class GLUser     // For testing when a class "has a" FuzzyGL object.
{
public:
  FuzzyGL myGL;
};


class TestParW : public Par
{
public:
  int window() const {return 4;}
  int overlap() const {return 2;}
};


class RefDataReaderTest4x4 : public RefDataReader
{
public:
  RefDataReaderTest4x4() : _position(0)
  {
    loadTestDataForRefData4x4(_samples, _refEmissionsData);
    _nRecs = _refEmissionsData.length();
  }

  bool canAdvanceWindow() const {return _position < _nRecs;}
  bool hasNextRec() const {return _position < _nRecs;}
  BitSetRefGT nextRec() const {return _refEmissionsData[_position];}
  void advanceRec() {_position++;}

private:
  QList<BitSetRefGT> _refEmissionsData;
  
  int _position;
  int _nRecs;
};

class TargDataReaderTest3x3 : public TargDataReader
{
public:
  TargDataReaderTest3x3(bool refReady) : _position(0)
  {
    if(refReady)
      loadTestDataForTargetData3x3(_samples, _targetEmissionsData, 1, true);
    else
      loadTestDataForTargetData3x3(_samples, _targetEmissionsData);

    _nRecs = _targetEmissionsData.length();
  }

  bool canAdvanceWindow() const {return _position < _nRecs;}
  bool hasNextRec() const {return _position < _nRecs;}
  BitSetGT nextRec() const {return _targetEmissionsData[_position];}
  void advanceRec() {_position++;}

private:
  QList<BitSetGT> _targetEmissionsData;

  int _position;
  int _nRecs;
};

class RefTargDataReaderTest6x3 : public TargDataReader
{
public:
  RefTargDataReaderTest6x3() : _position(0)
  {
    loadTestDataForTargetData3x3(_samples, _targetEmissionsData, 1, true, true);
    _nRecs = _targetEmissionsData.length();
  }

  bool canAdvanceWindow() const {return _position < _nRecs;}
  bool hasNextRec() const {return _position < _nRecs;}
  BitSetGT nextRec() const {return _targetEmissionsData[_position];}
  void advanceRec() {_position++;}

private:
  QList<BitSetGT> _targetEmissionsData;

  int _position;
  int _nRecs;
};

static QList<int> prevShpOverlapHapsTestList;
static QList<int> prevHpOverlapHapsTestList;
static QList<CurrentData> cdCopiesTestList;
static QList<int> shptargPairsTestList;
static QList<int> hptargPairsTestList;
static QList<int> postShpOverlapHapsTestList;
static QList<int> postHpOverlapHapsTestList;
static QList<int> overlapAmountsTestList;

void clearStaticTestLists()
{
  prevShpOverlapHapsTestList.clear();
  prevHpOverlapHapsTestList.clear();
  cdCopiesTestList.clear();
  shptargPairsTestList.clear();
  hptargPairsTestList.clear();
  postShpOverlapHapsTestList.clear();
  postHpOverlapHapsTestList.clear();
  overlapAmountsTestList.clear();
}


int testWindowDriverHelper(SampleHapPairs &overlapHaps, const CurrentData &cd, const SampleHapPairs &targetHapPairs)
{
  /*
  if (gv!=null)
    windowOut.printGV(cd, gv);        // (Except we won't be implementing GenotypeValues at this time.)
  else
  {
    Map<IntPair, List<IbdSegment>> ibd = mh.refinedIbd(cd, targetHapPairs);   // (Nor IBD.)
    AlleleProbs alProbs = mh.LSImpute(cd, targetHapPairs);
    printOutput(cd, targetHapPairs, alProbs, ibd);  // (We would need to add an output method parameter to this test setup....)
  }
  */

  // Testing....

  const SampleHapPairs& ovlhshp = *(&overlapHaps);
  prevShpOverlapHapsTestList.append(ovlhshp.nMarkers());
  for(int i=0; i < ovlhshp.nMarkers(); i++)
    prevShpOverlapHapsTestList.append(ovlhshp.allele(i, i));
  const HapPairs& ovlhhp = *(&overlapHaps);
  prevHpOverlapHapsTestList.append(ovlhhp.nMarkers());
  for(int i=0; i < ovlhhp.nMarkers(); i++)
    prevHpOverlapHapsTestList.append(ovlhhp.allele(i, i));

  cdCopiesTestList.append(cd);

  const SampleHapPairs& thpshp = *(&targetHapPairs);
  shptargPairsTestList.append(thpshp.nMarkers());
  for(int i=0; i < thpshp.nMarkers(); i++)
    shptargPairsTestList.append(thpshp.allele(i, i));
  const HapPairs& thphp = *(&targetHapPairs);
  hptargPairsTestList.append(thphp.nMarkers());
  for(int i=0; i < thphp.nMarkers(); i++)
    hptargPairsTestList.append(thphp.allele(i, i));

  // Normal code for driver-helper function:

  overlapHaps = ImputeDriver::overlapHaps(cd, targetHapPairs);
  return cd.nMarkers() - cd.nextOverlapStart();
}

void testWindowDriver(InputData &data, TargDataReader &targReader, RefDataReader &refReader, int windowSize, const Par &par)
{
  CurrentData cd;
  SampleHapPairs overlapHaps;
  int overlap = 0;
  while(data.canAdvanceWindow(targReader, refReader))
  {
    data.advanceWindow(overlap, par.window(), targReader, refReader);
    data.setCdData(cd, par, overlapHaps, targReader, refReader);

    if (cd.targetGL().isRefData())
      overlap = testWindowDriverHelper(overlapHaps, cd, GLSampleHapPairs(cd.targetGL(), true));
    // else
    // {
    //   // QList<HapPair> hapPairs = ImputeDriver::phase(cd);
    //   QList<HapPair> hapPairs; // For now, default to no haps, just to allow compilation.
    //   overlap = testWindowDriverHelper(overlapHaps, cd, SampleHapPairs(cd.targetSamples(), hapPairs));
    // }

    // Testing....

    SampleHapPairs& ovlhshp = *(&overlapHaps);
    postShpOverlapHapsTestList.append(ovlhshp.nMarkers());
    for(int i=0; i < ovlhshp.nMarkers(); i++)
      postShpOverlapHapsTestList.append(ovlhshp.allele(i, i));
    HapPairs& ovlhhp = *(&overlapHaps);
    postHpOverlapHapsTestList.append(ovlhhp.nMarkers());
    for(int i=0; i < ovlhhp.nMarkers(); i++)
      postHpOverlapHapsTestList.append(ovlhhp.allele(i, i));

    overlapAmountsTestList.append(overlap);
  }
}

