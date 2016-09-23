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

void loadTestDataForTargetData3x3(Samples &samplesT, QList<BitSetGT> &targetEmissions, int missingVal=-1, bool defaultPhasing=false)
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

  t0.setIdInfo(ChromeIds::getIndex((missingVal == 1) ? "1" : "17"), 12345, "RS12345");
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
}

class GLUser     // For testing when a class "has a" FuzzyGL object.
{
public:
  FuzzyGL myGL;
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

static QList<SampleHapPairs> prevOverlapHapsTestList;
static QList<CurrentData> cdCopiesTestList;
static QList<SampleHapPairs> targPairsTestList;
static QList<SampleHapPairs> postOverlapHapsTestList;
static QList<int> overlapAmountsTestList;

void clearStaticTestLists()
{
  prevOverlapHapsTestList.clear();
  cdCopiesTestList.clear();
  targPairsTestList.clear();
  postOverlapHapsTestList.clear();
  overlapAmountsTestList.clear();
}


int testWindowDriverHelper(SampleHapPairs &overlapHaps, /* const CurrentData &cd, */ const SampleHapPairs &targetHapPairs)
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
  prevOverlapHapsTestList.append(overlapHaps);
  // cdCopiesTestList.append(cd);
  targPairsTestList.append(targetHapPairs);

  // Normal code for driver-helper function:

  // overlapHaps = ImputeDriver::overlapHaps(cd, targetHapPairs);
  return 0; // return cd.nMarkers() - cd.nextOverlapStart();
}

void testWindowDriver(InputData &data, TargDataReader &targReader, RefDataReader &refReader, int windowSize)
{
  // CurrentData cd;
  SampleHapPairs overlapHaps;
  int overlap = 0;
  while(data.canAdvanceWindow(targReader, refReader))
  {
    data.advanceWindow(overlap, windowSize, targReader, refReader);
    // data.setCdData(cd, overlapHaps, targReader, refReader);

    // if (cd.targetGL().isRefData())
    //   overlap = testWindowDriverHelper(overlapHaps, cd, GLSampleHapPairs(cd.targetGL()));
    // else
    // {
    //   // QList<HapPair> hapPairs = ImputeDriver::phase(cd);
    //   QList<HapPair> hapPairs; // For now, default to no haps, just to allow compilation.
    //   overlap = testWindowDriverHelper(overlapHaps, cd, SampleHapPairs(cd.targetSamples(), hapPairs));
    // }

    // Testing....
    postOverlapHapsTestList.append(overlapHaps);
    overlapAmountsTestList.append(overlap);

  }
}

