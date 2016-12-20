// Test-dedicated utility include file.

void storePhasedAlleles(BitSetRefGT &rrec, const QList<int> &r1, const QList<int> &r2)
{
  QVector<int> rv1 = r1.toVector();
  QVector<int> rv2 = r2.toVector();

  rrec.storePhasedAlleles(rv1, rv2);
}

void loadTestDataForRefData4x4(Samples &samplesR, QList<BitSetRefGT> &refEmissions, bool useRef6x4=false)
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
  r0.addAllele("A");
  r0.addAllele("C");

  QList<int> r01; r01.append(1); r01.append(0); r01.append(0); r01.append(1);
  QList<int> r02; r02.append(1); r02.append(1); r02.append(0); r02.append(1);
  storePhasedAlleles(r0, r01, r02);

  refEmissions.append(r0);

  BitSetRefGT r1(samplesR);

  r1.setIdInfo(ChromeIds::getIndex("17"), 22345, "RS22345");
  r1.addAllele("G");
  r1.addAllele("T");

  QList<int> r11; r11.append(0); r11.append(0); r11.append(1); r11.append(0);
  QList<int> r12; r12.append(1); r12.append(1); r12.append(1); r12.append(1);
  storePhasedAlleles(r1, r11, r12);

  refEmissions.append(r1);

  BitSetRefGT r2(samplesR);  // This marker will be inbetween two markers that overlap the target data.

  r2.setIdInfo(ChromeIds::getIndex("17"), 22678, "RS22678");
  r2.addAllele("C");
  r2.addAllele("G");

  QList<int> r21; r21.append(0); r21.append(1); r21.append(1); r21.append(0);
  QList<int> r22; r22.append(1); r22.append(0); r22.append(1); r22.append(1);
  storePhasedAlleles(r2, r21, r22);

  refEmissions.append(r2);

  BitSetRefGT r3(samplesR);

  r3.setIdInfo(ChromeIds::getIndex("17"), 33345, "RS33345");
  r3.addAllele("A");
  r3.addAllele("T");

  QList<int> r31; r31.append(0); r31.append(1); r31.append(0); r31.append(0);
  QList<int> r32; r32.append(1); r32.append(0); r32.append(1); r32.append(1);
  storePhasedAlleles(r3, r31, r32);

  refEmissions.append(r3);

  if(useRef6x4)
  {
    BitSetRefGT r4(samplesR);

    r4.setIdInfo(ChromeIds::getIndex("17"), 52345, "RS52345");
    r4.addAllele("G");
    r4.addAllele("T");

    QList<int> r41; r41.append(0); r41.append(0); r41.append(1); r41.append(1);
    QList<int> r42; r42.append(1); r42.append(1); r42.append(1); r42.append(0);
    storePhasedAlleles(r4, r41, r42);

    refEmissions.append(r4);

    BitSetRefGT r5(samplesR);

    r5.setIdInfo(ChromeIds::getIndex("17"), 63345, "RS63345");
    r5.addAllele("A");
    r5.addAllele("T");

    QList<int> r51; r51.append(0); r51.append(1); r51.append(1); r51.append(0);
    QList<int> r52; r52.append(1); r52.append(0); r52.append(0); r52.append(0);
    storePhasedAlleles(r5, r51, r52);

    refEmissions.append(r5);
  }
}

void storeAlleles(BitSetGT &trec, const QList<int> &t1, const QList<int> &t2, const QList<bool> &arePhased)
{
  QVector<int> tv1 = t1.toVector();
  QVector<int> tv2 = t2.toVector();
  QVector<bool> aph = arePhased.toVector();

  trec.storeAlleles(tv1, tv2, aph);
}

void loadTestDataForTargetData3x3(Samples &samplesT, QList<BitSetGT> &targetEmissions,
                                  int missingVal=-1, bool defaultPhasing=false,
                                  bool read6x3=false, bool read4x3=false, bool read4x3B=false,
                                  bool read4x2=false, bool read4x2C=false)
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
  t0.addAllele("A");
  t0.addAllele("C");

  QList<int> t01; t01.append(1); t01.append(0); t01.append(missingVal);
  QList<int> t02; t02.append(1); t02.append(1); t02.append(0);
  // Have arePhased[1] be "true" on purpose....
  QList<bool> arePhased; arePhased.append(defaultPhasing); arePhased.append(true); arePhased.append(defaultPhasing);
  storeAlleles(t0, t01, t02, arePhased);

  targetEmissions.append(t0);

  if(!read4x2)
  {
    BitSetGT t1(samplesT);

    if(read4x3B  ||  read4x2C)
    {
      t1.setIdInfo(ChromeIds::getIndex("17"), 22678, "RS22678");
      t1.addAllele("C");
      t1.addAllele("G");
    }
    else
    {
      t1.setIdInfo(ChromeIds::getIndex("17"), 22345, "RS22345");
      t1.addAllele("G");
      t1.addAllele("T");
    }

    QList<int> t11; t11.append(0); t11.append(0); t11.append(1);
    QList<int> t12; t12.append(1); t12.append((read4x3B) ? 1 : missingVal); t12.append(1);
    arePhased[0] = true; arePhased[2] = true;
    storeAlleles(t1, t11, t12, arePhased);

    targetEmissions.append(t1);
  }

  if(!read4x2C)
  {
    BitSetGT t2(samplesT);

    t2.setIdInfo(ChromeIds::getIndex("17"), 33345, "RS33345");
    t2.addAllele("A");
    t2.addAllele("T");

    QList<int> t21; t21.append(0); t21.append(1); t21.append(1);
    QList<int> t22; t22.append(missingVal); t22.append(0); t22.append(0);
    arePhased[0] = defaultPhasing; arePhased[2] = defaultPhasing;
    storeAlleles(t2, t21, t22, arePhased);

    targetEmissions.append(t2);
  }

  if(read6x3)
  {
    BitSetGT t3(samplesT);

    t3.setIdInfo(ChromeIds::getIndex("17"), 42345, "RS42345");
    t3.addAllele("A");
    t3.addAllele("C");

    QList<int> t31; t31.append(1); t31.append(0); t31.append(missingVal);
    QList<int> t32; t32.append(1); t32.append(1); t32.append(0);
    // Have arePhased[1] be "true" on purpose....
    QList<bool> arePhased; arePhased.append(defaultPhasing); arePhased.append(true); arePhased.append(defaultPhasing);
    storeAlleles(t3, t31, t32, arePhased);

    targetEmissions.append(t3);

    BitSetGT t4(samplesT);

    t4.setIdInfo(ChromeIds::getIndex("17"), 52345, "RS52345");
    t4.addAllele("G");
    t4.addAllele("T");

    QList<int> t41; t41.append(0); t41.append(0); t41.append(1);
    QList<int> t42; t42.append(1); t42.append(missingVal); t42.append(1);
    arePhased[0] = true; arePhased[2] = true;
    storeAlleles(t4, t41, t42, arePhased);

    targetEmissions.append(t4);
  }

  if(read4x3  ||  read4x3B  ||  read6x3)
  {
    BitSetGT t5(samplesT);

    t5.setIdInfo(ChromeIds::getIndex("17"), 63345, "RS63345");
    t5.addAllele("A");
    t5.addAllele("T");

    QList<int> t51; t51.append(0); t51.append(1); t51.append(1);
    QList<int> t52; t52.append(missingVal); t52.append(0); t52.append(0);
    arePhased[0] = defaultPhasing; arePhased[2] = defaultPhasing;
    storeAlleles(t5, t51, t52, arePhased);

    targetEmissions.append(t5);
  }
}


class TestParW : public Par
{
public:
  int window() const {return 4;}
  int overlap() const {return 2;}
  int nThreads() const {return 1;}
  int nSamplingsPerIndividual() const { return 4; }
  int burnin_its() const { return 4; }
  int phase40_its() const { return 4; }
  int niterations() const { return 0; }
};


class RefDataReaderTest4x4 : public RefDataReader
{
public:
  RefDataReaderTest4x4(bool use6x4=false) : _position(0)
  {
    loadTestDataForRefData4x4(_samples, _refEmissionsData, use6x4);
    _nRecs = _refEmissionsData.length();
  }

  bool canAdvanceWindow() const {return _position < _nRecs;}
  bool hasNextRec() const {return _position < _nRecs;}
  BitSetRefGT nextRec() const {return _refEmissionsData[_position];}
  void advanceRec() {_position++;}

  // bool lastWindowOnChrom() const;  // Use default implementation here.

private:
  QList<BitSetRefGT> _refEmissionsData;
  
  int _position;
  int _nRecs;
};

static bool outputDumps = false;

void tdrDumpUtility(const QList<BitSetGT> &vma)
{
  if(outputDumps)
  {
    qDebug("");
    qDebug("Target Data Dump:");
    QString markers("Markers: ");
    for(int m=0; m < vma.length(); m++)
      markers.append("  " + vma[m].marker().id());
    qDebug("%s", (const char *) markers.toLatin1());

    Samples samples = vma[0].samples();
    int nSamples = samples.nSamples();
    for(int sampnum=0; sampnum < nSamples; sampnum++)
    {
      QString haps( QString("For Sample %1:").arg((QString) samples.name(sampnum)) );

      for(int m=0; m < vma.length(); m++)
      {
        int al1 = vma[m].allele1(sampnum);
        int al2 = vma[m].allele2(sampnum);
        QString al1s((al1 >= 0) ? QString("%1").arg(al1) : "?");
        QString al2s((al2 >= 0) ? QString("%1").arg(al2) : "?");
        if (vma[m].isPhased(sampnum))
          haps.append(QString("   %1 | %2").arg(al1s).arg(al2s));
        else
          haps.append(QString("   %1 _ %2").arg(al1s).arg(al2s));
      }
      qDebug("%s", (const char *) haps.toLatin1());
    }
    qDebug("");
  }
}

class TargDataReaderTest3x3 : public TargDataReader
{
public:
 TargDataReaderTest3x3(bool refReady, bool use4x3=false, bool use4x3B=false,
                       bool use4x2=false, bool use4x2C=false) : _position(0)
  {
    if(refReady)
      loadTestDataForTargetData3x3(_samples, _targetEmissionsData, 1, true);
    else if(use4x3)
      loadTestDataForTargetData3x3(_samples, _targetEmissionsData, -1, false, false, true);
    else if(use4x3B)
      loadTestDataForTargetData3x3(_samples, _targetEmissionsData, -1, false, false, false, true);
    else if(use4x2)
      loadTestDataForTargetData3x3(_samples, _targetEmissionsData, -1, false, false, false, false, true);
    else if(use4x2C)
      loadTestDataForTargetData3x3(_samples, _targetEmissionsData, -1, false, false, false, false, false, true);
    else
      loadTestDataForTargetData3x3(_samples, _targetEmissionsData);

    _nRecs = _targetEmissionsData.length();
  }

  bool canAdvanceWindow() const {return _position < _nRecs;}
  bool hasNextRec() const {return _position < _nRecs;}
  BitSetGT nextRec() const {return _targetEmissionsData[_position];}
  void advanceRec() {_position++;}

  // bool lastWindowOnChrom() const;  // Use default implementation here.

  void tdrDump() {tdrDumpUtility(_targetEmissionsData);}

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

  // bool lastWindowOnChrom() const;  // Use default implementation here.

  void tdrDump() {tdrDumpUtility(_targetEmissionsData);}

private:
  QList<BitSetGT> _targetEmissionsData;

  int _position;
  int _nRecs;
};

void hpDump(const HapPairs &hp)
{
  if(outputDumps)
  {
    QString markers("Markers: ");
    for(int m=0; m < hp.nMarkers(); m++)
      markers.append("  " + hp.marker(m).id());
    qDebug("%s", (const char *) markers.toLatin1());

    for(int hpnum=0; hpnum < hp.nHapPairs(); hpnum++)
    {
      QString haps( QString("For Sample %1:").arg((QString) hp.samples(hpnum).name(hp.sampleIndex(hpnum))) );

      for(int m=0; m < hp.nMarkers(); m++)
        haps.append(QString("   %1 | %2").arg(hp.allele1(m, hpnum)).arg(hp.allele2(m, hpnum)));

      qDebug("%s", (const char *) haps.toLatin1());
    }
    qDebug("");
  }
}

void shpDump(const SampleHapPairs &shp)
{
  if(outputDumps)
  {
    qDebug("");
    qDebug("SampleHapPairs:");
    hpDump(shp);
  }
}

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


int testWindowDriverHelper(SampleHapPairs &overlapHaps, const CurrentData &cd, const Par &par, const SampleHapPairs &targetHapPairs)
{
  /*
  if (gv!=null)
    windowOut.printGV(cd, par, gv);        // (Except we won't be implementing GenotypeValues at this time.)
  else
  {
    Map<IntPair, List<IbdSegment>> ibd = mh.refinedIbd(cd, par, targetHapPairs);   // (Nor IBD.)
    AlleleProbs alProbs = mh.LSImpute(cd, par, targetHapPairs);
    printOutput(cd, par, targetHapPairs, alProbs, ibd);  // (We would need to add an output method parameter to this test setup....)
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
      overlap = testWindowDriverHelper(overlapHaps, cd, par, GLSampleHapPairs(cd.targetGL(), true));
    else
      // "Politically incorrect" usage of GLSampleHapPairs....
      overlap = testWindowDriverHelper(overlapHaps, cd, par, GLSampleHapPairs(cd.targetGL(), false));
    // else
    // {
    //   // QList<HapPair> hapPairs = ImputeDriver::phase(cd, par);
    //   QList<HapPair> hapPairs; // For now, default to no haps, just to allow compilation.
    //   overlap = testWindowDriverHelper(overlapHaps, cd, par, SampleHapPairs(cd.targetSamples(), hapPairs));
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

class TestDataWriter : public ImputeDataWriter
{
public:
  TestDataWriter(const Samples &samples) : ImputeDataWriter(samples) {}

  void writeHeader();
  void writeEOF() {}

protected:
  void initializeWindowBuffering(const int initSize, const int nMarkers);
  void appendPhasedVariantData();
  void finishAndWriteRec();

private:
  void outputMarker();
  void outputInfo();
  void outputFormat();

  QString _variantRec;
  QString _outRec;
};

void TestDataWriter::writeHeader()
{
  _outRec = "Samples:";

  for(int s=0, n=_samples.nSamples(); s<n; s++)
    _outRec.append( QString("  %1").arg((QString) _samples.name(s)) );

  qDebug("%s", (const char *) _outRec.toLatin1());
}

void TestDataWriter::initializeWindowBuffering(const int initSize, int nMarkers)
{
  _variantRec.clear();
}

void TestDataWriter::appendPhasedVariantData()
{
  _variantRec.append( QString("  %1|%2").arg(_allele1).arg(_allele2) );

  if (_printDS)
  {
    for (int j=1; j < _nAlleles; ++j)
    {
      _variantRec.append( QString("%1%2").arg((j==1) ? ":" : "," ).arg(_dose[j], 4, 'f', 2) );
    }
  }

  if (_printGP)
  {
    for (int j=0; j < _gtProbs.length(); ++j)
    {
      _variantRec.append( QString("%1%2").arg((j==0) ? ":" : "," ).arg(_gtProbs[j], 4, 'f', 2) );
    }
  }
}

void TestDataWriter::finishAndWriteRec()
{
  outputMarker();
  _outRec.append("  ?  PASS  ");   // QUAL  FILTER
  outputInfo();                    // INFO
  outputFormat();                  // FORMAT

  qDebug("%s  %s", (const char *) _outRec.toLatin1(), (const char *) _variantRec.toLatin1());

  _variantRec.clear();
}

void TestDataWriter::outputMarker()
{
  _outRec = QString("%1 %2 %3 %4 %5")
    .arg( (QString) _marker.chrom() )
    .arg( _marker.pos() )
    .arg( (QString) _marker.id() )
    .arg( (QString) _marker.allele(0) )
    .arg( (_marker.nAlleles() == 1) ? "?" : (QString) _marker.allele(1) );
}

void TestDataWriter::outputInfo()
{
  if (_printDS || _printGP)
  {
    _outRec.append( QString("AR2=%1;DR2=%2")
                    .arg(_r2Est.allelicR2(), 4, 'f', 2)
                    .arg(_r2Est.doseR2(), 4, 'f', 2) );

    for (int j=1; j < _nAlleles; ++j)
    {
      _outRec.append( QString("%1%2")
                      .arg((j==1) ? ";AF=" : "," )
                      .arg(_cumAlleleProbs[j]/(2*_r2Est.nGenotypes()), 4, 'f', 2) );
    }

    if (_isImputed[_mNum])
      _outRec.append(";IMP");
  }
  else
    _outRec.append("?");
}

void TestDataWriter::outputFormat()
{
  if (_printDS)
    _outRec.append( QString("  %1").arg((_printGP) ? "GT:DS:GP" : "GT:DS") );
  else
    _outRec.append("  GT");
}

/*
QList<HapPair> testPhase(CurrentData &cd, const Par &par)
{
  QList<HapPair> hapPairs = ImputeDriver::initialHaps(cd, par);

  //// Phasing happens....then "consensus phasing" happens....

  // HapPairs hp1(hapPairs, false);
  // hpDump(hp1);

  if (par.burnin_its()>0)
  {
    // runStats.println(Const.nl + "Starting burn-in iterations");
    hapPairs = ImputeDriver::runBurnin1(cd, par, hapPairs);

    // HapPairs hp2(hapPairs, false);
    // hpDump(hp2);
  }

  if (par.phase40_its()>0)
  {
    hapPairs = ImputeDriver::runBurnin2(cd, par, hapPairs);

    // HapPairs hp3(hapPairs, false);
    // hpDump(hp3);
  }

  if (par.niterations()>0)
  {
    // runStats.println(Const.nl + "Starting phasing iterations");
    // hapPairs = ImputeDriver::runRecomb(cd, par, hapPairs, gv);
  }
  else
    hapPairs = ConsensusPhaser::consensusPhase(hapPairs);

  HapPairs hp4(hapPairs, false);
  hpDump(hp4);

  return hapPairs;
}

int testPhaseDriverHelper(SampleHapPairs &overlapHaps, const CurrentData &cd, const Par &par,
                          ImputeDataWriter &impWriter, const GLSampleHapPairs &targetHapPairs)
{
  impWriter.printWindowOutput(cd, targetHapPairs, ConstrainedGLAlleleProbs(targetHapPairs), par);
  
  overlapHaps = ImputeDriver::overlapHaps(cd, targetHapPairs);
  return cd.nMarkers() - cd.nextOverlapStart();
}

int testPhaseDriverHelper(SampleHapPairs &overlapHaps, const CurrentData &cd, const Par &par,
                          ImputeDataWriter &impWriter, const SampleHapPairs &targetHapPairs)
{
  // Neither a "printGV" method nor a "refinedIbd" method is invoked here at this time.
  // (COULD RE-COMBINE THE printOutput CALLS (ALWAYS CALLING LSImpute()), SINCE WE ALWAYS NOW USE ConstrainedAlleleProbs.)

  if (cd.nMarkers()==cd.nTargetMarkers() || par.impute() == false)
    impWriter.printWindowOutput(cd, targetHapPairs, ConstrainedAlleleProbs(targetHapPairs), par);
  else
  {
    impWriter.printWindowOutput(cd, targetHapPairs, ImputeDriver::LSImpute(cd, par, targetHapPairs), par);
  }
  
  overlapHaps = ImputeDriver::overlapHaps(cd, targetHapPairs);
  return cd.nMarkers() - cd.nextOverlapStart();
}

void testPhaseDriver(InputData &data, TargDataReader &targReader, RefDataReader &refReader,
                     ImputeDataWriter &impWriter, int windowSize, const Par &par)
{
  CurrentData cd;
  SampleHapPairs overlapHaps;
  int overlap = 0;
  impWriter.writeHeader();
  while(data.canAdvanceWindow(targReader, refReader))
  {
    data.advanceWindow(overlap, par.window(), targReader, refReader);
    data.setCdData(cd, par, overlapHaps, targReader, refReader);

    if (cd.targetGL().isRefData())
      overlap = testPhaseDriverHelper(overlapHaps, cd, par, impWriter, GLSampleHapPairs(cd.targetGL(), true));
    else
    {
      QList<HapPair> hapPairs = testPhase(cd, par);  //   QList<HapPair> hapPairs = ImputeDriver::phase(cd, par);
      overlap = testPhaseDriverHelper(overlapHaps, cd, par, impWriter, SampleHapPairs(cd.targetSamples(), hapPairs, false));
    }
  }
  impWriter.writeEOF();
  } */

