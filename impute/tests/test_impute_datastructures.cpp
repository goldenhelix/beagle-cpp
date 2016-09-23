/* Copyright 2016 Golden Helix, Inc. */

#include "impute/haplotypepair.h"
#include "impute/iointerface.h"
#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"

#include "impute/tests/haptests.h"

#include <QObject>
#include <QtTest/QtTest>

class TestImputeDataStructures : public QObject
{
  Q_OBJECT;
private slots:

  void testSamples();
  void testMarkers();
  void testVcfEmissions();
  void testHaplotypePairs();
  void testRefTargetData3x3();
  void testTargetData3x3();
  void testAllData4x4and3x3();
};

void TestImputeDataStructures::testSamples()
{
  int n1 = SampleNames::getIndex("SAMP001");
  int n2 = SampleNames::getIndex("SAMP003");
  int n3 = SampleNames::getIndex("SAMP005");
  int n4 = SampleNames::getIndexIfIndexed("SAMP003");
  int n5 = SampleNames::getIndexIfIndexed("SAMP072");

  QCOMPARE(SampleNames::nNames(), 3);
  QCOMPARE(n1, 0);
  QCOMPARE(n2, 1);
  QCOMPARE(n3, 2);
  QCOMPARE(n4, 1);
  QCOMPARE(n5, -1);
  QCOMPARE(SampleNames::name(0).constData(), "SAMP001");
  QCOMPARE(SampleNames::name(1).constData(), "SAMP003");
  QCOMPARE(SampleNames::name(2).constData(), "SAMP005");

  Samples samples;
  samples.setSamp(1);
  samples.setSamp(2);

  QCOMPARE(samples.nSamples(), 2);
  QCOMPARE(samples.idIndex(0), 1);
  QCOMPARE(samples.idIndex(1), 2);
  QCOMPARE(samples.findLocalIndex(1), 0);
  QCOMPARE(samples.findLocalIndex(2), 1);
  QCOMPARE(samples.findLocalIndex(0), -1);
  QCOMPARE(samples.name(0).constData(), "SAMP003");
  QCOMPARE(samples.name(1).constData(), "SAMP005");

  Samples samples2;
  samples2.setSamp(SampleNames::getIndex("SAMP073"));
  samples2.setSamp(SampleNames::getIndex("SAMP085"));

  QCOMPARE(SampleNames::name(3).constData(), "SAMP073");
  QCOMPARE(SampleNames::name(4).constData(), "SAMP085");

  QCOMPARE(samples2.nSamples(), 2);
  QCOMPARE(samples2.idIndex(0), 3);
  QCOMPARE(samples2.idIndex(1), 4);
  QCOMPARE(samples2.findLocalIndex(4), 1);
  QCOMPARE(samples2.findLocalIndex(3), 0);
  QCOMPARE(samples2.findLocalIndex(2), -1);
  QCOMPARE(samples2.name(0).constData(), "SAMP073");
  QCOMPARE(samples2.name(1).constData(), "SAMP085");

  Samples samples3 = samples2;

  QCOMPARE(samples2.nSamples(), 2);
  QCOMPARE(samples2.idIndex(0), 3);
  QCOMPARE(samples2.idIndex(1), 4);
  QCOMPARE(samples2.findLocalIndex(4), 1);
  QCOMPARE(samples2.findLocalIndex(3), 0);
  QCOMPARE(samples2.findLocalIndex(2), -1);
  QCOMPARE(samples2.name(0).constData(), "SAMP073");
  QCOMPARE(samples2.name(1).constData(), "SAMP085");
}

void TestImputeDataStructures::testMarkers()
{
  int c1 = ChromeIds::getIndex("1");
  int c17 = ChromeIds::getIndex("17");
  int cx = ChromeIds::getIndex("X");
  int c1a = ChromeIds::getIndexIfIndexed("1");
  int c2 = ChromeIds::getIndexIfIndexed("2");
  QCOMPARE(ChromeIds::nIds(), 3);
  QCOMPARE(ChromeIds::chromeId(0).constData(), "1");
  QCOMPARE(ChromeIds::chromeId(1).constData(), "17");
  QCOMPARE(ChromeIds::chromeId(2).constData(), "X");
  QCOMPARE(c1, 0);
  QCOMPARE(c17, 1);
  QCOMPARE(cx, 2);
  QCOMPARE(c1a, 0);
  QCOMPARE(c2, -1);

  // To use objects of a class (such as Marker) in a QList, etc., the
  // class must have:
  //
  // * a default constructor,
  // * a copy constructor, and
  // * an assignment operator.

  Marker m;
  m.setIdInfo(c1, 7632, "RS72351");
  m.setAllele("A");
  m.setAllele("C");

  QCOMPARE(m.chrom().constData(), "1");
  QCOMPARE(m.chromIndex(), 0);
  QCOMPARE(m.pos(), 7632);
  QCOMPARE(m.nAlleles(), 2);
  QCOMPARE(m.allele(0).constData(), "A");
  QCOMPARE(m.allele(1).constData(), "C");

  Marker m2;
  m2.setIdInfo(c1, 7632, "RS72351");
  m2.setAllele("A");
  m2.setAllele("C");

  QCOMPARE(m == m2, true);

  Marker m3;
  m3.setIdInfo(c17, 5432, "RS89351");
  m3.setAllele("G");
  m3.setAllele("T");

  QCOMPARE(m3.chrom().constData(), "17");
  QCOMPARE(m3.pos(), 5432);
  QCOMPARE(m3.nAlleles(), 2);
  QCOMPARE(m3.allele(0).constData(), "G");
  QCOMPARE(m3.allele(1).constData(), "T");
  QCOMPARE(m == m3, false);

  Marker m4 = m2;
  QCOMPARE(m4 == m2, true);

  Marker m5(m2);
  QCOMPARE(m5.pos(), 7632);

  // But, additionally, we want Marker to be a "Moveable" type of
  // class, with its data remotely located (and managed by Qt) as a
  // "QSharedData" class. (Now, the following test would pass if we
  // just had explicit data and copied it around, but since we're
  // using "QSharedData", let's test that shared data is being used
  // successfully.)

  Marker *m6 = new Marker();
  m6->setIdInfo(cx, 23654, "RS93756");
  m6->setAllele("T");
  m6->setAllele("C");

  Marker *m7 = new Marker;
  *m7 = *m6;  // Assign a "copy of m7" itself, rather than assigning a pointer.
  delete m6;  // The original m6 data should still be around after
              // deleting this "Marker" instance.

  QCOMPARE(m7->chrom().constData(), "X");
  QCOMPARE(m7->pos(), 23654);
  QCOMPARE(m7->nAlleles(), 2);
  QCOMPARE(m7->allele(0).constData(), "T");
  QCOMPARE(m7->allele(1).constData(), "C");

  QList<Marker> lm;
  lm.append(Marker());
  lm.append(Marker());
  lm[0].setIdInfo(c17, 23465, "RS76542");
  lm[0].setAllele("G");
  lm[0].setAllele("C");
  lm[1].setIdInfo(cx, 5462, "RS6541");
  lm[1].setAllele("T");
  lm[1].setAllele("A");

  QCOMPARE(lm[0].chrom().constData(), "17");
  QCOMPARE(lm[0].pos(), 23465);
  QCOMPARE(lm[0].id().constData(), "RS76542");
  QCOMPARE(lm[0].nAlleles(), 2);
  QCOMPARE(lm[0].allele(0).constData(), "G");
  QCOMPARE(lm[0].allele(1).constData(), "C");

  QList<Marker> lmb = lm;

  QCOMPARE(lmb[1].chrom().constData(), "X");
  QCOMPARE(lmb[1].pos(), 5462);
  QCOMPARE(lmb[1].id().constData(), "RS6541");
  QCOMPARE(lmb[1].nAlleles(), 2);
  QCOMPARE(lmb[1].allele(0).constData(), "T");
  QCOMPARE(lmb[1].allele(1).constData(), "A");

  Marker m8;
  m8.setIdInfo(c1, 7635, "RS72351");
  m8.setAllele("A");
  m8.setAllele("C");

  Marker m9;
  m9.setIdInfo(c1, 7632, "RS72355");  // Should not make it any "different" from marker m.
  m9.setAllele("A");
  m9.setAllele("C");

  Marker m10;
  m10.setIdInfo(c1, 7632, "RS72351");
  m10.setAllele("A");

  Marker m11;
  m11.setIdInfo(c1, 7632, "RS72351");
  m11.setAllele("G");
  m11.setAllele("C");

  Marker m12;
  m12.setIdInfo(c1, 7632, "RS72351");
  m12.setAllele("A");
  m12.setAllele("G");

  QCOMPARE(m == m8, false);
  QCOMPARE(m == m9, true);
  QCOMPARE(m == m10, false);
  QCOMPARE(m == m11, false);
  QCOMPARE(m == m12, false);

  QList<Marker> biglist;
  biglist.append(m);
  biglist.append(m8);
  biglist.append(m3);
  biglist.append(lm[0]);
  biglist.append(lm[1]);
  biglist.append(*m7);

  delete m7;

  Markers marksfwd(biglist);
  Markers marksfwcopy(marksfwd, false);
  Markers marksrev(marksfwd, true);

  QCOMPARE(marksfwd.nMarkers(), 6);
  QCOMPARE(marksfwcopy.nMarkers(), 6);
  QCOMPARE(marksrev.nMarkers(), 6);
  QCOMPARE(marksfwd.marker(2).pos(), 5432);
  QCOMPARE(marksfwcopy.marker(2).pos(), 5432);
  QCOMPARE(marksrev.marker(2).pos(), 23465);
  QCOMPARE(marksfwd.markers()[1].id().constData(), "RS72351");
  QCOMPARE(marksfwcopy.markers()[1].id().constData(), "RS72351");
  QCOMPARE(marksrev.markers()[1].id().constData(), "RS6541");

  QCOMPARE(marksfwd.restrict(3, 5).markers() == lm, true);
  QList<Marker> restrictedList = marksfwcopy.restrict(3, 5).markers();
  QCOMPARE(restrictedList.length(), 2);
  QCOMPARE(restrictedList[0] == lm[0], true);
  QCOMPARE(restrictedList[1] == lm[1], true);
  QCOMPARE(restrictedList == lm, true);

  QCOMPARE(marksfwd.sumAlleles(2), 4);
  QCOMPARE(marksfwcopy.sumAlleles(), 12);
  QCOMPARE(marksfwd.sumGenotypes(2), 6);
  QCOMPARE(marksfwcopy.sumGenotypes(), 18);
  QCOMPARE(marksfwd.sumHaplotypeBits(2), 2);
  QCOMPARE(marksfwcopy.sumHaplotypeBits(), 6);

  Markers markers2(marksfwd);
  QCOMPARE(markers2.marker(3) == marksfwd.marker(3), true);

  Markers markers3 = markers2;
  QCOMPARE(markers3.marker(2) == marksfwd.marker(2), true);
}

void TestImputeDataStructures::testVcfEmissions()
{
  Samples samplesR;
  Samples samplesT;
  Samples samplesT2;
  QList<BitSetRefGT> refEmissions;
  QList<BitSetGT> unphasedTargetEmissions;
  QList<BitSetGT> phasedTargetEmissions;

  loadTestDataForRefData4x4(samplesR, refEmissions);

  QCOMPARE(SampleNames::getIndexIfIndexed("SAMP073"), 3);
  QCOMPARE(SampleNames::getIndexIfIndexed("SAMP007"), 5);

  QCOMPARE(ChromeIds::getIndexIfIndexed("X"), 2);   // This chromosome name should exist already.
  QCOMPARE(ChromeIds::getIndexIfIndexed("2"), -1);  // This one should not.

  QCOMPARE(refEmissions[0].allele2(1), 1);
  QCOMPARE(refEmissions[1].allele1(0), 0);
  QCOMPARE(refEmissions[2].allele1(2), 1);
  QCOMPARE(refEmissions[3].allele2(3), 1);

  loadTestDataForTargetData3x3(samplesT, unphasedTargetEmissions);

  QCOMPARE(ChromeIds::getIndexIfIndexed("X"), 2);   // This chromosome name should exist already.
  QCOMPARE(ChromeIds::getIndexIfIndexed("2"), -1);  // This one still should not.

  QCOMPARE(refEmissions.length(), 4);
  QCOMPARE(unphasedTargetEmissions.length(), 3);

  QCOMPARE(unphasedTargetEmissions[1].isPhased(1), true);
  QCOMPARE(unphasedTargetEmissions[0].isRefData(), false);
  QCOMPARE(unphasedTargetEmissions[1].isRefData(), false);
  QCOMPARE(unphasedTargetEmissions[0].allele1(2), -1);
  QCOMPARE(unphasedTargetEmissions[0].allele2(1), 1);
  QCOMPARE(unphasedTargetEmissions[1].allele1(0), 0);
  QCOMPARE(unphasedTargetEmissions[2].allele1(2), 1);
  QCOMPARE(unphasedTargetEmissions[1].allele2(1), -1);
  QCOMPARE(unphasedTargetEmissions[2].allele1(0), 0);
  QCOMPARE(unphasedTargetEmissions[2].allele2(0), -1);
  QCOMPARE(unphasedTargetEmissions[0].gl(0, 1, 1), 1.0);
  QCOMPARE(unphasedTargetEmissions[0].gl(1, 0, 1), 1.0);
  QCOMPARE(unphasedTargetEmissions[0].gl(1, 1, 0), 0.0);
  QCOMPARE(unphasedTargetEmissions[0].gl(2, 1, 0), 1.0);
  QCOMPARE(unphasedTargetEmissions[0].gl(2, 0, 1), 1.0);
  QCOMPARE(unphasedTargetEmissions[0].gl(2, 1, 1), 0.0);
  QCOMPARE(unphasedTargetEmissions[2].gl(1, 1, 0), 1.0);
  QCOMPARE(unphasedTargetEmissions[2].gl(1, 0, 1), 0.0);
  QCOMPARE(unphasedTargetEmissions[2].gl(2, 1, 0), 1.0);
  QCOMPARE(unphasedTargetEmissions[2].gl(2, 0, 1), 1.0);

  loadTestDataForTargetData3x3(samplesT2, phasedTargetEmissions, 1, true);

  QCOMPARE(ChromeIds::getIndexIfIndexed("X"), 2);   // This chromosome name should exist already.
  QCOMPARE(ChromeIds::getIndexIfIndexed("2"), -1);  // This one still should not.

  QCOMPARE(phasedTargetEmissions.length(), 3);

  QCOMPARE(phasedTargetEmissions[1].isPhased(1), true);
  QCOMPARE(phasedTargetEmissions[0].isRefData(), true);
  QCOMPARE(phasedTargetEmissions[1].isRefData(), true);
  QCOMPARE(phasedTargetEmissions[0].allele1(2), 1);
  QCOMPARE(phasedTargetEmissions[0].allele2(1), 1);
  QCOMPARE(phasedTargetEmissions[1].allele1(0), 0);
  QCOMPARE(phasedTargetEmissions[2].allele1(2), 1);
  QCOMPARE(phasedTargetEmissions[1].allele2(1), 1);
  QCOMPARE(phasedTargetEmissions[2].allele1(0), 0);
  QCOMPARE(phasedTargetEmissions[2].allele2(0), 1);
  QCOMPARE(phasedTargetEmissions[0].gl(0, 1, 1), 1.0);
  QCOMPARE(phasedTargetEmissions[0].gl(1, 0, 1), 1.0);
  QCOMPARE(phasedTargetEmissions[0].gl(1, 1, 0), 0.0);
  QCOMPARE(phasedTargetEmissions[0].gl(2, 1, 0), 1.0);
  QCOMPARE(phasedTargetEmissions[0].gl(2, 0, 1), 0.0);
  QCOMPARE(phasedTargetEmissions[0].gl(2, 1, 1), 0.0);
  QCOMPARE(phasedTargetEmissions[2].gl(1, 1, 0), 1.0);
  QCOMPARE(phasedTargetEmissions[2].gl(1, 0, 1), 0.0);
  QCOMPARE(phasedTargetEmissions[2].gl(2, 1, 0), 1.0);
  QCOMPARE(phasedTargetEmissions[2].gl(2, 0, 1), 0.0);
}

void TestImputeDataStructures::testHaplotypePairs()
{
  Samples samplesR;
  QList<BitSetRefGT> refEmissions;
  loadTestDataForRefData4x4(samplesR, refEmissions);

  int nummarks = refEmissions.length();

  QList<Marker> biglist;
  QList<int> als10;
  QList<int> als20;
  QList<int> als11;
  QList<int> als21;
  QList<int> als12;
  QList<int> als22;
  QList<int> als13;
  QList<int> als23;

  for (int mnum = 0; mnum < nummarks; mnum++) {
    biglist.append(refEmissions[mnum].marker());
    als10.append(refEmissions[mnum].allele1(0));
    als20.append(refEmissions[mnum].allele2(0));
    als11.append(refEmissions[mnum].allele1(1));
    als21.append(refEmissions[mnum].allele2(1));
    als12.append(refEmissions[mnum].allele1(2));
    als22.append(refEmissions[mnum].allele2(2));
    als13.append(refEmissions[mnum].allele1(3));
    als23.append(refEmissions[mnum].allele2(3));
  }

  Markers marks(biglist);
  HapPair pair1(marks, samplesR, 1, als11, als21);
  HapPair pair3(marks, samplesR, 3, als13, als23);

  HapPair pair1copy(pair1, false);  // Copy but don't reverse.
  HapPair pair3rev(pair3, true);    // Copy and reverse this one.

  QCOMPARE(pair1.allele1(1), 0);
  QCOMPARE(pair1.allele2(1), 1);
  QCOMPARE(pair1.allele1(2), 1);
  QCOMPARE(pair1.allele2(2), 0);

  QCOMPARE(pair3.allele1(0), 1);
  QCOMPARE(pair3.allele2(0), 1);
  QCOMPARE(pair3.allele1(3), 0);
  QCOMPARE(pair3.allele2(3), 1);

  QCOMPARE(pair1copy.allele1(1), 0);
  QCOMPARE(pair1copy.allele2(1), 1);
  QCOMPARE(pair1copy.allele1(2), 1);
  QCOMPARE(pair1copy.allele2(2), 0);

  QCOMPARE(pair3rev.allele1(0), 0);
  QCOMPARE(pair3rev.allele2(0), 1);
  QCOMPARE(pair3rev.allele1(3), 1);
  QCOMPARE(pair3rev.allele2(3), 1);

  QCOMPARE(pair1.sampleIndex(), 1);
  QCOMPARE(pair3.sampleIndex(), 3);
  QCOMPARE(pair1copy.sampleIndex(), 1);
  QCOMPARE(pair3rev.sampleIndex(), 3);

  QCOMPARE(pair1.samples().name(0).constData(), "SAMP001");
  QCOMPARE(pair3.samples().name(1).constData(), "SAMP003");
  QCOMPARE(pair1copy.samples().name(2).constData(), "SAMP005");
  QCOMPARE(pair3rev.samples().name(3).constData(), "SAMP007");

  QCOMPARE(pair1.markers().marker(1).id().constData(), "RS22345");
  QCOMPARE(pair3.marker(2).pos(), 22678);
  QCOMPARE(pair3rev.marker(2).pos(), 22345);
  QCOMPARE(pair1copy.nMarkers(), 4);

  QList<HapPair> hpsList;
  hpsList.append(pair1);
  hpsList.append(pair3);
  hpsList.append(pair1copy);
  HapPairs hps(hpsList, false);
  HapPairs hpsr(hpsList, true);

  QCOMPARE(hps.allele1(1, 0), 0);
  QCOMPARE(hps.allele2(0, 1), 1);
  QCOMPARE(hps.allele1(3, 2), 1);
  QCOMPARE(hps.allele2(3, 2), 0);

  QCOMPARE(hps.allele2(1, 1), 1);
  QCOMPARE(hps.allele1(1, 2), 0);

  QCOMPARE(hps.allele(1, 3), 1);
  QCOMPARE(hps.allele(1, 4), 0);
  QCOMPARE(hps.nMarkers(), 4);
  Markers mks = hps.markers();
  Marker m2 = hps.marker(2);
  Marker m3 = hps.marker(3);

  QCOMPARE(hps.nHaps(), 6);
  QCOMPARE(hps.nHapPairs(), 3);
  Samples shps1 = hps.samples(1);
  QCOMPARE(hps.sampleIndex(0), 1);
  QCOMPARE(hps.sampleIndex(1), 3);
  QCOMPARE(hps.sampleIndex(2), 1);

  QCOMPARE(mks == marks, true);
  QCOMPARE(m2 == refEmissions[2].marker(), true);
  QCOMPARE(shps1 == refEmissions[1].samples(), true);
  QCOMPARE(shps1 == samplesR, true);

  QCOMPARE(hpsr.allele1(2, 0), 0);
  QCOMPARE(hpsr.allele2(3, 1), 1);
  QCOMPARE(hpsr.allele1(0, 2), 1);
  QCOMPARE(hpsr.allele2(0, 2), 0);

  QCOMPARE(hpsr.allele2(2, 1), 1);
  QCOMPARE(hpsr.allele1(2, 2), 0);

  QCOMPARE(hpsr.allele(2, 3), 1);
  QCOMPARE(hpsr.allele(2, 4), 0);
  QCOMPARE(hpsr.nMarkers(), 4);
  Markers mksr = hpsr.markers();
  Marker m2r = hpsr.marker(2);

  QCOMPARE(hpsr.nHaps(), 6);
  QCOMPARE(hpsr.nHapPairs(), 3);
  Samples shpsr1 = hpsr.samples(1);
  QCOMPARE(hpsr.sampleIndex(0), 1);
  QCOMPARE(hpsr.sampleIndex(1), 3);
  QCOMPARE(hpsr.sampleIndex(2), 1);

  QCOMPARE(mksr == marks, false);
  QCOMPARE(m2r == refEmissions[1].marker(), true);
  QCOMPARE(shpsr1 == refEmissions[1].samples(), true);
  QCOMPARE(shpsr1 == samplesR, true);

  Samples samplesComplete;
  samplesComplete.setSamp(SampleNames::getIndex("SAMP001"));
  samplesComplete.setSamp(SampleNames::getIndex("SAMP003"));
  samplesComplete.setSamp(SampleNames::getIndex("SAMP005"));
  samplesComplete.setSamp(SampleNames::getIndex("SAMP007"));

  // SampleHapPairs shpsWrongLengthSamples(samplesComplete, hpsList, false);  // This should
  // ASSERT-crash, and does.

  // Samples samplesSHPMMS;
  // samplesSHPMMS.setSamp(SampleNames::getIndex("SAMP001"));
  // samplesSHPMMS.setSamp(SampleNames::getIndex("SAMP003"));
  // samplesSHPMMS.setSamp(SampleNames::getIndex("SAMP005"));

  // SampleHapPairs shpsMisMatchSamples(samplesSHPMMS, hpsList, false);  // This should
  // ASSERT-crash, and does.

  HapPair pair0(marks, samplesR, 0, als10, als20);
  HapPair pair2(marks, samplesR, 2, als12, als22);

  QList<HapPair> shpscList;
  shpscList.append(pair0);
  shpscList.append(pair1);
  shpscList.append(pair2);
  shpscList.append(pair3);
  SampleHapPairs shpsc(samplesComplete, shpscList, false);
  SampleHapPairs shpscr(samplesComplete, shpscList, true);

  QCOMPARE(shpsc.allele1(0, 0), 1);
  QCOMPARE(shpsc.allele2(1, 1), 1);
  QCOMPARE(shpsc.allele1(2, 2), 1);
  QCOMPARE(shpsc.allele2(3, 3), 1);
  QCOMPARE(shpsc.allele1(3, 0), 0);
  QCOMPARE(shpsc.allele2(2, 1), 0);
  QCOMPARE(shpsc.allele1(1, 2), 1);
  QCOMPARE(shpsc.allele2(0, 3), 1);
  QCOMPARE(shpscr.allele1(3, 0), 1);
  QCOMPARE(shpscr.allele2(2, 1), 1);
  QCOMPARE(shpscr.allele1(1, 2), 1);
  QCOMPARE(shpscr.allele2(0, 3), 1);
  QCOMPARE(shpscr.allele1(0, 0), 0);
  QCOMPARE(shpscr.allele2(1, 1), 0);
  QCOMPARE(shpscr.allele1(2, 2), 1);
  QCOMPARE(shpscr.allele2(3, 3), 1);

  // The following code will test having HapPair objects associated
  // with different Samples objects work together both in a HapPairs
  // object and in a SampleHapPairs object.

  Samples samplesX;
  samplesX.setSamp(SampleNames::getIndex("SAMP029"));

  HapPair pairX(marks, samplesX, 0, als13, als21);
  QList<HapPair> xList;
  xList.append(pairX);

  HapPairs hpsx(xList, true);  // Test reversing, while we're here.

  QCOMPARE(hpsx.allele2(0, 0), 0);
  QCOMPARE(hpsx.allele1(2, 0), 0);
  QCOMPARE(hpsx.allele1(3, 0), 1);
  QCOMPARE(hpsx.allele2(2, 0), 1);

  Samples samplesCompleteX;
  // Notice the order of samples 005 and 003 is switched from the
  // order in samplesComplete.
  samplesCompleteX.setSamp(SampleNames::getIndex("SAMP001"));
  samplesCompleteX.setSamp(SampleNames::getIndex("SAMP005"));
  samplesCompleteX.setSamp(SampleNames::getIndex("SAMP003"));
  samplesCompleteX.setSamp(SampleNames::getIndex("SAMP007"));
  samplesCompleteX.setSamp(SampleNames::getIndex("SAMP029"));

  shpscList.append(pairX);
  SampleHapPairs shpscx(samplesCompleteX, shpscList, false);
  SampleHapPairs shpscxr(samplesCompleteX, shpscList, true);

  QCOMPARE(shpscx.allele1(0, 2), 0);
  QCOMPARE(shpscx.allele2(1, 1), 1);
  QCOMPARE(shpscx.allele1(2, 3), 0);
  QCOMPARE(shpscx.allele2(3, 0), 1);
  QCOMPARE(shpscx.allele2(2, 4), 0);
  QCOMPARE(shpscxr.allele1(3, 2), 0);
  QCOMPARE(shpscxr.allele2(2, 1), 1);
  QCOMPARE(shpscxr.allele1(1, 3), 0);
  QCOMPARE(shpscxr.allele2(0, 0), 1);
  QCOMPARE(shpscxr.allele2(1, 4), 0);
  QCOMPARE(shpscx.allele1(3, 2), 1);
  QCOMPARE(shpscx.allele2(2, 1), 1);
  QCOMPARE(shpscx.allele1(1, 3), 0);
  QCOMPARE(shpscx.allele2(0, 0), 1);
  QCOMPARE(shpscx.allele2(1, 4), 1);

  QCOMPARE(shpscx.allele(1, 4), 0);
  QCOMPARE(shpscx.allele(1, 9), 1);
  QCOMPARE(shpscx.allele(2, 4), 1);
  QCOMPARE(shpscx.allele(2, 9), 0);
  QCOMPARE(shpscx.nMarkers(), 4);
  Markers mksscx = shpscx.markers();
  Marker m2scx = shpscx.marker(2);

  QCOMPARE(shpscx.nHaps(), 10);
  QCOMPARE(shpscx.nHapPairs(), 5);
  Samples sshpscx4 = shpscx.samples(4);
  Samples sshpscx = shpscx.samples();
  QCOMPARE(shpscx.sampleIndex(0), 0);
  QCOMPARE(shpscx.sampleIndex(1), 2);
  QCOMPARE(shpscx.sampleIndex(2), 1);
  QCOMPARE(shpscx.sampleIndex(3), 3);
  QCOMPARE(shpscx.sampleIndex(4), 0);

  QCOMPARE(mksscx == marks, true);
  QCOMPARE(m2scx == m2, true);
  QCOMPARE(sshpscx4 == samplesX, true);
  QCOMPARE(sshpscx == samplesCompleteX, true);
  QCOMPARE(shpscx.nSamples(), 5);

  QCOMPARE(shpscxr.allele(0, 4), 1);
  QCOMPARE(shpscxr.allele(0, 9), 0);
  QCOMPARE(shpscxr.allele(3, 4), 0);
  QCOMPARE(shpscxr.allele(3, 9), 1);
  QCOMPARE(shpscxr.nMarkers(), 4);
  Markers mksscxr = shpscxr.markers();
  Marker m1scxr = shpscxr.marker(1);

  QCOMPARE(shpscxr.nHaps(), 10);
  QCOMPARE(shpscxr.nHapPairs(), 5);
  Samples sshpscxr4 = shpscxr.samples(4);
  Samples sshpscxr = shpscxr.samples();
  QCOMPARE(shpscxr.sampleIndex(0), 0);
  QCOMPARE(shpscxr.sampleIndex(1), 2);
  QCOMPARE(shpscxr.sampleIndex(2), 1);
  QCOMPARE(shpscxr.sampleIndex(3), 3);
  QCOMPARE(shpscxr.sampleIndex(4), 0);

  QCOMPARE(mksscxr == marks, false);
  QCOMPARE(m1scxr == m2, true);
  QCOMPARE(sshpscxr4 == samplesX, true);
  QCOMPARE(sshpscxr == samplesCompleteX, true);
  QCOMPARE(shpscxr.nSamples(), 5);

  Samples resamples = refEmissions[1].samples();
  QCOMPARE(resamples == samplesR, true);

  RefHapPairs rhp(resamples, refEmissions);
  QCOMPARE(rhp.markers().marker(1) == marks.marker(1), true);
  QCOMPARE(rhp.samples(1) == resamples, true);
  QCOMPARE(rhp.samples(1).idIndex(2), samplesComplete.idIndex(2));
  QCOMPARE(rhp.allele1(1, 2), 1);
  QCOMPARE(rhp.allele2(2, 1), 0);
  QCOMPARE(rhp.allele(3, 5), 1);
  QCOMPARE(rhp.allele(0, 6), 1);
  QCOMPARE(rhp.sampleIndex(2), 2);
  QCOMPARE(rhp.nAlleles(3), 2);

  Samples samplesT;
  QList<BitSetGT> phasedTargetEmissions;
  loadTestDataForTargetData3x3(samplesT, phasedTargetEmissions, 1, true);
  QCOMPARE(phasedTargetEmissions[1].samples() == samplesT, true);

  SplicedGL phasedtarggl(samplesT, phasedTargetEmissions);
  GLSampleHapPairs glshp(phasedtarggl);

  QCOMPARE(phasedtarggl.allele1(1, 1), 0);
  QCOMPARE(phasedtarggl.allele2(1, 1), 1);
  QCOMPARE(phasedtarggl.allele2(2, 0), 1);
  QCOMPARE(phasedtarggl.allele1(0, 2), 1);
  QCOMPARE(phasedtarggl.allele(0, 4), 1);
  QCOMPARE(phasedtarggl.allele(2, 1), 1);
  QCOMPARE(phasedtarggl.nMarkers(), 3);
  QCOMPARE(phasedtarggl.nAlleles(1), 2);
  QCOMPARE(phasedtarggl.isRefData(), true);

  QCOMPARE(phasedtarggl.nHaps(), 6);
  QCOMPARE(phasedtarggl.nHapPairs(), 3);
  Samples sphasedtarggl2 = phasedtarggl.samples(2);
  Samples sphasedtarggl = phasedtarggl.samples();
  QCOMPARE(phasedtarggl.sampleIndex(0), 0);
  QCOMPARE(phasedtarggl.sampleIndex(2), 2);

  QCOMPARE(sphasedtarggl2 == samplesT, true);
  QCOMPARE(phasedtarggl.nSamples(), 3);

  QCOMPARE(phasedtarggl.gl(0, 0, 1, 1), 1.0);
  QCOMPARE(phasedtarggl.gl(0, 0, 1, 0), 0.0);
  QCOMPARE(phasedtarggl.gl(0, 0, 0, 1), 0.0);
  QCOMPARE(phasedtarggl.gl(0, 0, 0, 0), 0.0);
  QCOMPARE(phasedtarggl.gl(0, 1, 1, 1), 0.0);
  QCOMPARE(phasedtarggl.gl(0, 1, 1, 0), 0.0);
  QCOMPARE(phasedtarggl.gl(0, 1, 0, 1), 1.0);
  QCOMPARE(phasedtarggl.gl(0, 1, 0, 0), 0.0);
  QCOMPARE(phasedtarggl.gl(0, 2, 1, 1), 0.0);
  QCOMPARE(phasedtarggl.gl(0, 2, 1, 0), 1.0);
  QCOMPARE(phasedtarggl.gl(0, 2, 0, 1), 0.0);
  QCOMPARE(phasedtarggl.gl(0, 2, 0, 0), 0.0);
  QCOMPARE(phasedtarggl.gl(1, 1, 1, 1), 0.0);
  QCOMPARE(phasedtarggl.gl(1, 1, 1, 0), 0.0);
  QCOMPARE(phasedtarggl.gl(1, 1, 0, 1), 1.0);
  QCOMPARE(phasedtarggl.gl(1, 1, 0, 0), 0.0);
  QCOMPARE(phasedtarggl.gl(2, 0, 1, 1), 0.0);
  QCOMPARE(phasedtarggl.gl(2, 0, 1, 0), 0.0);
  QCOMPARE(phasedtarggl.gl(2, 0, 0, 1), 1.0);
  QCOMPARE(phasedtarggl.gl(2, 0, 0, 0), 0.0);
  QCOMPARE(phasedtarggl.isPhased(0, 0), true);
  QCOMPARE(phasedtarggl.isPhased(0, 1), true);
  QCOMPARE(phasedtarggl.isPhased(0, 2), true);
  QCOMPARE(phasedtarggl.isPhased(1, 0), true);
  QCOMPARE(phasedtarggl.isPhased(1, 1), true);
  QCOMPARE(phasedtarggl.isPhased(1, 2), true);
  QCOMPARE(phasedtarggl.isPhased(2, 0), true);
  QCOMPARE(phasedtarggl.isPhased(2, 1), true);
  QCOMPARE(phasedtarggl.isPhased(2, 2), true);

  QCOMPARE(glshp.allele1(1, 1), 0);
  QCOMPARE(glshp.allele2(1, 1), 1);
  QCOMPARE(glshp.allele2(2, 0), 1);
  QCOMPARE(glshp.allele1(0, 2), 1);
  QCOMPARE(glshp.allele(0, 4), 1);
  QCOMPARE(glshp.allele(2, 1), 1);
  QCOMPARE(glshp.nMarkers(), 3);
  QCOMPARE(glshp.nAlleles(1), 2);
  QCOMPARE(glshp.isRefData(), true);

  QCOMPARE(glshp.nHaps(), 6);
  QCOMPARE(glshp.nHapPairs(), 3);
  Samples sglshp2 = glshp.samples(2);
  Samples sglshp = glshp.samples();
  QCOMPARE(glshp.sampleIndex(0), 0);
  QCOMPARE(glshp.sampleIndex(2), 2);

  QCOMPARE(sglshp2 == samplesT, true);
  QCOMPARE(glshp.nSamples(), 3);

  Samples samplesT2;
  QList<BitSetGT> unphasedTargetEmissions;
  loadTestDataForTargetData3x3(samplesT2, unphasedTargetEmissions);

  SplicedGL unphasedtarggl(samplesT2, unphasedTargetEmissions);

  // This one should not work, and indeed doesn't:
  //   GLSampleHapPairs crash(unphasedtarggl);

  QCOMPARE(unphasedtarggl.allele1(1, 1), 0);
  QCOMPARE(unphasedtarggl.allele2(1, 1), -1);
  QCOMPARE(unphasedtarggl.allele2(2, 0), -1);
  QCOMPARE(unphasedtarggl.allele1(0, 2), -1);
  QCOMPARE(unphasedtarggl.allele(0, 4), -1);
  QCOMPARE(unphasedtarggl.allele(2, 1), -1);
  QCOMPARE(unphasedtarggl.nMarkers(), 3);
  QCOMPARE(unphasedtarggl.nAlleles(1), 2);
  QCOMPARE(unphasedtarggl.isRefData(), false);

  QCOMPARE(unphasedtarggl.nHaps(), 6);
  QCOMPARE(unphasedtarggl.nHapPairs(), 3);
  Samples sunphasedtarggl2 = unphasedtarggl.samples(2);
  Samples sunphasedtarggl = unphasedtarggl.samples();
  QCOMPARE(unphasedtarggl.sampleIndex(0), 0);
  QCOMPARE(unphasedtarggl.sampleIndex(2), 2);

  QCOMPARE(sunphasedtarggl2 == samplesT2, true);
  QCOMPARE(unphasedtarggl.nSamples(), 3);

  QCOMPARE(unphasedtarggl.gl(0, 0, 1, 1), 1.0);
  QCOMPARE(unphasedtarggl.gl(0, 0, 1, 0), 0.0);
  QCOMPARE(unphasedtarggl.gl(0, 0, 0, 1), 0.0);
  QCOMPARE(unphasedtarggl.gl(0, 0, 0, 0), 0.0);
  QCOMPARE(unphasedtarggl.gl(0, 1, 1, 1), 0.0);
  QCOMPARE(unphasedtarggl.gl(0, 1, 1, 0), 0.0);
  QCOMPARE(unphasedtarggl.gl(0, 1, 0, 1), 1.0);
  QCOMPARE(unphasedtarggl.gl(0, 1, 0, 0), 0.0);
  QCOMPARE(unphasedtarggl.gl(0, 2, 1, 1), 0.0);
  QCOMPARE(unphasedtarggl.gl(0, 2, 1, 0), 1.0);
  QCOMPARE(unphasedtarggl.gl(0, 2, 0, 1), 1.0);
  QCOMPARE(unphasedtarggl.gl(0, 2, 0, 0), 1.0);
  QCOMPARE(unphasedtarggl.gl(1, 1, 1, 1), 0.0);
  QCOMPARE(unphasedtarggl.gl(1, 1, 1, 0), 0.0);
  QCOMPARE(unphasedtarggl.gl(1, 1, 0, 1), 1.0);
  QCOMPARE(unphasedtarggl.gl(1, 1, 0, 0), 1.0);
  QCOMPARE(unphasedtarggl.gl(2, 0, 1, 1), 0.0);
  QCOMPARE(unphasedtarggl.gl(2, 0, 1, 0), 1.0);
  QCOMPARE(unphasedtarggl.gl(2, 0, 0, 1), 1.0);
  QCOMPARE(unphasedtarggl.gl(2, 0, 0, 0), 1.0);
  QCOMPARE(unphasedtarggl.gl(2, 1, 1, 1), 0.0);
  QCOMPARE(unphasedtarggl.gl(2, 1, 1, 0), 1.0);
  QCOMPARE(unphasedtarggl.gl(2, 1, 0, 1), 0.0);
  QCOMPARE(unphasedtarggl.gl(2, 1, 0, 0), 0.0);
  QCOMPARE(unphasedtarggl.gl(2, 2, 1, 1), 0.0);
  QCOMPARE(unphasedtarggl.gl(2, 2, 1, 0), 1.0);
  QCOMPARE(unphasedtarggl.gl(2, 2, 0, 1), 1.0);
  QCOMPARE(unphasedtarggl.gl(2, 2, 0, 0), 0.0);
  QCOMPARE(unphasedtarggl.isPhased(0, 0), false);
  QCOMPARE(unphasedtarggl.isPhased(0, 1), true);
  QCOMPARE(unphasedtarggl.isPhased(0, 2), false);
  QCOMPARE(unphasedtarggl.isPhased(1, 0), true);
  QCOMPARE(unphasedtarggl.isPhased(1, 1), true);
  QCOMPARE(unphasedtarggl.isPhased(1, 2), true);
  QCOMPARE(unphasedtarggl.isPhased(2, 0), false);
  QCOMPARE(unphasedtarggl.isPhased(2, 1), true);
  QCOMPARE(unphasedtarggl.isPhased(2, 2), false);

  QList<Marker> smalllist;
  QList<int> alss10;
  QList<int> alss20;
  QList<int> alss11;
  QList<int> alss21;
  QList<int> alss12;
  QList<int> alss22;

  for (int mnum = 0; mnum < 2; mnum++) {
    smalllist.append(refEmissions[mnum].marker());
    alss10.append(phasedTargetEmissions[mnum].allele1(0));
    alss20.append(phasedTargetEmissions[mnum].allele2(0));
    alss11.append(phasedTargetEmissions[mnum].allele1(1));
    alss21.append(phasedTargetEmissions[mnum].allele2(1));
    alss12.append(phasedTargetEmissions[mnum].allele1(2));
    alss22.append(phasedTargetEmissions[mnum].allele2(2));
  }

  Markers markss(smalllist);
  Markers mksptgl = phasedtarggl.markers();
  Marker m2ptgl = phasedtarggl.marker(2);
  QCOMPARE(mksptgl == markss, false);
  QCOMPARE(m2ptgl == m2, false);
  QCOMPARE(m2ptgl == m3, true);
  Markers mksglshp = glshp.markers();
  Marker m2glshp = glshp.marker(2);
  QCOMPARE(mksglshp == markss, false);
  QCOMPARE(m2glshp == m3, true);

  HapPair pair0s(markss, samplesT2, 0, alss10, alss20);
  HapPair pair1s(markss, samplesT2, 1, alss11, alss21);
  HapPair pair2s(markss, samplesT2, 2, alss12, alss22);
  QList<HapPair> pteList;
  pteList.append(pair0s);
  pteList.append(pair1s);
  pteList.append(pair2s);
  SampleHapPairs shppte(samplesT2, pteList, false);
  SplicedGL spliced(shppte, unphasedtarggl);

  QCOMPARE(spliced.allele1(1, 1), 0);
  QCOMPARE(spliced.allele2(1, 1), 1);
  QCOMPARE(spliced.allele2(2, 0), -1);
  QCOMPARE(spliced.allele1(0, 2), 1);
  QCOMPARE(spliced.allele(0, 4), 1);
  QCOMPARE(spliced.allele(2, 1), -1);
  QCOMPARE(spliced.nMarkers(), 3);
  QCOMPARE(spliced.nAlleles(1), 2);
  QCOMPARE(spliced.isRefData(), false);

  QCOMPARE(spliced.nHaps(), 6);
  QCOMPARE(spliced.nHapPairs(), 3);
  Samples sspliced2 = spliced.samples(2);
  Samples sspliced = spliced.samples();
  QCOMPARE(spliced.sampleIndex(0), 0);
  QCOMPARE(spliced.sampleIndex(2), 2);

  QCOMPARE(sspliced2 == samplesT2, true);
  QCOMPARE(spliced.nSamples(), 3);

  QCOMPARE(spliced.gl(0, 0, 1, 1), 1.0);
  QCOMPARE(spliced.gl(0, 0, 1, 0), 0.0);
  QCOMPARE(spliced.gl(0, 0, 0, 1), 0.0);
  QCOMPARE(spliced.gl(0, 0, 0, 0), 0.0);
  QCOMPARE(spliced.gl(0, 1, 1, 1), 0.0);
  QCOMPARE(spliced.gl(0, 1, 1, 0), 0.0);
  QCOMPARE(spliced.gl(0, 1, 0, 1), 1.0);
  QCOMPARE(spliced.gl(0, 1, 0, 0), 0.0);
  QCOMPARE(spliced.gl(0, 2, 1, 1), 0.0);
  QCOMPARE(spliced.gl(0, 2, 1, 0), 1.0);
  QCOMPARE(spliced.gl(0, 2, 0, 1), 0.0);
  QCOMPARE(spliced.gl(0, 2, 0, 0), 0.0);
  QCOMPARE(spliced.gl(1, 1, 1, 1), 0.0);
  QCOMPARE(spliced.gl(1, 1, 1, 0), 0.0);
  QCOMPARE(spliced.gl(1, 1, 0, 1), 1.0);
  QCOMPARE(spliced.gl(1, 1, 0, 0), 0.0);
  QCOMPARE(spliced.gl(2, 0, 1, 1), 0.0);
  QCOMPARE(spliced.gl(2, 0, 1, 0), 1.0);
  QCOMPARE(spliced.gl(2, 0, 0, 1), 1.0);
  QCOMPARE(spliced.gl(2, 0, 0, 0), 1.0);
  QCOMPARE(spliced.gl(2, 1, 1, 1), 0.0);
  QCOMPARE(spliced.gl(2, 1, 1, 0), 1.0);
  QCOMPARE(spliced.gl(2, 1, 0, 1), 0.0);
  QCOMPARE(spliced.gl(2, 1, 0, 0), 0.0);
  QCOMPARE(spliced.gl(2, 2, 1, 1), 0.0);
  QCOMPARE(spliced.gl(2, 2, 1, 0), 1.0);
  QCOMPARE(spliced.gl(2, 2, 0, 1), 1.0);
  QCOMPARE(spliced.gl(2, 2, 0, 0), 0.0);
  QCOMPARE(spliced.isPhased(0, 0), true);
  QCOMPARE(spliced.isPhased(0, 1), true);
  QCOMPARE(spliced.isPhased(0, 2), true);
  QCOMPARE(spliced.isPhased(1, 0), true);
  QCOMPARE(spliced.isPhased(1, 1), true);
  QCOMPARE(spliced.isPhased(1, 2), true);
  QCOMPARE(spliced.isPhased(2, 0), false);
  QCOMPARE(spliced.isPhased(2, 1), true);
  QCOMPARE(spliced.isPhased(2, 2), false);

  SplicedGL revphased(phasedtarggl, true);

  QCOMPARE(revphased.allele1(1, 1), 0);
  QCOMPARE(revphased.allele2(1, 1), 1);
  QCOMPARE(revphased.allele2(0, 0), 1);
  QCOMPARE(revphased.allele1(2, 2), 1);
  QCOMPARE(revphased.allele(2, 4), 1);
  QCOMPARE(revphased.allele(0, 1), 1);
  QCOMPARE(revphased.nMarkers(), 3);
  QCOMPARE(revphased.nAlleles(1), 2);
  QCOMPARE(revphased.isRefData(), true);

  QCOMPARE(revphased.nHaps(), 6);
  QCOMPARE(revphased.nHapPairs(), 3);
  Samples srevphased2 = revphased.samples(2);
  Samples srevphased = revphased.samples();
  QCOMPARE(revphased.sampleIndex(0), 0);
  QCOMPARE(revphased.sampleIndex(2), 2);

  QCOMPARE(srevphased2 == samplesT, true);
  QCOMPARE(revphased.nSamples(), 3);

  Markers mksrph = revphased.markers();
  Marker m2rph = revphased.marker(2);
  QCOMPARE(mksrph == markss, false);
  QCOMPARE(m2rph == m2, false);
  QCOMPARE(m2rph == m3, false);
  QCOMPARE(revphased.marker(0) == m3, true);

  QCOMPARE(revphased.gl(2, 0, 1, 1), 1.0);
  QCOMPARE(revphased.gl(2, 0, 1, 0), 0.0);
  QCOMPARE(revphased.gl(2, 0, 0, 1), 0.0);
  QCOMPARE(revphased.gl(2, 0, 0, 0), 0.0);
  QCOMPARE(revphased.gl(2, 1, 1, 1), 0.0);
  QCOMPARE(revphased.gl(2, 1, 1, 0), 0.0);
  QCOMPARE(revphased.gl(2, 1, 0, 1), 1.0);
  QCOMPARE(revphased.gl(2, 1, 0, 0), 0.0);
  QCOMPARE(revphased.gl(2, 2, 1, 1), 0.0);
  QCOMPARE(revphased.gl(2, 2, 1, 0), 1.0);
  QCOMPARE(revphased.gl(2, 2, 0, 1), 0.0);
  QCOMPARE(revphased.gl(2, 2, 0, 0), 0.0);
  QCOMPARE(revphased.gl(1, 1, 1, 1), 0.0);
  QCOMPARE(revphased.gl(1, 1, 1, 0), 0.0);
  QCOMPARE(revphased.gl(1, 1, 0, 1), 1.0);
  QCOMPARE(revphased.gl(1, 1, 0, 0), 0.0);
  QCOMPARE(revphased.gl(0, 0, 1, 1), 0.0);
  QCOMPARE(revphased.gl(0, 0, 1, 0), 0.0);
  QCOMPARE(revphased.gl(0, 0, 0, 1), 1.0);
  QCOMPARE(revphased.gl(0, 0, 0, 0), 0.0);
  QCOMPARE(revphased.isPhased(2, 0), true);
  QCOMPARE(revphased.isPhased(2, 1), true);
  QCOMPARE(revphased.isPhased(2, 2), true);
  QCOMPARE(revphased.isPhased(1, 0), true);
  QCOMPARE(revphased.isPhased(1, 1), true);
  QCOMPARE(revphased.isPhased(1, 2), true);
  QCOMPARE(revphased.isPhased(0, 0), true);
  QCOMPARE(revphased.isPhased(0, 1), true);
  QCOMPARE(revphased.isPhased(0, 2), true);

  SplicedGL revspliced(spliced, true);

  QCOMPARE(revspliced.allele1(1, 1), 0);
  QCOMPARE(revspliced.allele2(1, 1), 1);
  QCOMPARE(revspliced.allele2(0, 0), -1);
  QCOMPARE(revspliced.allele1(2, 2), 1);
  QCOMPARE(revspliced.allele(2, 4), 1);
  QCOMPARE(revspliced.allele(0, 1), -1);
  QCOMPARE(revspliced.nMarkers(), 3);
  QCOMPARE(revspliced.nAlleles(1), 2);
  QCOMPARE(revspliced.isRefData(), false);

  QCOMPARE(revspliced.nHaps(), 6);
  QCOMPARE(revspliced.nHapPairs(), 3);
  Samples srevspliced2 = revspliced.samples(2);
  Samples srevspliced = revspliced.samples();
  QCOMPARE(revspliced.sampleIndex(0), 0);
  QCOMPARE(revspliced.sampleIndex(2), 2);

  QCOMPARE(srevspliced2 == samplesT2, true);
  QCOMPARE(revspliced.nSamples(), 3);

  Markers mksrsp = revspliced.markers();
  Marker m2rsp = revspliced.marker(2);
  QCOMPARE(mksrsp == markss, false);
  QCOMPARE(m2rsp == m2, false);
  QCOMPARE(m2rsp == m3, false);
  QCOMPARE(revspliced.marker(0) == m3, true);

  QCOMPARE(revspliced.gl(2, 0, 1, 1), 1.0);
  QCOMPARE(revspliced.gl(2, 0, 1, 0), 0.0);
  QCOMPARE(revspliced.gl(2, 0, 0, 1), 0.0);
  QCOMPARE(revspliced.gl(2, 0, 0, 0), 0.0);
  QCOMPARE(revspliced.gl(2, 1, 1, 1), 0.0);
  QCOMPARE(revspliced.gl(2, 1, 1, 0), 0.0);
  QCOMPARE(revspliced.gl(2, 1, 0, 1), 1.0);
  QCOMPARE(revspliced.gl(2, 1, 0, 0), 0.0);
  QCOMPARE(revspliced.gl(2, 2, 1, 1), 0.0);
  QCOMPARE(revspliced.gl(2, 2, 1, 0), 1.0);
  QCOMPARE(revspliced.gl(2, 2, 0, 1), 0.0);
  QCOMPARE(revspliced.gl(2, 2, 0, 0), 0.0);
  QCOMPARE(revspliced.gl(1, 1, 1, 1), 0.0);
  QCOMPARE(revspliced.gl(1, 1, 1, 0), 0.0);
  QCOMPARE(revspliced.gl(1, 1, 0, 1), 1.0);
  QCOMPARE(revspliced.gl(1, 1, 0, 0), 0.0);
  QCOMPARE(revspliced.gl(0, 0, 1, 1), 0.0);
  QCOMPARE(revspliced.gl(0, 0, 1, 0), 1.0);
  QCOMPARE(revspliced.gl(0, 0, 0, 1), 1.0);
  QCOMPARE(revspliced.gl(0, 0, 0, 0), 1.0);
  QCOMPARE(revspliced.gl(0, 1, 1, 1), 0.0);
  QCOMPARE(revspliced.gl(0, 1, 1, 0), 1.0);
  QCOMPARE(revspliced.gl(0, 1, 0, 1), 0.0);
  QCOMPARE(revspliced.gl(0, 1, 0, 0), 0.0);
  QCOMPARE(revspliced.gl(0, 2, 1, 1), 0.0);
  QCOMPARE(revspliced.gl(0, 2, 1, 0), 1.0);
  QCOMPARE(revspliced.gl(0, 2, 0, 1), 1.0);
  QCOMPARE(revspliced.gl(0, 2, 0, 0), 0.0);
  QCOMPARE(revspliced.isPhased(2, 0), true);
  QCOMPARE(revspliced.isPhased(2, 1), true);
  QCOMPARE(revspliced.isPhased(2, 2), true);
  QCOMPARE(revspliced.isPhased(1, 0), true);
  QCOMPARE(revspliced.isPhased(1, 1), true);
  QCOMPARE(revspliced.isPhased(1, 2), true);
  QCOMPARE(revspliced.isPhased(0, 0), false);
  QCOMPARE(revspliced.isPhased(0, 1), true);
  QCOMPARE(revspliced.isPhased(0, 2), false);

  FuzzyGL fzgl(revspliced, .25, false);
  QCOMPARE(fzgl.gl(0, 0, 1, 1), 0.0);
  QCOMPARE(fzgl.gl(0, 0, 1, 0), 1.0);
  QCOMPARE(fzgl.gl(0, 0, 0, 1), 1.0);
  QCOMPARE(fzgl.gl(0, 0, 0, 0), 1.0);
  QCOMPARE(fzgl.gl(2, 0, 1, 1), 0.5625);
  QCOMPARE(fzgl.gl(2, 0, 1, 0), 0.1875);
  QCOMPARE(fzgl.gl(2, 0, 0, 1), 0.1875);
  QCOMPARE(fzgl.gl(2, 0, 0, 0), 0.0625);
  QCOMPARE(fzgl.gl(0, 1, 1, 1), 0.1875);
  QCOMPARE(fzgl.gl(0, 1, 1, 0), 0.5625);
  QCOMPARE(fzgl.gl(0, 1, 0, 1), 0.0625);
  QCOMPARE(fzgl.gl(0, 1, 0, 0), 0.1875);
  QCOMPARE(fzgl.gl(0, 2, 1, 1), 0.375);
  QCOMPARE(fzgl.gl(0, 2, 1, 0), 0.625);
  QCOMPARE(fzgl.gl(0, 2, 0, 1), 0.625);
  QCOMPARE(fzgl.gl(0, 2, 0, 0), 0.375);

  GLUser gluser;

  gluser.myGL = fzgl;
  QCOMPARE(gluser.myGL.gl(0, 0, 1, 1), 0.0);
  QCOMPARE(gluser.myGL.gl(0, 0, 1, 0), 1.0);
  QCOMPARE(gluser.myGL.gl(0, 0, 0, 1), 1.0);
  QCOMPARE(gluser.myGL.gl(0, 0, 0, 0), 1.0);
  QCOMPARE(gluser.myGL.gl(2, 0, 1, 1), 0.5625);
  QCOMPARE(gluser.myGL.gl(2, 0, 1, 0), 0.1875);
  QCOMPARE(gluser.myGL.gl(2, 0, 0, 1), 0.1875);
  QCOMPARE(gluser.myGL.gl(2, 0, 0, 0), 0.0625);
  QCOMPARE(gluser.myGL.gl(0, 1, 1, 1), 0.1875);
  QCOMPARE(gluser.myGL.gl(0, 1, 1, 0), 0.5625);
  QCOMPARE(gluser.myGL.gl(0, 1, 0, 1), 0.0625);
  QCOMPARE(gluser.myGL.gl(0, 1, 0, 0), 0.1875);
  QCOMPARE(gluser.myGL.gl(0, 2, 1, 1), 0.375);
  QCOMPARE(gluser.myGL.gl(0, 2, 1, 0), 0.625);
  QCOMPARE(gluser.myGL.gl(0, 2, 0, 1), 0.625);
  QCOMPARE(gluser.myGL.gl(0, 2, 0, 0), 0.375);
}

void TestImputeDataStructures::testRefTargetData3x3()
{
  clearStaticTestLists();

  RefDataReader rr;
  TargDataReaderTest3x3 tr(true);  // "true" setting gives reference-ready data.

  TargetData td;

  testWindowDriver(td, tr, rr, 4);

  QCOMPARE(overlapAmountsTestList.length(), 2);
  // QCOMPARE(targPairsTestList.length(), 2);
}

void TestImputeDataStructures::testTargetData3x3()
{
  clearStaticTestLists();

  RefDataReader rr;
  TargDataReaderTest3x3 tr(false);  // "false" setting gives unphased data.

  TargetData td;

  testWindowDriver(td, tr, rr, 4);

  QCOMPARE(overlapAmountsTestList.length(), 1);
  // QCOMPARE(targPairsTestList.length(), 1);
}

void TestImputeDataStructures::testAllData4x4and3x3()
{
  /*
  clearStaticTestLists();

  RefDataReaderTest4x4 rr;
  TargDataReaderTest3x3 tr(false); // "false" setting gives unphased data.

  AllData ad;

  testWindowDriver(ad, tr, rr, 4);

  QCOMPARE(overlapAmountsTestList.length(), 1);
  QCOMPARE(targPairsTestList.length(), 1);
  */
}

QTEST_MAIN(TestImputeDataStructures)
#include "test_impute_datastructures.moc"
