/* Copyright 2016 Golden Helix, Inc. */

#include "impute/haplotypepair.h"
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
  QList<BitSetRefGT> refEmissions;
  QList<BitSetGT> targetEmissions;

  loadTestDataForRefData(refEmissions);

  QCOMPARE(SampleNames::getIndexIfIndexed("SAMP073"), 3);
  QCOMPARE(SampleNames::getIndexIfIndexed("SAMP007"), 5);

  QCOMPARE(ChromeIds::getIndexIfIndexed("X"), 2);   // This chromosome name should exist already.
  QCOMPARE(ChromeIds::getIndexIfIndexed("2"), -1);  // This one should not.

  QCOMPARE(refEmissions[0].allele2(1), 1);
  QCOMPARE(refEmissions[1].allele1(0), 0);
  QCOMPARE(refEmissions[2].allele1(2), 1);
  QCOMPARE(refEmissions[3].allele2(3), 1);

  loadTestDataForTargetData(targetEmissions);

  QCOMPARE(ChromeIds::getIndexIfIndexed("X"), 2);   // This chromosome name should exist already.
  QCOMPARE(ChromeIds::getIndexIfIndexed("2"), -1);  // This one still should not.

  QCOMPARE(refEmissions.length(), 4);
  QCOMPARE(targetEmissions.length(), 3);

  QCOMPARE(targetEmissions[1].isPhased(1), true);
  QCOMPARE(targetEmissions[0].isRefData(), false);
  QCOMPARE(targetEmissions[1].isRefData(), true);
  QCOMPARE(targetEmissions[0].allele1(2), -1);
  QCOMPARE(targetEmissions[0].allele2(1), 1);
  QCOMPARE(targetEmissions[1].allele1(0), 0);
  QCOMPARE(targetEmissions[2].allele1(2), 1);

  refEmissions.clear();
  targetEmissions.clear();

  QCOMPARE(refEmissions.length(), 0);
  QCOMPARE(targetEmissions.length(), 0);
}

void TestImputeDataStructures::testHaplotypePairs()
{
  QList<BitSetRefGT> refEmissions;
  loadTestDataForRefData(refEmissions);

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
  HapPair pair1(marks, refEmissions[0].samples(), 1, als11, als21);
  HapPair pair3(marks, refEmissions[0].samples(), 3, als13, als23);

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

  QCOMPARE(hps.nHaps(), 6);
  QCOMPARE(hps.nHapPairs(), 3);
  Samples shps1 = hps.samples(1);
  QCOMPARE(hps.sampleIndex(0), 1);
  QCOMPARE(hps.sampleIndex(1), 3);
  QCOMPARE(hps.sampleIndex(2), 1);

  QCOMPARE(mks == marks, true);
  QCOMPARE(m2 == refEmissions[2].marker(), true);
  QCOMPARE(shps1 == refEmissions[1].samples(), true);

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

  HapPair pair0(marks, refEmissions[0].samples(), 0, als10, als20);
  HapPair pair2(marks, refEmissions[0].samples(), 2, als12, als22);

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
  samplesCompleteX.setSamp(SampleNames::getIndex("SAMP001"));
  samplesCompleteX.setSamp(
      SampleNames::getIndex("SAMP005"));  // Notice the order is switched, here....
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

  QList<BitSetGT> phasedTargetEmissions;
  loadTestDataForTargetData(phasedTargetEmissions, 1, true);
  Samples samplesT = phasedTargetEmissions[1].samples();

  /*
  Samples samplesT;
  samplesT.setSamp(SampleNames::getIndex("SAMP073"));
  samplesT.setSamp(SampleNames::getIndex("SAMP087"));
  samplesT.setSamp(SampleNames::getIndex("SAMP095"));
  */

  SplicedGL phasedtarggl(samplesT, phasedTargetEmissions);
  GLSampleHapPairs glshp(phasedtarggl);

  QList<BitSetGT> unphasedTargetEmissions;
  loadTestDataForTargetData(unphasedTargetEmissions);
  Samples samplesU = unphasedTargetEmissions[2].samples();
  SplicedGL unphasedtarggl(samplesU, unphasedTargetEmissions);
  // This one should not work, and doesn't. GLSampleHapPairs crash(unphasedtarggl);

  // FOR THE NEXT ONE, NEED TO SET UP PHASED HAPS FOR TWO OF THE THREE MARKERS....
  // SplicedGL spliced(SampleHapPairs &haps, GLSampleHapPairs &otherGL);

  SplicedGL revphased(phasedtarggl, true);
}

/*
void TestImputeDataStructures::testDataDrivers()
{
  NullDataReader nr;
  TargReaderTest tr;
  RefReaderTest rr;

  QList<BitSetRefGT> refEmissions;
  QList<BitSetGT> targetEmissions;

  rr.initialize(refEmissions);
  tr.initialize(targetEmissions);

  VcfWindow refWind;
  RestrictedVcfWindow targetWind;
  AllData ad(refWind, targetWind);

  int overlap = 0;
  ad.advanceWindow(overlap, tr, rr);
}
*/

QTEST_MAIN(TestImputeDataStructures)
#include "test_impute_datastructures.moc"
