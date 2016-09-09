/* Copyright 2016 Golden Helix, Inc. */

#include "impute/samples.h"
#include "impute/markers.h"
#include "impute/vcfemission.h"
#include "impute/haplotypepair.h"

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
  QCOMPARE(samples.index(0), 1);
  QCOMPARE(samples.index(1), 2);
  QCOMPARE(samples.findIndex(1), 0);
  QCOMPARE(samples.findIndex(2), 1);
  QCOMPARE(samples.findIndex(0), -1);
  QCOMPARE(samples.name(0).constData(), "SAMP003");
  QCOMPARE(samples.name(1).constData(), "SAMP005");

  Samples samples2;
  samples2.setSamp(SampleNames::getIndex("SAMP073"));
  samples2.setSamp(SampleNames::getIndex("SAMP085"));

  QCOMPARE(SampleNames::name(3).constData(), "SAMP073");
  QCOMPARE(SampleNames::name(4).constData(), "SAMP085");

  QCOMPARE(samples2.nSamples(), 2);
  QCOMPARE(samples2.index(0), 3);
  QCOMPARE(samples2.index(1), 4);
  QCOMPARE(samples2.findIndex(4), 1);
  QCOMPARE(samples2.findIndex(3), 0);
  QCOMPARE(samples2.findIndex(2), -1);
  QCOMPARE(samples2.name(0).constData(), "SAMP073");
  QCOMPARE(samples2.name(1).constData(), "SAMP085");

  Samples samples3 = samples2;

  QCOMPARE(samples2.nSamples(), 2);
  QCOMPARE(samples2.index(0), 3);
  QCOMPARE(samples2.index(1), 4);
  QCOMPARE(samples2.findIndex(4), 1);
  QCOMPARE(samples2.findIndex(3), 0);
  QCOMPARE(samples2.findIndex(2), -1);
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
  m9.setIdInfo(c1, 7632, "RS72355");   // Should not make it any "different" from marker m.
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

  QCOMPARE(ChromeIds::getIndexIfIndexed("X"), 2);       // This chromosome name should exist already.
  QCOMPARE(ChromeIds::getIndexIfIndexed("2"), -1);      // This one should not.

  loadTestDataForTargetData(targetEmissions);

  QCOMPARE(ChromeIds::getIndexIfIndexed("X"),  2);      // This chromosome name should exist already.
  QCOMPARE(ChromeIds::getIndexIfIndexed("2"), -1);      // This one still should not.

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
  QList<int> als11;
  QList<int> als21;
  QList<int> als13;
  QList<int> als23;

  for (int mnum=0; mnum < nummarks; mnum++)
  {
    biglist.append(refEmissions[mnum].marker());
    als11.append(refEmissions[mnum].allele1(1));
    als21.append(refEmissions[mnum].allele2(1));
    als13.append(refEmissions[mnum].allele1(3));
    als23.append(refEmissions[mnum].allele2(3));
  }

  Markers marks(biglist);
  HapPair pair1(marks, refEmissions[0].samples(), 1, als11, als21);
  HapPair pair3(marks, refEmissions[0].samples(), 3, als13, als23);

  HapPair pair1copy(pair1, false); // Copy but don't reverse.
  HapPair pair3rev(pair3, true);   // Copy and reverse this one.

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
