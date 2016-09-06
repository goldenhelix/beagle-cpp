/* Copyright 2016 Golden Helix, Inc. */
#include "impute/haplotypepair.h"
#include "impute/samples.h"
#include "impute/markers.h"
#include "impute/vcfEmission.h"

#include <QObject>
#include <QtTest/QtTest>

class TestImputeDataStructures : public QObject
{
  Q_OBJECT;
private slots:

  void testHaplotypePairs();
  void testSamples();
  void testMarkers();
  void testVcfEmissions();
};

void TestImputeDataStructures::testHaplotypePairs()
{
  // HaplotypePair pair;
  // pair.setData(qMakePair<QByteArray, QByteArray>("ABC", "DEF"));

  // QCOMPARE(pair.a(), QString("ABC"));
  // QCOMPARE(pair.b(), QString("DEF"));
}

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
  Markers marksrev = marksfwd.reverse();

  QCOMPARE(marksfwd.nMarkers(), 6);
  QCOMPARE(marksrev.nMarkers(), 6);
  QCOMPARE(marksfwd.marker(2).pos(), 5432);
  QCOMPARE(marksrev.marker(2).pos(), 23465);
  QCOMPARE(marksfwd.markers()[1].id().constData(), "RS72351");

  QCOMPARE(marksfwd.restrict(3, 5).markers() == lm, true);
  QList<Marker> restrictedList = marksfwd.restrict(3, 5).markers();
  QCOMPARE(restrictedList.length(), 2);
  QCOMPARE(restrictedList[0] == lm[0], true);
  QCOMPARE(restrictedList[1] == lm[1], true);
  QCOMPARE(restrictedList == lm, true);

  QCOMPARE(marksfwd.sumAlleles(2), 4);
  QCOMPARE(marksfwd.sumAlleles(), 12);
  QCOMPARE(marksfwd.sumGenotypes(2), 6);
  QCOMPARE(marksfwd.sumGenotypes(), 18);
  QCOMPARE(marksfwd.sumHaplotypeBits(2), 2);
  QCOMPARE(marksfwd.sumHaplotypeBits(), 6);

  Markers markers2(marksfwd);
  QCOMPARE(markers2.marker(3) == marksfwd.marker(3), true);

  Markers markers3 = markers2;
  QCOMPARE(markers3.marker(2) == marksfwd.marker(2), true);
}

void TestImputeDataStructures::testVcfEmissions()
{
  QList<BitSetRefGT> _refEmissions;
  QList<BitSetGT> _targetEmissions;

  QCOMPARE(SampleNames::getIndexIfIndexed("SAMP073"), 3);       // This name should exist already.
  QCOMPARE(SampleNames::getIndexIfIndexed("SAMP007"), -1);      // This name should not.

  // Construct "samplesR" and set its data before there are any other
  // references to the object.

  Samples samplesR;
  samplesR.setSamp(SampleNames::getIndex("SAMP001"));   // The first three names and global sample
  samplesR.setSamp(SampleNames::getIndex("SAMP003"));   // indexes already exist (globally).
  samplesR.setSamp(SampleNames::getIndex("SAMP005"));
  samplesR.setSamp(SampleNames::getIndex("SAMP007"));   // This name is new.

  QCOMPARE(SampleNames::getIndexIfIndexed("SAMP073"), 3);
  QCOMPARE(SampleNames::getIndexIfIndexed("SAMP007"), 5);

  // Build some VcfEmissions objects as "reference data" (containing
  // phased genotypes without missing data). Be careful to load the
  // data before doing any "copy" or "assignment" operations (which
  // operations are actually reference copy operations).

  QCOMPARE(ChromeIds::getIndexIfIndexed("X"),  2);      // This chromosome name should exist already.
  QCOMPARE(ChromeIds::getIndexIfIndexed("2"), -1);      // This one should not.

  BitSetRefGT r1(samplesR);

  r1.setIdInfo(ChromeIds::getIndex("1"), 12345, "RS12345");
  r1.setAllele("A");
  r1.setAllele("C");

  QList<int> r11; r11.append(1); r11.append(0); r11.append(0); r11.append(1);
  QList<int> r12; r12.append(1); r12.append(1); r12.append(0); r12.append(1);
  QList<bool> arePhased; arePhased.append(true); arePhased.append(true); arePhased.append(true); arePhased.append(true);
  r1.storeAlleles(r11, r12, arePhased);

  _refEmissions.append(r1);

  BitSetRefGT r2(samplesR);

  r2.setIdInfo(ChromeIds::getIndex("17"), 22345, "RS22345");
  r2.setAllele("G");
  r2.setAllele("T");

  QList<int> r21; r21.append(0); r21.append(0); r21.append(1); r21.append(0);
  QList<int> r22; r22.append(1); r22.append(1); r22.append(1); r22.append(1);
  r2.storeAlleles(r21, r22, arePhased);

  _refEmissions.append(r2);

  BitSetRefGT r2a(samplesR);  // This marker will be inbetween two markers that overlap the target data.

  r2a.setIdInfo(ChromeIds::getIndex("17"), 22678, "RS22678");
  r2a.setAllele("C");
  r2a.setAllele("G");

  QList<int> r21a; r21a.append(0); r21a.append(1); r21a.append(1); r21a.append(0);
  QList<int> r22a; r22a.append(1); r22a.append(0); r22a.append(1); r22a.append(1);
  r2a.storeAlleles(r21a, r22a, arePhased);

  _refEmissions.append(r2a);

  BitSetRefGT r3(samplesR);

  r3.setIdInfo(ChromeIds::getIndex("17"), 33345, "RS33345");
  r3.setAllele("A");
  r3.setAllele("T");

  QList<int> r31; r31.append(0); r31.append(1); r31.append(0); r31.append(0);
  QList<int> r32; r32.append(1); r32.append(0); r32.append(1); r32.append(1);
  r3.storeAlleles(r31, r32, arePhased);

  _refEmissions.append(r3);


  // Construct "samplesT" and set its data before there are any other
  // references to the object.

  Samples samplesT;
  samplesT.setSamp(SampleNames::getIndex("SAMP073"));
  samplesT.setSamp(SampleNames::getIndex("SAMP087"));
  samplesT.setSamp(SampleNames::getIndex("SAMP095"));

  // Build some VcfEmissions objects as "target data" (containing
  // normal unphased genotypes). Be careful to load the data before
  // doing any "copy" or "assignment" operations (which are actually
  // reference copy operations).

  BitSetGT t1(samplesT);

  t1.setIdInfo(ChromeIds::getIndex("1"), 12345, "RS12345");
  t1.setAllele("A");
  t1.setAllele("C");

  QList<int> t11; t11.append(1); t11.append(0); t11.append(-1);
  QList<int> t12; t12.append(1); t12.append(1); t12.append(0);
  arePhased[0] = false; arePhased[2] = false; arePhased.removeLast();   // Leave arePhased[1] at "true" on purpose.
  t1.storeAlleles(t11, t12, arePhased);

  _targetEmissions.append(t1);

  BitSetGT t2(samplesT);

  t2.setIdInfo(ChromeIds::getIndex("17"), 22345, "RS22345");
  t2.setAllele("G");
  t2.setAllele("T");

  QList<int> t21; t21.append(0); t21.append(0); t21.append(1);
  QList<int> t22; t22.append(1); t22.append(1); t22.append(1);
  arePhased[0] = true; arePhased[2] = true;
  t2.storeAlleles(t21, t22, arePhased);

  _targetEmissions.append(t2);

  BitSetGT t3(samplesT);

  t3.setIdInfo(ChromeIds::getIndex("17"), 33345, "RS33345");
  t3.setAllele("A");
  t3.setAllele("T");

  QList<int> t31; t31.append(0); t31.append(1); t31.append(1);
  QList<int> t32; t32.append(1); t32.append(0); t32.append(0);
  arePhased[0] = false; arePhased[2] = false;
  t3.storeAlleles(t31, t32, arePhased);

  _targetEmissions.append(t3);

  QCOMPARE(_refEmissions.length(), 4);
  QCOMPARE(_targetEmissions.length(), 3);

  QCOMPARE(_refEmissions[3].isPhased(2), true);
  QCOMPARE(_targetEmissions[1].isPhased(1), true);
  QCOMPARE(_targetEmissions[0].isRefData(), false);
  QCOMPARE(t2.isRefData(), true);
  QCOMPARE(_targetEmissions[1].isRefData(), true);
  QCOMPARE(_refEmissions[2].isRefData(), true);
  QCOMPARE(_targetEmissions[0].allele1(2), -1);

  _refEmissions.clear();
  _targetEmissions.clear();

  QCOMPARE(_refEmissions.length(), 0);
  QCOMPARE(_targetEmissions.length(), 0);
}

QTEST_MAIN(TestImputeDataStructures)
#include "test_impute_datastructures.moc"
