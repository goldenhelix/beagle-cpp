/* Copyright 2016 Golden Helix, Inc. */
#include "impute/haplotypepair.h"
#include "impute/samples.h"
#include "impute/markers.h"

#include <QObject>
#include <QtTest/QtTest>

class TestImputeDataStructures : public QObject
{
  Q_OBJECT;
private slots:

  void testHaplotypePairs();
  void testSamples();
  void testMarkers();
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
  SampleNames names;
  int n1 = names.getIndex("SAMP001");
  int n2 = names.getIndex("SAMP003");
  int n3 = names.getIndex("SAMP005");
  int n4 = names.getIndexIfIndexed("SAMP003");
  int n5 = names.getIndexIfIndexed("SAMP072");

  QCOMPARE(names.nNames(), 3);
  QCOMPARE(n1, 0);
  QCOMPARE(n2, 1);
  QCOMPARE(n3, 2);
  QCOMPARE(n4, 1);
  QCOMPARE(n5, -1);
  QCOMPARE(names.name(0).constData(), "SAMP001");
  QCOMPARE(names.name(1).constData(), "SAMP003");
  QCOMPARE(names.name(2).constData(), "SAMP005");

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
  samples2.setSamp(names.getIndex("SAMP073"));
  samples2.setSamp(names.getIndex("SAMP085"));

  QCOMPARE(names.name(3).constData(), "SAMP073");
  QCOMPARE(names.name(4).constData(), "SAMP085");

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
  ChromeIds chIds;
  int c1 = chIds.getIndex("1");
  int c17 = chIds.getIndex("17");
  int cx = chIds.getIndex("X");
  int c1a = chIds.getIndexIfIndexed("1");
  int c2 = chIds.getIndexIfIndexed("2");
  QCOMPARE(chIds.nIds(), 3);
  QCOMPARE(chIds.chromeId(0).constData(), "1");
  QCOMPARE(chIds.chromeId(1).constData(), "17");
  QCOMPARE(chIds.chromeId(2).constData(), "X");
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

QTEST_MAIN(TestImputeDataStructures)
#include "test_impute_datastructures.moc"
