/* Copyright 2016 Golden Helix, Inc. */
#include "impute/haplotypepair.h"
#include "impute/samples.h"

#include <QObject>
#include <QtTest/QtTest>

class TestImputeDataStructures : public QObject
{
  Q_OBJECT;
private slots:

  void testHaplotypePairs();
  void testSamples();
  void testMarker();
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

void TestImputeDataStructures::testMarker()
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

  QList<CString> als;
  als.append("A");
  als.append("C");

  // To use objects of a class (such as Marker) in a QList, etc., the
  // class must have:
  //
  // * a default constructor,
  // * a copy constructor, and
  // * an assignment operator.

  Marker m();
  m.setIdInfo(c1, 7632, "RS72351");
  m.setAllele("A");
  m.setAllele("C");

  QCOMPARE(m.chrom().constData(), "1");
  QCOMPARE(m.chromIndex(), 0);
  QCOMPARE(m.pos(), 7632);
  QCOMPARE(m.nAlleles(), 2);
  QCOMPARE(m.allele(0).constData(), "A");
  QCOMPARE(m.allele(1).constData(), "C");

  Marker m2();
  m2.setIdInfo(c1, 7632, "RS72351");
  m2.setAllele("A");
  m2.setAllele("C");

  QCOMPARE(m == m2, true);
  QCOMPARE(m, m2);

  Marker m3();
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
  QCOMPARE(m4, m2);

  Marker m5(m2);
  QCOMPARE(m5.pos(), 7632);
  QCOMPARE(m5, m2);
}

QTEST_MAIN(TestImputeDataStructures)
#include "test_impute_datastructures.moc"
