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
}

QTEST_MAIN(TestImputeDataStructures)
#include "test_impute_datastructures.moc"
