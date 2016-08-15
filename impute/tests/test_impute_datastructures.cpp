/* Copyright 2016 Golden Helix, Inc. */
#include "impute/haplotypepair.h"

#include <QObject>
#include <QtTest/QtTest>


class TestImputeDataStructures : public QObject
{
  Q_OBJECT;
private slots:

  void testHaplotypePairs();
};

void TestImputeDataStructures::testHaplotypePairs()
{
  HaplotypePair pair;
  pair.setData(qMakePair<QByteArray, QByteArray>("ABC", "DEF"));

  QCOMPARE(pair.a(), QString("ABC"));
  QCOMPARE(pair.b(), QString("DEF"));
}

QTEST_MAIN(TestImputeDataStructures)
#include "test_impute_datastructures.moc"
