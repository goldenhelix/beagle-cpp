#include "impute/consensusphaser.h"

#include <QVector>

static void checkMarkers(const QList<HapPair> &hapList)
{
  Markers markers = hapList[0].markers();
  Samples samples = hapList[0].samples();

  for (int j = 1; j < hapList.size(); ++j) {
    Q_ASSERT_X(hapList[j].markers() == markers,
               "checkMarkers (consensusphaser.cpp)",
               "inconsistent markers");

    Q_ASSERT_X(hapList[j].samples() == samples,
               "checkMarkers (consensusphaser.cpp)",
               "inconsistent samples");
  }
}

static bool hapsComparator(const HapPair &hp1, const HapPair &hp2)
{
  return (hp1.sampleIndex() < hp2.sampleIndex());
}

static ConsensusPhaser::Phase flip(const ConsensusPhaser::Phase &phase)
{
  if (phase == ConsensusPhaser::IDENTICAL)
    return ConsensusPhaser::OPPOSITE;
  else if (phase == ConsensusPhaser::OPPOSITE)
    return ConsensusPhaser::IDENTICAL;
  else
  {
    Q_ASSERT_X(false, "flip (consensusphaser.cpp)", "bad phase");
    return ConsensusPhaser::INCONSISTENT;
  }
}

static void storePhase(const QList<HapPair> &hapList, int marker, int a1, int a2,
                       QVector<int> &phaseArray)
{
  Q_ASSERT_X(phaseArray.length() == hapList.size(),
             "storePhase (consensusphaser.cpp)",
             "sizes of hapList and phaseArray do not match");

  for (int j = 0; j < phaseArray.length(); ++j) {
    int b1 = hapList[j].allele1(marker);
    int b2 = hapList[j].allele2(marker);
    if ((a1 == b1 && a2 == b2) || (a1 == b2 && a2 == b1)) {
      phaseArray[j] = (b1 < b2) ? ConsensusPhaser::IDENTICAL : ConsensusPhaser::OPPOSITE;
    } else {
      phaseArray[j] = ConsensusPhaser::INCONSISTENT;
    }
  }
}

static int relPhase(QVector<int> ph1, QVector<int> ph2 /*, Random rand */)
{
  Q_ASSERT_X(
      ph1.length() == ph2.length(), "relPhase (consensusphaser.cpp)", "ph1.length != ph2.length");

  int identCnt = 0;
  int oppCnt = 0;

  for (int j = 0; j < ph1.length(); ++j) {
    if (ph1[j] == ConsensusPhaser::IDENTICAL) {
      switch (ph2[j]) {
        case ConsensusPhaser::IDENTICAL:
          ++identCnt;
          break;
        case ConsensusPhaser::OPPOSITE:
          ++oppCnt;
          break;
      }
    } else if (ph1[j] == ConsensusPhaser::OPPOSITE) {
      switch (ph2[j]) {
        case ConsensusPhaser::IDENTICAL:
          ++oppCnt;
          break;
        case ConsensusPhaser::OPPOSITE:
          ++identCnt;
          break;
      }
    }
  }

  if (identCnt > oppCnt)
    return ConsensusPhaser::IDENTICAL;
  else if (oppCnt > identCnt)
    return ConsensusPhaser::OPPOSITE;
  else {
    /// return rand.nextBoolean() ? ConsensusPhaser::IDENTICAL : ConsensusPhaser::OPPOSITE;
    return ConsensusPhaser::IDENTICAL;
  }
}

/**
 * Returns the VCF genotype index for the specified pair of alleles.
 * @param a1 the first allele
 * @param a2 the second allele
 */
static int gtIndex(int a1, int a2)
{
  Q_ASSERT_X(a1 >= 0, "gtIndex (consensusphaser.cpp)", "a1 < 0");

  Q_ASSERT_X(a2 >= 0, "gtIndex (consensusphaser.cpp)", "a2 < 0");

  if (a1 < a2)
    return (a2 * (a2 + 1)) / 2 + a1;
  else
    return (a1 * (a1 + 1)) / 2 + a2;
}

static QVector<int> gtCounts(const QList<HapPair> &hapList, const Markers &markers, int marker)
{
  int nGt = markers.marker(marker).nGenotypes();
  QVector<int> gtc(nGt, 0);
  for (int j = 0, n = hapList.size(); j < n; ++j) {
    const HapPair &hp = hapList[j];
    int a1 = hp.allele1(marker);
    int a2 = hp.allele2(marker);
    int gt = gtIndex(a1, a2);
    ++gtc[gt];
  }

  return gtc;
}

static int consensusGT(const QList<HapPair> &hapList, const Markers &markers,
                       int marker /*, Random random */)
{
  QVector<int> gtc = gtCounts(hapList, markers, marker);

  int gtclen = gtc.length();
  /// int start = random.nextInt(gtclen);
  int start = gtclen / 2;
  int bestGt = start;
  for (int j = 1; j < gtclen; ++j) {
    int gt = start + j;

    if (gt >= gtclen)
      gt -= gtclen;

    if (gtc[gt] > gtc[bestGt])
      bestGt = gt;
  }

  return bestGt;
}

static int hapPairWithConsensusGT(const QList<HapPair> &hapList, const Markers &markers,
                                  int marker /*, Random random */)
{
  int bestGT = consensusGT(hapList, markers, marker /*, random */);

  for (int j = 0, n = hapList.size(); j < n; ++j) {
    const HapPair &hp = hapList[j];
    int a1 = hp.allele1(marker);
    int a2 = hp.allele2(marker);
    if (gtIndex(a1, a2) == bestGT)
      return j;
  }

  Q_ASSERT_X(false, "hapPairWithConsensusGT (consensusphaser.cpp)", "no sample with consensus GT");
  return 0;
}

static HapPair consensus(const QList<HapPair> &hapList /*, Random rand */)
{
  const HapPair &firstHP = hapList[0];
  int sampleIndex = firstHP.sampleIndex();
  const Samples &samples = firstHP.samples();
  const Markers &markers = firstHP.markers();
  int nMarkers = markers.nMarkers();
  ConsensusPhaser::Phase lastConsensus = ConsensusPhaser::UNKNOWN;
  QVector<int> lastPhase(hapList.size());
  QVector<int> currentPhase(hapList.size());
  QList<int> alleles1;
  QList<int> alleles2;

  for (int m = 0; m < nMarkers; ++m) {
    int hp = hapPairWithConsensusGT(hapList, markers, m /*, rand */);
    // retrieve actual allele order to match input phased data
    int a1 = hapList[hp].allele1(m);
    int a2 = hapList[hp].allele2(m);
    if (a1 != a2) {
      storePhase(hapList, m, a1, a2, currentPhase);
      if (lastConsensus != ConsensusPhaser::UNKNOWN)
      {
        ConsensusPhaser::Phase thisConsensus;
		if (relPhase(lastPhase, currentPhase /*, rand */) == ConsensusPhaser::IDENTICAL)
          thisConsensus = lastConsensus;
        else
          thisConsensus = flip(lastConsensus);

        if ((thisConsensus == ConsensusPhaser::IDENTICAL && a1 > a2) ||
            (thisConsensus == ConsensusPhaser::OPPOSITE && a1 < a2)) {
          int tmpa = a1;
          a1 = a2;
          a2 = tmpa;
        }
      }

      lastConsensus = a1 < a2 ? ConsensusPhaser::IDENTICAL : ConsensusPhaser::OPPOSITE;
      lastPhase = currentPhase;
    }

    alleles1.append(a1);
    alleles2.append(a2);
  }

  return HapPair(markers, samples, sampleIndex, alleles1, alleles2);
}

/**
 * Returns a list of consensus haplotype pairs (one pair per individual)
 * sorted in order of increasing sample index. The specified list of
 * haplotype pairs may contain multiple haplotype pairs for each individual.
 *
 * @param hapPairs a list of haplotype pairs
 */
QList<HapPair> ConsensusPhaser::consensusPhase(const QList<HapPair> &hapPairs)
{
  QList<HapPair> copy = hapPairs;

  if (copy.isEmpty())
    return copy;

  checkMarkers(copy);
  qStableSort(copy.begin(), copy.end(), hapsComparator);

  /// Random random = new Random(copy.size());

  QList<HapPair> consensusPairs;
  int start = 0;
  while (start < copy.size()) {
    int end = start + 1;
    while (end < copy.size() && copy[end].sampleIndex() == copy[start].sampleIndex()) {
      ++end;
    }

    if (end - start == 1)
      consensusPairs.append(copy[start]);
    else
      consensusPairs.append(consensus(copy.mid(start, end - start) /*, random */));

    start = end;
  }

  return consensusPairs;
}
