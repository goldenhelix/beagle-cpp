#include "impute/imputedriver.h"

#include "impute/imputationdata.h"

#include <QVector>

#define NON_REFERENCE_WEIGHT 1.0
#define N_COPIES               4

void ImputeDriver::phaseAndImpute(InputData &data, TargDataReader &targReader,
                                  RefDataReader &refReader, ImputeDataWriter &impWriter,
                                  int windowSize, const Par &par)
{
  CurrentData cd;
  SampleHapPairs overlapHaps;
  int overlap = 0;
  impWriter.writeHeader();
  while (data.canAdvanceWindow(targReader, refReader))
  {
    data.advanceWindow(overlap, par.window(), targReader, refReader);
    data.setCdData(cd, par, overlapHaps, targReader, refReader);

    if (cd.targetGL().isRefData())
      overlap = ImputeDriver::finishWindow(overlapHaps, cd, par, impWriter,
                                           GLSampleHapPairs(cd.targetGL(), true));
    else {
      QList<HapPair> hapPairs = ImputeDriver::phase(cd, par);
      overlap = ImputeDriver::finishWindow(overlapHaps, cd, par, impWriter,
                                           SampleHapPairs(cd.targetSamples(), hapPairs, false));
    }
  }
  impWriter.writeEOF();
}

int ImputeDriver::finishWindow(SampleHapPairs &overlapHaps, const CurrentData &cd, const Par &par,
                               ImputeDataWriter &impWriter, const SampleHapPairs &targetHapPairs)
{
  // Neither a "printGV" method nor a "refinedIbd" method is invoked here at this time.

  if(cd.nTargetMarkers())
  {
    if (cd.nMarkers() == cd.nTargetMarkers() || par.impute() == false)
      impWriter.printWindowOutput(cd, targetHapPairs, ConstrainedAlleleProbs(targetHapPairs), par);
    else
      impWriter.printWindowOutput(cd, targetHapPairs, ImputeDriver::LSImpute(cd, par, targetHapPairs),
                                  par);
  }

  overlapHaps = ImputeDriver::overlapHaps(cd, targetHapPairs);
  return cd.nMarkers() - cd.nextOverlapStart();
}

QList<HapPair> ImputeDriver::phase(CurrentData &cd, const Par &par)
{
  QList<HapPair> hapPairs = ImputeDriver::initialHaps(cd, par);

  if (par.burnin_its() > 0)
  {
    /// runStats.println(Const.nl + "Starting burn-in iterations");
    hapPairs = ImputeDriver::runBurnin1(cd, par, hapPairs);
  }

  if (par.phase40_its() > 0)
  {
    hapPairs = ImputeDriver::runBurnin2(cd, par, hapPairs);
  }

  //// if (par.niterations() > 0)
  //// {
  ////   /// runStats.println(Const.nl + "Starting phasing iterations");
  ////   hapPairs = ImputeDriver::runRecomb(cd, par, hapPairs);
  //// } else
    hapPairs = ConsensusPhaser::consensusPhase(hapPairs);

  return hapPairs;
}

SampleHapPairs ImputeDriver::overlapHaps(const CurrentData &cd,
                                         const SampleHapPairs &targetHapPairs)
{
  int nextOverlap = cd.nextTargetOverlapStart();
  int nextSplice = cd.nextTargetSpliceStart();
  if (cd.nextOverlapStart() == cd.nextSpliceStart()) {
    return SampleHapPairs();  // Default constructor creates an empty result.
  }

  Markers markers = targetHapPairs.markers().restrict(nextOverlap, nextSplice);
  Samples samples = targetHapPairs.samples();

  QList<int> mapping;
  for (int i = nextOverlap; i < nextSplice; i++)
    mapping.append(i);
  QList<HapPair> list = HapUtility::createHapPairList(markers, targetHapPairs, mapping);

  return SampleHapPairs(targetHapPairs.samples(), list, false);
}

QList<HapPair> ImputeDriver::initialHaps(CurrentData &cd, const Par &par)
{
  SplicedGL freqGL = cd.targetGL();
  SplicedGL emitGL = cd.targetGL();
  double minAlleleFreq = 0.0001;
  LinkageEquilibriumDag dag(freqGL, minAlleleFreq);
  QList<HapPair> sampledHaps;
  ImputeDriver::sample(dag, emitGL, par.seed(), false, par.nSamplingsPerIndividual(), sampledHaps,
                       par.nThreads(), par.lowMem());
  return sampledHaps;
}

QList<HapPair> ImputeDriver::runBurnin1(const CurrentData &cd, const Par &par,
                                        QList<HapPair> hapPairs)
{
  for (int j = 0; j < par.burnin_its(); ++j) {
    bool useRevDag = (j & 1) == 1;
    hapPairs = ImputeDriver::sample(cd, par, hapPairs, useRevDag);
    // runStats.printIterationUpdate(cd.window(), j+1);
  }
  return hapPairs;
}
static int b2j;
QList<HapPair> ImputeDriver::runBurnin2(const CurrentData &cd, const Par &par,
                                        QList<HapPair> hapPairs)
{
  QList<HapPair> cumHapPairs;
  int start = par.burnin_its();
  int end = start + par.phase40_its();
  for (b2j = start; b2j < end; ++b2j) {
    bool useRevDag = (b2j & 1) == 1;
    hapPairs = ImputeDriver::sample(cd, par, hapPairs, useRevDag);
    // runStats.printIterationUpdate(cd.window(), b2j+1);
    cumHapPairs.append(hapPairs);
  }
  return cumHapPairs;
}

QList<HapPair> ImputeDriver::runRecomb(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs)
{
  hapPairs = ConsensusPhaser::consensusPhase(hapPairs);
  QList<HapPair> cumHapPairs;
  int start = par.burnin_its() + par.phase40_its();
  int end = start + par.niterations();
  for (int j=start; j<end; ++j)
  {
    bool useRevDag = (j & 1)==1;
    hapPairs = ImputeDriver::recombSample(cd, par, hapPairs, useRevDag);
    // runStats.printIterationUpdate(cd.window(), j+1);
    cumHapPairs.append(hapPairs);
  }
  hapPairs = ConsensusPhaser::consensusPhase(cumHapPairs);
  hapPairs = correctGenotypes(cd, par, hapPairs);
  return hapPairs;
}

QList<HapPair> ImputeDriver::sample(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs,
                                    bool useRevDag)
{
  Q_ASSERT_X(!hapPairs.isEmpty(), "ImputeDriver::sample", "hapPairs.isEmpty()");

  int nThreads = par.nThreads();
  int nSampledHaps = par.nSamplingsPerIndividual() * cd.nTargetSamples();
  SplicedGL gl(cd.targetGL(), useRevDag);

  // if(b2j == 6)
  //   HapUtility::dumpGl(gl);

  // if(b2j == 7)
  //   HapUtility::dumpHp(hapPairs);

  cd.addRestrictedRefHapPairs(hapPairs);
  HapPairs dagHaps(hapPairs, useRevDag);

  if(b2j == 7)
    globalMgd.setOkToDump(true);
  else if(b2j == 8)
    globalMgd.setOkToDump(false);
  ImmutableDag dag(dagHaps, ImputeDriver::getHapWeights(dagHaps, cd), par.modelScale(),
                   par.dagInitLevels());

  // if(b2j == 7)
  //   DagDump::dagDump(dag);

  QList<HapPair> sampledHaps;

  sample(dag, gl, par.seed(), useRevDag, par.nSamplingsPerIndividual(), sampledHaps, nThreads,
         par.lowMem());

  return sampledHaps;
}

QVector<double> ImputeDriver::getHapWeights(HapPairs haps, const CurrentData &cd)
{
  //// Samples samples = families().samples();
  Samples samples = cd.targetSamples();
  QVector<double> fa(haps.nHaps());

  int nHapPairs = haps.nHapPairs();
  QMap<int, int> cntMap;
  for (int j = 0; j < nHapPairs; ++j) {
    int idIndex = haps.samples(j).idIndex(haps.sampleIndex(j));

    // "cntMap" defaults to zero if there is no value already assigned
    // for idIndex.

    cntMap[idIndex] = cntMap[idIndex] + 1;
  }

  int hapIndex = 0;
  for (int j = 0, n = nHapPairs; j < n; ++j) {
    int idIndex = haps.samples(j).idIndex(haps.sampleIndex(j));
    int sampleIndex = samples.findLocalIndex(idIndex);

    //// int parentCnt = 0;
    //// if (sampleIndex != -1)
    //// {
    ////   // sample is a non-reference sample
    ////   if (families().father(sampleIndex)>=0) {
    ////     ++parentCnt;
    ////   }
    ////   if (families().mother(sampleIndex)>=0) {
    ////     ++parentCnt;
    ////   }
    //// }

    double sampleWeight = (sampleIndex == -1) ? 1.0 : NON_REFERENCE_WEIGHT;
    int cnt = cntMap[idIndex];
    double wt = sampleWeight / cnt;

    //// float MIN_SAMPLE_WEIGHT = 0.01f;
    //// float minWt = MIN_SAMPLE_WEIGHT/cnt;
    //// fa[hapIndex++] = parentCnt>0  ? minWt : wt;
    //// fa[hapIndex++] = parentCnt==2 ? minWt : wt;

    fa[hapIndex++] = wt;
    fa[hapIndex++] = wt;
  }

  return fa;
}

void ImputeDriver::sample(Dag &dag, SplicedGL &gl, int seed, bool markersAreReversed,
                          int nSamplingsPerIndividual, QList<HapPair> &sampledHaps, int nThreads,
                          bool lowmem)
{
  SingleBaum baum(dag, gl, seed, nSamplingsPerIndividual, lowmem);

  int nSamples = gl.nSamples();
  for (int single = 0; single < nSamples; single++)
  {
    QList<HapPair> newHaps = baum.randomSample(single);

    if (markersAreReversed)
    {
      for (int h = 0; h < newHaps.length(); h++)
        sampledHaps.append(HapPair(newHaps[h], true));
    } else
      sampledHaps.append(newHaps);
  }

  // runStats.sampleNanos(System.nanoTime() - t0);
}

QList<HapPair> ImputeDriver::recombSample(const CurrentData &cd, const Par &par,
                                          const QList<HapPair> &hapPairs,
                                          bool useRevDag)
{
  SamplerData samplerData(par, cd, hapPairs, useRevDag /* , runStats */ );
  int nSampledHaps = N_COPIES * cd.nTargetSamples();
  QList<HapPair> sampledHaps;
  ImputeDriver::recombSample(samplerData, par, sampledHaps);
  return sampledHaps;
}

void ImputeDriver::recombSample(const SamplerData &samplerData, const Par &par, QList<HapPair> &sampledHaps)
{
  // long t0 = System.nanoTime();
  // int nThreads = samplerData.par().nthreads();
  bool markersAreReversed = samplerData.markersAreReversed();
  // Random rand = new Random(par.seed());

  RecombSingleBaum baum(samplerData, /* rand.nextLong(), */ N_COPIES, true /* par.lowmem() */);

  int nSamples = samplerData.nSamples();
  for (int single = 0; single < nSamples; single++)
  {
    QList<HapPair> newHaps = baum.randomSample(single);

    if (markersAreReversed)
    {
      for (int h = 0; h < newHaps.length(); h++)
        sampledHaps.append(HapPair(newHaps[h], true));
    } else
      sampledHaps.append(newHaps);
  }

  // runStats.sampleNanos(System.nanoTime() - t0);
}

//// Fake definition standing in for Version 4.1 phasing:
QList<HapPair> GenotypeCorrection::correct(QList<HapPair> hapPairs, const MaskedEndsGL &gl, int seed)
{
  QList<HapPair> qlhp; return qlhp;
}

QList<HapPair> ImputeDriver::correctGenotypes(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs)
{
  int start = cd.prevTargetSpliceStart();
  int end = cd.nextTargetSpliceStart();
  MaskedEndsGL modGL(cd.targetGL(), start, end);
  GenotypeCorrection::correct(hapPairs, modGL, par.seed());
  return hapPairs;
}

ConstrainedAlleleProbs ImputeDriver::LSImpute(const CurrentData &cd, const Par &par,
                                              const SampleHapPairs &targetHapPairs)
{
  double scaleFactor = 1e-6;
  PositionMap imputationMap(scaleFactor);

  ImputationData impData = ImputationData(par, cd, targetHapPairs, imputationMap);

  LSHapBaum hb(impData, true);

  QList<HapAlleleProbs> hapAlProbList;

  for (int sample = 0, n = targetHapPairs.nSamples(); sample < n; ++sample) {
    int hap1 = 2 * sample;
    int hap2 = 2 * sample + 1;
    hapAlProbList.append(hb.randomHapSample(hap1));
    hapAlProbList.append(hb.randomHapSample(hap2));
  }

  return ConstrainedAlleleProbs(targetHapPairs, hapAlProbList, cd.markerIndices());
}
