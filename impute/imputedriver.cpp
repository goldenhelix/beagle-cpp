#include "impute/imputedriver.h"

#include "impute/imputationdata.h"

#include <QVector>

#define NON_REFERENCE_WEIGHT 1.0

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

QList<HapPair> ImputeDriver::phase(CurrentData &cd, const Par &par)
{
  QList<HapPair> hapPairs = ImputeDriver::initialHaps(cd, par);

  if (par.burnin_its() > 0) {
    /// runStats.println(Const.nl + "Starting burn-in iterations");
    hapPairs = ImputeDriver::runBurnin1(cd, par, hapPairs);
  }

  if (par.phase40_its() > 0) {
    hapPairs = ImputeDriver::runBurnin2(cd, par, hapPairs);
  }

  if (par.niterations() > 0) {
    /// runStats.println(Const.nl + "Starting phasing iterations");
    /// hapPairs = ImputeDriver::runRecomb(cd, par, hapPairs, gv);
  } else
    hapPairs = ConsensusPhaser::consensusPhase(hapPairs);

  return hapPairs;
}

int ImputeDriver::finishWindow(SampleHapPairs &overlapHaps, const CurrentData &cd, const Par &par,
                               ImputeDataWriter &impWriter, const GLSampleHapPairs &targetHapPairs)
{
  impWriter.printWindowOutput(cd, targetHapPairs, ConstrainedGLAlleleProbs(targetHapPairs), par);

  overlapHaps = ImputeDriver::overlapHaps(cd, targetHapPairs);
  return cd.nMarkers() - cd.nextOverlapStart();
}

int ImputeDriver::finishWindow(SampleHapPairs &overlapHaps, const CurrentData &cd, const Par &par,
                               ImputeDataWriter &impWriter, const SampleHapPairs &targetHapPairs)
{
  // Neither a "printGV" method nor a "refinedIbd" method is invoked here at this time.

  if (cd.nMarkers() == cd.nTargetMarkers() || par.impute() == false)
    impWriter.printWindowOutput(cd, targetHapPairs, ConstrainedAlleleProbs(targetHapPairs), par);
  else {
    impWriter.printWindowOutput(cd, targetHapPairs, ImputeDriver::LSImpute(cd, par, targetHapPairs),
                                par);
  }

  overlapHaps = ImputeDriver::overlapHaps(cd, targetHapPairs);
  return cd.nMarkers() - cd.nextOverlapStart();
}

void ImputeDriver::phaseAndImpute(InputData &data, TargDataReader &targReader,
                                  RefDataReader &refReader, ImputeDataWriter &impWriter,
                                  int windowSize, const Par &par)
{
  CurrentData cd;
  SampleHapPairs overlapHaps;
  int overlap = 0;
  impWriter.writeHeader();
  while (data.canAdvanceWindow(targReader, refReader)) {
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

QList<HapPair> ImputeDriver::runBurnin2(const CurrentData &cd, const Par &par,
                                        QList<HapPair> hapPairs)
{
  QList<HapPair> cumHapPairs;
  int start = par.burnin_its();
  int end = start + par.phase40_its();
  for (int j = start; j < end; ++j) {
    bool useRevDag = (j & 1) == 1;
    hapPairs = ImputeDriver::sample(cd, par, hapPairs, useRevDag);
    // runStats.printIterationUpdate(cd.window(), j+1);
    cumHapPairs.append(hapPairs);
  }
  return cumHapPairs;
}

QList<HapPair> ImputeDriver::sample(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs,
                                    bool useRevDag)
{
  Q_ASSERT_X(!hapPairs.isEmpty(), "ImputeDriver::sample", "hapPairs.isEmpty()");

  int nThreads = par.nThreads();
  int nSampledHaps = par.nSamplingsPerIndividual() * cd.nTargetSamples();
  SplicedGL gl(cd.targetGL(), useRevDag);

  //// Dag dag = getDagsAndUpdatePos(cd, hapPairs, useRevDag);
  //// Dag dag(cd, par.modelscale(), hapPairs, useRevDag));
  ////   private Dag ImputeDriver::getDagsAndUpdatePos(CurrentData cd, List<HapPair> hapPairs,
  ////           boolean useRevDag) {

  cd.addRestrictedRefHapPairs(hapPairs);
  HapPairs dagHaps(hapPairs, useRevDag);

  //// float[] wts = cd.weights().get(dagHaps, cd);
  //// Dag dag = makeDag(dagHaps, wts, par.modelscale());
  //// private Dag ImputeDriver::makeDag(HapPairs hapPairs, float[] weights, float scale) {
  ////   long t0 = System.nanoTime();
  ////   int nInitLevels = 500;
  ////   Dag dag = MergeableDag.dag(hapPairs, weights, scale, nInitLevels);
  ////   runStats.buildNanos(System.nanoTime() - t0);
  ////   return dag;
  //// }
  ////       // runStats.setDagStats(dag);
  ////       return dag;
  ////   }

  ImmutableDag dag(dagHaps, ImputeDriver::getHapWeights(dagHaps, cd), par.modelScale(),
                   par.dagInitLevels());

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
  for (int single = 0; single < nSamples; single++) {
    QList<HapPair> newHaps = baum.randomSample(single);

    if (markersAreReversed) {
      for (int h = 0; h < newHaps.length(); h++)
        sampledHaps.append(HapPair(newHaps[h], true));
    } else
      sampledHaps.append(newHaps);
  }
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
