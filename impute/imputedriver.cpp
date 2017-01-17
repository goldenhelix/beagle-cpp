#include "impute/imputedriver.h"

#include "impute/imputationdata.h"
#include "impute/genotypecorrection.h"
#include "impute/workerrunner.h"

#include <QVector>
#include <QThread>

#include <stdio.h>

#define NON_REFERENCE_WEIGHT 1.0
#define N_COPIES               4

#define SEND_PROG_MSG(msg, ...) { \
  fprintf(stderr, "PROGRESS_OUTPUT\t" msg "\n", ##__VA_ARGS__); \
  fflush(stderr);  \
}

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
    SEND_PROG_MSG("Loading (marker window) data...");

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
    {
      SEND_PROG_MSG("Outputting (marker window) data...");

      impWriter.printWindowOutput(cd, targetHapPairs, ConstrainedAlleleProbs(targetHapPairs), par);
    }
    else
      impWriter.printWindowOutput(cd, targetHapPairs, ImputeDriver::LSImpute(cd, par, targetHapPairs),
                                  par);
  }

  SEND_PROG_MSG("Finishing marker window...");

  overlapHaps = ImputeDriver::overlapHaps(cd, targetHapPairs);
  return cd.nMarkers() - cd.nextOverlapStart();
}

QList<HapPair> ImputeDriver::phase(CurrentData &cd, const Par &par)
{
  SEND_PROG_MSG("Preparing for initialization...");

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
  float minAlleleFreq = (float) 0.0001;
  LinkageEquilibriumDag dag(freqGL, minAlleleFreq);
  QList<HapPair> sampledHaps;
  ImputeDriver::sample(dag, emitGL, par.seed(), false, par.nSamplingsPerIndividual(), sampledHaps,
                       par.nThreads(), par.lowMem(), "Initializing");
  return sampledHaps;
}

QList<HapPair> ImputeDriver::runBurnin1(const CurrentData &cd, const Par &par,
                                        QList<HapPair> hapPairs)
{
  for (int j = 0; j < par.burnin_its(); ++j) {
    bool useRevDag = (j & 1) == 1;

    char progressBuff[27];
	sprintf(progressBuff, "Burn-in iteration %d of %d", j + 1, par.burnin_its());

    hapPairs = ImputeDriver::sample(cd, par, hapPairs, useRevDag, progressBuff);
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

    char progressBuff[33];
	sprintf(progressBuff, "Phasing (4.0) iteration %d of %d", j + 1 - start, par.phase40_its());

    hapPairs = ImputeDriver::sample(cd, par, hapPairs, useRevDag, progressBuff);
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

    char progressBuff[33];
	sprintf(progressBuff, "Phasing (4.1) iteration %d of %d", j + 1 - start, par.niterations());

    hapPairs = ImputeDriver::recombSample(cd, par, hapPairs, useRevDag, progressBuff);
    cumHapPairs.append(hapPairs);
  }
  hapPairs = ConsensusPhaser::consensusPhase(cumHapPairs);
  hapPairs = correctGenotypes(cd, par, hapPairs);
  return hapPairs;
}

QList<HapPair> ImputeDriver::sample(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs,
                                    bool useRevDag, char *whichIteration)
{
  Q_ASSERT_X(!hapPairs.isEmpty(), "ImputeDriver::sample", "hapPairs.isEmpty()");

  SEND_PROG_MSG("%s: Preparing directed acyclic graph...", whichIteration);

  int nThreads = par.nThreads();

  SplicedGL gl(cd.targetGL(), useRevDag);

  cd.addRestrictedRefHapPairs(hapPairs);
  HapPairs dagHaps(hapPairs, useRevDag);

  ImmutableDag dag(dagHaps, ImputeDriver::getHapWeights(dagHaps, cd), par.modelScale(),
                   par.dagInitLevels());

  QList<HapPair> sampledHaps;

  sample(dag, gl, par.seed(), useRevDag, par.nSamplingsPerIndividual(), sampledHaps, nThreads,
         par.lowMem(), whichIteration);

  return sampledHaps;
}

QVector<float> ImputeDriver::getHapWeights(HapPairs haps, const CurrentData &cd)
{
  //// Samples samples = families().samples();
  Samples samples = cd.targetSamples();
  QVector<float> fa(haps.nHaps());

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

    float sampleWeight = (sampleIndex == -1) ? 1.0 : NON_REFERENCE_WEIGHT;
    int cnt = cntMap[idIndex];
    float wt = sampleWeight / cnt;

    //// float MIN_SAMPLE_WEIGHT = 0.01f;
    //// float minWt = MIN_SAMPLE_WEIGHT/cnt;
    //// fa[hapIndex++] = parentCnt>0  ? minWt : wt;
    //// fa[hapIndex++] = parentCnt==2 ? minWt : wt;

    fa[hapIndex++] = wt;
    fa[hapIndex++] = wt;
  }

  return fa;
}

    /*
struct WorkerArgs {
  Dag* dag;
  SplicedGL* gl;
  int seed;
  int nSamplingsPerIndividual;
  bool lowmem;
  int sampleIdx;

  WorkerArgs(Dag* dag_, SplicedGL* gl_, int seed_, int nSamplingsPerIndividual_, bool lowmem_, int sampleIdx_)
    : dag(dag_), gl(gl_), seed(seed_), nSamplingsPerIndividual(nSamplingsPerIndividual_), lowmem(lowmem_), sampleIdx(sampleIdx_)
  {}
};

QList<HapPair> runSingleBaumWork(WorkerArgs args)
{
 SingleBaum baum(*args.dag, *args.gl, args.seed, args.nSamplingsPerIndividual, args.lowmem);
 return baum.randomSample(single);
}
    */

QList<QThread*> gThreads; // Global thread pool

void ImputeDriver::sample(Dag &dag, SplicedGL &gl, int seed, bool markersAreReversed,
                          int nSamplingsPerIndividual, QList<HapPair> &sampledHaps, int nThreads,
                          bool lowmem, char *whichIteration)
{
  for (int i = gThreads.size(); i < nThreads; i++) {
    QThread *t = new QThread();
    t->setObjectName("Imputation Worker Thread");
    t->start();
    gThreads << t;
  }

  SingleBaumRunner runner;

  // One worker per thread
  for (int i = 0; i < nThreads; i++) {
    SingleBaumWorker *worker = new SingleBaumWorker(dag, gl, seed, nSamplingsPerIndividual, lowmem);
    runner.moveWorkerToThread(worker, gThreads[i]);
  }

  int nSamples = gl.nSamples();
  runner.start(nSamples);

  while (runner.hasNextWorkerResult()) {
    QList<HapPair> newHaps = runner.nextWorkerResult();
    int curSample = runner.nextResultIdx();
    SEND_PROG_MSG("%s: sample %d of %d...", whichIteration, curSample, nSamples);

    if (markersAreReversed) {
      for (int h = 0; h < newHaps.length(); h++)
        sampledHaps.append(HapPair(newHaps[h], true));
    } else
      sampledHaps.append(newHaps);
  }
  // runStats.sampleNanos(System.nanoTime() - t0);
}

QList<HapPair> ImputeDriver::recombSample(const CurrentData &cd, const Par &par,
                                          const QList<HapPair> &hapPairs,
                                          bool useRevDag, char *whichIteration)
{
  /*
  SEND_PROG_MSG("%s: Preparing DAG and IBS data...", whichIteration);

  QList<HapPair> haps = ConsensusPhaser.consensusPhase(hapPairs);
  cd.addRestrictedRefHapPairs(haps);

  SampleHapPairs dagHaps(cd.allSamples(), haps, useRevDag);

  ImmutableDag dag(dagHaps, ImputeDriver::getHapWeights(dagHaps, cd), par.modelScale(),
                   par.dagInitLevels());

  RestrictedDag rdag(dag, dagHaps, par.ibdlength(), par.ibdextend());

  SamplerData samplerData(rdag, par, cd, useRevDag /+ , runStats +/ );

  QList<HapPair> sampledHaps;
  ImputeDriver::recombSample(samplerData, par, sampledHaps, whichIteration);
  return sampledHaps;
  */
  QList<HapPair> retHaps = hapPairs;
  return retHaps;
}
/*
void ImputeDriver::recombSample(const SamplerData &samplerData, const Par &par,
                                QList<HapPair> &sampledHaps, char *whichIteration)
{
  // long t0 = System.nanoTime();
  // int nThreads = par().nthreads();
  bool markersAreReversed = samplerData.markersAreReversed();
  // Random rand = new Random(par.seed());

  RecombSingleBaum baum(samplerData, /+ rand.nextLong(), +/ N_COPIES, par.lowmem());

  int nSamples = samplerData.nSamples();
  for (int single = 0; single < nSamples; single++)
  {
    SEND_PROG_MSG("%s: sample %d of %d...", whichIteration, single + 1, nSamples);

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
*/
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
  SEND_PROG_MSG("Preparing to impute (marker window) data...");

  // Should expect our gThreads initialized by previous ::sample calls
  Q_ASSERT(gThreads.size() > 0 && gThreads.size() == par.nThreads());

  LSImputeRunner runner;

  double scaleFactor = 1e-6;

  // One worker per thread
  for (int i = 0; i < par.nThreads(); i++) {
    LSImputeWorker* worker = new LSImputeWorker(par, cd, targetHapPairs, scaleFactor);
    runner.moveWorkerToThread(worker, gThreads[i]);
  }

  int nSamples = targetHapPairs.nSamples();
  runner.start(nSamples);

  QList<HapAlleleProbs> hapAlProbList;
  while (runner.hasNextWorkerResult()) {
    QList<HapAlleleProbs> haps = runner.nextWorkerResult();
    int curSample = runner.nextResultIdx();
    SEND_PROG_MSG("Imputing for sample %d of %d...", curSample, nSamples);

    hapAlProbList << haps;
  }
/*
  double scaleFactor = 1e-6;
  PositionMap imputationMap(scaleFactor);

  ImputationData impData(par, cd, targetHapPairs, imputationMap);

  LSHapBaum hb(impData, par.lowMem());

  QList<HapAlleleProbs> hapAlProbList;

  for (int sample = 0, n = targetHapPairs.nSamples(); sample < n; ++sample) {
    int hap1 = 2 * sample;
    int hap2 = 2 * sample + 1;

    SEND_PROG_MSG("Imputing for sample %d of %d...", sample + 1, n);

    hapAlProbList.append(hb.randomHapSample(hap1));
    hapAlProbList.append(hb.randomHapSample(hap2));
  }
*/
  SEND_PROG_MSG("Outputting (marker window) data...");

  return ConstrainedAlleleProbs(targetHapPairs, hapAlProbList, cd.markerIndices());
}
