#include "impute/imputedriver.h"

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

//// int ImputeDriver::windowDriverHelper(SampleHapPairs &overlapHaps, const CurrentData &cd,
////                                      const Par &par, const SampleHapPairs &targetHapPairs)
//// { etc.; }

//// void ImputeDriver::windowDriver(InputData &data, TargDataReader &targReader,
////                                 RefDataReader &refReader, int windowSize, const Par &par)
//// { etc.; }

//// QList<HapPair> ImputeDriver::phase(const CurrentData &cd, const Par &par) { etc.; }

QList<HapPair> ImputeDriver::initialHaps(CurrentData &cd, const Par &par)
{
  SplicedGL freqGL = cd.targetGL();
  SplicedGL emitGL = cd.targetGL();
  double minAlleleFreq = 0.0001;
  LinkageEquilibriumDag dag(freqGL, minAlleleFreq);
  QList<HapPair> sampledHaps;
  ImputeDriver::sample(dag, emitGL, par.seed(), false,
                       par.nSamplingsPerIndividual(), sampledHaps,
                       par.nThreads(), par.lowMem());
  return sampledHaps;
}

QList<HapPair> ImputeDriver::runBurnin1(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs)
{
  for (int j=0; j < par.burnin_its(); ++j)
  {
    bool useRevDag = (j & 1)==1;
    hapPairs = ImputeDriver::sample(cd, par, hapPairs, useRevDag);
    // runStats.printIterationUpdate(cd.window(), j+1);
  }
  return hapPairs;
}

QList<HapPair> ImputeDriver::runBurnin2(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs)
{
  QList<HapPair> cumHapPairs;
  int start = par.burnin_its();
  int end = start + par.phase40_its();
  for (int j=start; j<end; ++j)
  {
    bool useRevDag = (j & 1)==1;
    hapPairs = ImputeDriver::sample(cd, par, hapPairs, useRevDag);
    // runStats.printIterationUpdate(cd.window(), j+1);
    cumHapPairs.append(hapPairs);
  }
  return cumHapPairs;
}

QList<HapPair> ImputeDriver::sample(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs,
                                    bool useRevDag)
{
  Q_ASSERT_X(!hapPairs.isEmpty(),
             "ImputeDriver::sample",
             "hapPairs.isEmpty()");

  int nThreads =  par.nThreads();
  int nSampledHaps = par.nSamplingsPerIndividual() * cd.nTargetSamples();
  SplicedGL gl(cd.targetGL(), useRevDag);
  /////////// Dag dag = getDagsAndUpdatePos(cd, hapPairs, useRevDag);
  /////////// Dag dag(cd, par.modelscale(), hapPairs, useRevDag));
  QList<HapPair> sampledHaps;
  // sample(dag, gl, par.seed(), useRevDag, par.nSamplingsPerIndividual(),
  //        sampledHaps, nThreads, par.lowMem());

  return sampledHaps;
}

void ImputeDriver::sample(Dag &dag, SplicedGL &gl, int seed,
                          bool markersAreReversed, int nSamplingsPerIndividual,
                          QList<HapPair> &sampledHaps, int nThreads, bool lowmem)
{
  SingleBaum baum(dag, gl, seed, nSamplingsPerIndividual, lowmem);

  int nSamples = gl.nSamples();
  for(int single=0; single < nSamples; single++)
  {
    QList<HapPair> newHaps = baum.randomSample(single);

    if (markersAreReversed)
    {
      for(int h=0; h < newHaps.length(); h++)
        sampledHaps.append(HapPair(newHaps[h], true));
    }
    else
      sampledHaps.append(newHaps);
  }
}
