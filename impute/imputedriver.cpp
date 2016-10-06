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
  // bool useRevDag = false;
  double minAlleleFreq = 0.0001;
  LinkageEquilibriumDag dag(freqGL, minAlleleFreq);
  QList<HapPair> sampledHaps;
  ImputeDriver::sample(dag, emitGL, par.seed(), false /* useRevDag */,
                       par.nSamplingsPerIndividual(), sampledHaps,
                       par.nThreads(), par.lowMem());
  return sampledHaps;
}

/*
    private List<HapPair> runBurnin1(CurrentData cd, const Par &par, List<HapPair> hapPairs) {
        GenotypeValues gv = null;
        for (int j=0; j<par.burnin_its(); ++j) {
            boolean useRevDag = (j & 1)==1;
            hapPairs = hapSampler.sample(cd, hapPairs, useRevDag, gv);
            runStats.printIterationUpdate(cd.window(), j+1);
        }
        return hapPairs;
    }
*/

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
