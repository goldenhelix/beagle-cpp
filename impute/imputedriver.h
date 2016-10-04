#ifndef IMPUTEDRIVER_H
#define IMPUTEDRIVER_H

#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"
#include "impute/haplotypepair.h"
#include "impute/iointerface.h"
#include "impute/dag.h"
#include "impute/phase.h"

namespace ImputeDriver
{
  SampleHapPairs overlapHaps(const CurrentData &cd, const SampleHapPairs &targetHapPairs);

  //// int windowDriverHelper(SampleHapPairs &overlapHaps, const CurrentData &cd,
  ////                        const SampleHapPairs &targetHapPairs);

  //// void windowDriver(InputData &data, TargDataReader &targReader,
  ////                   RefDataReader &refReader, int windowSize, const Par &par);

  //// QList<HapPair> ImputeDriver::phase(const CurrentData &cd);

  //// QList<HapPair> ImputeDriver::initialHaps(const CurrentData &cd);

  void sample(Dag &dag, SplicedGL &gl, int seed, bool markersAreReversed,
              int nSamples, QList<HapPair> &sampledHaps, int nThreads, bool lowmem);
};

#endif
