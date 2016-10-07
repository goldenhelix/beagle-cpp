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


  /**
   * Phases the current window of genotype data. Returns the phased
   * genotype data.
   * @param cd the current window of data
   * @param par the current running parameters
   */
  //// QList<HapPair> ImputeDriver::phase(const CurrentData &cd, const Par &par);

  /**
   * "Lower-level utilities" for phasing.
   */
  QList<HapPair> runBurnin1(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs);
  QList<HapPair> runBurnin2(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs);

  /**
   * Returns a list of sampled haplotype pairs.  Haplotype pairs are
   * sampled conditional on the observed genotype data and a haplotype
   * frequency model in which all markers are in linkage equilibrium.
   *
   * @param cd the input data for the current marker window
   * @param par the current running parameters
   */
  QList<HapPair> initialHaps(CurrentData &cd, const Par &par);

  /**
   * Returns a list of sampled haplotype pairs. Haplotype pairs are
   * sampled conditional on the observed genotype and a haplotype
   * frequency model constructed from the specified {@code hapPairs}.
   * The contract for this method is undefined if the specified
   * {@code hapPairs} is inconsistent with the input data
   * contained in the {@code cd} parameter.
   *
   * @param cd the input data for the current marker window
   * @param par the current running parameters
   * @param hapPairs the haplotype pairs used to build the haplotype
   * frequency model
   * @param useRevDag {@code true} if the order of markers should
   * be reversed when building the haplotype frequency model, and
   * {@code false} otherwise
   */
  QList<HapPair> sample(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs,
                        bool useRevDag);

  /**
   * "Lower-level utility" for performing sampling.
   */
  void sample(Dag &dag, SplicedGL &gl, int seed, bool markersAreReversed,
              int nSamplingsPerIndividual, QList<HapPair> &sampledHaps, int nThreads, bool lowmem);
};

#endif
