#ifndef IMPUTEDRIVER_H
#define IMPUTEDRIVER_H

#include "impute/dag.h"
#include "impute/hapalleleprobs.h"
#include "impute/haplotypepair.h"
#include "impute/iointerface.h"

class SamplerData;

namespace ImputeDriver
{
  void phaseAndImpute(InputData &data, TargDataReader &targReader, RefDataReader &refReader,
                      ImputeDataWriter &impWriter, int windowSize, const Par &par);

  int finishWindow(SampleHapPairs &overlapHaps, const CurrentData &cd, const Par &par,
                   ImputeDataWriter &impWriter, const SampleHapPairs &targetHapPairs);

  /**
   * Phases the current window of genotype data. Returns the phased
   * genotype data.
   * @param cd the current window of data
   * @param par the current running parameters
   */
  QList<HapPair> phase(CurrentData &cd, const Par &par);

  SampleHapPairs overlapHaps(const CurrentData &cd, const SampleHapPairs &targetHapPairs);

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
   * "Lower-level" utilities for phasing.
   */
  QList<HapPair> runBurnin1(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs);
  QList<HapPair> runBurnin2(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs);
  QList<HapPair> runRecomb(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs);

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
                        bool useRevDag, char *whichIteration);

  /**
   * Returns an array of length {@code haps.nHaps()} with
   * per-haplotype weights. Array elements {@code 2*j} and {@code 2*j
   * + 1} are the weights for the first and second haplotype in the
   * {@code j}-th haplotype pair.
   *
   * From Weights.get(): "Reference haplotypes are assigned a weight
   * of {@code 1.0f}.  Non-reference haplotypes are assigned a weight
   * of {@code this.nonRefWt()} if the haplotype is not inherited from
   * a parent in the sample, and a weight of {@code 0.01f} if the
   * haplotype is inherited from a parent in the sample.  The first
   * haplotype in the offspring is required to be the transmitted
   * transmitted haplotype for a parent-offspring duo." (Should have
   * then said that each weight is afterward divided by one plus the
   * count of how many other haplotype pairs are assigned to the same
   * sample.)
   *
   * Here, reference haplotypes are simply assigned a tentative weight
   * of {@code 1.0f} and non-reference haplotypes are assigned a
   * tentative weight of NON_REFERENCE_WEIGHT. Then, each is divided
   * by one plus the count of how many other haplotype pairs are
   * assigned to the same sample.
   *
   * @param haps an array of haplotype pairs
   * @param cd the input data for the current marker window
   */
  QVector<float> getHapWeights(HapPairs haps, const CurrentData &cd);

  /**
   * "Lower-level utility" for performing sampling.
   */
  void sample(const Dag &dag, const SplicedGL &gl, int seed, bool markersAreReversed,
              int nSamplingsPerIndividual, QList<HapPair> &sampledHaps,
              int nThreads, bool lowmem, char *whichIteration);

  /**
   * Returns a list of sampled haplotype pairs. Haplotype pairs are
   * sampled conditional on the observed genotype and a haplotype
   * frequency model constructed from the specified {@code hapPairs}.
   * The contract for this method is undefined if the specified
   * {@code hapPairs} and {@code gv} are inconsistent with the input data
   * contained in the {@code cd} parameter.
   *
   * @param cd the input data for the current marker window
   * @param par the current parameters
   * @param hapPairs the target haplotype pairs used to build the haplotype
   * frequency model
   * @param useRevDag {@code true} if the order of markers should
   * be reversed when building the haplotype frequency model, and
   * {@code false} otherwise
   */
  QList<HapPair> recombSample(const CurrentData &cd, const Par &par,
                              const QList<HapPair> &hapPairs,
                              bool useRevDag, char *whichIteration);

  /**
   * "Lower-level utility" for performing (phasing 4.1) sampling.
   */
  void recombSample(const SamplerData &samplerData, const Par &par,
                    QList<HapPair> &sampledHaps, char *whichIteration);

  QList<HapPair> correctGenotypes(const CurrentData &cd, const Par &par, QList<HapPair> hapPairs);

  ConstrainedAlleleProbs LSImpute(const CurrentData &cd, const Par &par,
                                  const SampleHapPairs &targetHapPairs);
};

#endif
