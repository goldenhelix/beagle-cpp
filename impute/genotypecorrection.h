/* Copyright notice. */
#ifndef GENOTYPECORRECTION_H
#define GENOTYPECORRECTION_H

#include "impute/haplotypepair.h"

/**
 * {@code GenotypeCorrection::correct()} removes any inconsistencies between
 * haplotype pairs and genotypes that determine genotype likelihoods.
 */
namespace GenotypeCorrection
{
  /**
   * Removes any inconsistencies between the specified list of
   * haplotype pairs and the genotypes determined by the {@code allele1()}
   * and {@code allele2()} methods of the specified genotype likelihoods.
   * Inconsistencies are resolved by changing the minimum number
   * of alleles in the haplotype pairs.
   *
   * @param hapPairs a list of haplotype pairs
   * @param gl genotype likelihoods
   * @param seed a seed for generating random numbers
   *
   */
  void correct(QList<HapPair> &hapPairs, const MaskedEndsGL &gl, int seed);
};

#endif
