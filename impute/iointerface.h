/* Copyright notice.... */
#ifndef IOUTILITIES_H
#define IOUTILITIES_H

#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"
#include "impute/haplotypepair.h"

#include <QVector>

class Par
{
public:
  virtual int window() const { return 50000; }
  virtual int overlap() const { return 3000; }
  virtual int nThreads() const { return 4; }
  virtual int nSamplingsPerIndividual() const { return 4; }
  virtual int seed() const { return 12345; }
  virtual bool lowMem() const { return true; }
  virtual int burnin_its() const { return 5; }
  virtual int phase40_its() const { return 5; }
  virtual int niterations() const { return 5; }
  virtual int dagInitLevels() const { return 500; }
  virtual float modelScale() const { return (float) 0.8; }
  virtual bool impute() const { return true; }
  virtual bool gprobs() const { return false; }
  virtual float cluster() const { return (float) 0.005; }
  virtual float ne() const { return 1000000.0; }
  virtual float err() const { return (float) 0.0001; }
};

class GenericDataReader
{
public:
  virtual bool canAdvanceWindow() const = 0;
  virtual void makeNewWindow(int overlap) = 0;
  virtual void addNewDataToNewWindow(int windowSize) = 0;
  virtual int windowSize() const = 0;
  virtual bool lastWindowOnChrom() const = 0;

  Samples samples() const { return _samples; }
  Markers markers() const { return _markers; }
protected:
  Samples _samples;
  Markers _markers;
};

class RefDataReader : public GenericDataReader
{
public:
  void makeNewWindow(int overlap);

  int windowSize() const { return _vcfRefRecs.length(); }
  QList<BitSetRefGT> refRecs() const { return _vcfRefRecs; }

  virtual bool canAdvanceWindow() const { return false; }

  virtual void addNewDataToNewWindow(int windowSize);      // Re-implementation of this method is optional.
  virtual bool lastWindowOnChrom() const;                  // Re-implementation of this method is optional.

protected:
  virtual bool hasNextRec() const { return false; }
  virtual BitSetRefGT nextRec() const { return BitSetRefGT(); }
  virtual void advanceRec() {}

  bool sameChrom(const BitSetRefGT &a, const BitSetRefGT &b) const;
  bool samePosition(const BitSetRefGT &a, const BitSetRefGT &b) const;
  int currentChromIndex() const;

  QList<BitSetRefGT> _oldVcfRefRecs;
  QList<BitSetRefGT> _vcfRefRecs;
};

class TargDataReader : public GenericDataReader
{
public:
  TargDataReader();

  void makeNewWindow(int overlap);

  QList<int> restrictedMakeNewWindow(const QList<int> &oldRefIndices, int overlap);
  void restrictedAdvanceWindow(QList<int> &refIndices, int refOverlap, const Markers &nextMarkers);

  int windowSize() const { return _vcfEmissions.length(); }
  QList<BitSetGT> vcfRecs() const { return _vcfEmissions; }
  int restrictedCumMarkerCount() { return _restrictedCumMarkerCnt; }

  virtual bool canAdvanceWindow() const = 0;

  virtual void addNewDataToNewWindow(int windowSize);      // Re-implementation of this method is optional.
  virtual bool lastWindowOnChrom() const;                  // Re-implementation of this method is optional.

  int restrictedSingleMarkerCnt() { return _restrictedSingleMarkerCnt; }
  int restrictedZeroMarkerCnt() { return _restrictedZeroMarkerCnt; }

protected:
  virtual bool hasNextRec() const = 0;
  virtual BitSetGT nextRec() const = 0;
  virtual void advanceRec() = 0;

  bool sameChrom(const BitSetGT &a, const BitSetGT &b) const;
  bool samePosition(const BitSetGT &a, const BitSetGT &b) const;
  int currentChromIndex() const;

  QList<BitSetGT> _oldVcfEmissions;
  QList<BitSetGT> _vcfEmissions;

  int _restrictedCumMarkerCnt;
  int _restrictedSingleMarkerCnt;
  int _restrictedZeroMarkerCnt;
};

class VcfWindow
{
public:
  VcfWindow() : _overlap(0), _cumMarkerCnt(0) {}
  /**
   * Advances the sliding window of VCF records. The advanced window
   * (a list of {@code BitSetRefGT} or of {@code BitSetGT} objects) is
   * kept in the DataReader object.  The size of the advanced window
   * and the number of markers of overlap between the marker window
   * immediately before method invocation and the marker window
   * immediately after method invocation may differ from the requested
   * values.  If the advanced window size or overlap is less than the
   * requested value, the actual value will be as large as
   * possible. If {@code this.lastWindowOnChrom() == true} before
   * method invocation, then there will be no overlap between the
   * advanced window and the previous window.
   *
   * @param overlap the number of markers of overlap
   * @param windowSize the requested number of the markers in the window
   * immediately after the method returns
   * @param the data reader--either {@code RefDataReader} or {@code
   * TargetDataReader}. Access is through virtual methods of their
   * common base class, {@code GenericDataReader}.
   */
  void advanceWindow(int overlap, int windowSize, GenericDataReader &dr);

  int overlap() { return _overlap; }
private:
  void checkParameters(int overlap, int windowSize, GenericDataReader &dr);

  int _overlap;
  int _cumMarkerCnt;
};

class CurrentData;

class InputData
{
public:
  InputData() : _window(0) {}
  virtual bool canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const = 0;
  virtual void advanceWindow(int overlap, int windowSize, TargDataReader &tr,
                             RefDataReader &rr) = 0;
  virtual void setCdData(CurrentData &cd, const Par &par, const SampleHapPairs &overlapHaps,
                         const TargDataReader &tr, const RefDataReader &rr) = 0;

protected:
  /* Returns the index of the first marker in the overlap */
  int nextOverlapStart(int targetOverlap, const GenericDataReader &gr);

  VcfWindow _vcfWindow;
  int _window;
  SplicedGL _gl;
};

class TargetData : public InputData
{
public:
  TargetData() : InputData() {}
  bool canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const;
  void advanceWindow(int overlap, int windowSize, TargDataReader &tr, RefDataReader &rr);
  void setCdData(CurrentData &cd, const Par &par, const SampleHapPairs &overlapHaps,
                 const TargDataReader &tr, const RefDataReader &rr);
};

class AllData : public InputData
{
public:
  AllData() : InputData() { _refIndices.clear(); }
  bool canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const;
  void advanceWindow(int overlap, int windowSize, TargDataReader &tr, RefDataReader &rr);
  void setCdData(CurrentData &cd, const Par &par, const SampleHapPairs &overlapHaps,
                 const TargDataReader &tr, const RefDataReader &rr);

private:
  Samples allSamples(const TargDataReader &tr, const RefDataReader &rr);
  void checkSampleOverlap(Samples ref, Samples nonRef);

  RefHapPairs _refSampleHapPairs;

  QList<int> _refIndices;
  // QList<int> _targetIndices;
  Samples _allSamples;
  QList<HapPair> _refHapPairs;
  QList<HapPair> _targetRefHapPairs;
};

class CurrentData
{
  friend class TargetData;  // For initialization purposes.
  friend class AllData;

public:
  CurrentData() {}
  /**
   * Returns the marker window index.
   */
  int window() const { return _window; }
  /**
   * Returns the first marker index in the overlap between this
   * marker window and the next marker window, or
   * returns {@code this.nMarkers()} there is no overlap.
   */
  int nextOverlapStart() const { return _nextOverlapStart; }
  /**
   * Returns the first target marker index in the overlap between this
   * marker window and the next marker window, or
   * returns {@code this.nMarkers()} if there is no overlap or if there are
   * no target markers in the overlap.
   */
  int nextTargetOverlapStart() const { return _nextTargetOverlapStart; }
  /**
   * Returns the first marker index after the splice point with
   * the previous marker window. Returns 0 if the current marker window
   * is the first marker window.
   */
  int prevSpliceStart() const { return _prevSpliceStart; }
  /**
   * Returns the first marker index after the splice point between this
   * marker window and the next marker window, or returns
   * {@code this.nMarkers()} if there is no overlap or if there are
   * no markers after the splice point.
   */
  int nextSpliceStart() const { return _nextSpliceStart; }
  /**
   * Returns the first target marker index after the splice point with
   * the previous marker window. Returns 0 if the current marker window
   * is the first marker window.
   */
  int prevTargetSpliceStart() const { return (_initHaps.nHaps() == 0) ? 0 : _initHaps.nMarkers(); }
  /**
   * Returns the first target marker index after the splice point between this
   * marker window and the next marker window, or returns
   * {@code this.nTargetMarkers()} if there is no overlap or if there are
   * no target markers after the splice point
   */
  int nextTargetSpliceStart() const { return _nextTargetSpliceStart; }
  /**
   * Returns the target data haplotype pairs in the segment of the current
   * marker window preceding the splice point with the previous marker window:
   * {@code this.targetMarkers().restrict(0, this.prevTargetSplice())}
   */
  SampleHapPairs initHaps() const { return _initHaps; }
  /**
   * Returns the parent-offspring relationships.
   */
  // NuclearFamilies families() {
  //     return _families;
  // }

  /**
   * Returns the per-haplotype weights.
   */
  // Weights weights() {
  //     return _weights;
  // }

  /**
   * Returns the number of reference samples. Might be zero.
   */
  int nRefSamples() const { return _nRefSamples; }
  /**
   * Returns the list of reference samples.
   */
  Samples refSamples() const { return _refSamples; }
  /**
   * Returns the number of target samples.
   */
  int nTargetSamples() const { return _targetSamples.nSamples(); }
  /**
   * Returns the list of target samples.
   */
  Samples targetSamples() const { return _targetSamples; }
  /**
   * Returns the number of reference and target samples.
   */
  int nAllSamples() const { return _allSamples.nSamples(); }
  /**
   * Returns a list of all target and reference samples.
   * Target samples are listed first in the same order as the list returned
   * by {@code this.targetSamples()}. Reference samples are listed last
   * in the same order as the list returned by {@code this.refSamples()}.
   */
  Samples allSamples() const { return _allSamples; }
  /**
   * Returns the number of target data markers.
   */
  int nTargetMarkers() const { return _targetMarkers.nMarkers(); }
  /**
   * Returns the list of target data markers.
   */
  Markers targetMarkers() const { return _targetMarkers; }
  /**
   * Returns the number of reference data markers.
   */
  int nMarkers() const { return _markers.nMarkers(); }
  /**
   * Returns the list of reference data markers.
   */
  Markers markers() const { return _markers; }
  /**
   * Returns the index of the specified marker in the reference data markers.
   * @param targetMarker index of a marker in the list of target data markers
   */
  int markerIndex(int targetMarker) const { return _markerIndex[targetMarker]; }
  /**
   * Returns an array of length {@code this.nTargetMarkers()} which maps
   * the {@code k}-th marker in the list of target data markers to the
   * index of the marker in the list of reference data markers.
   */
  QList<int> markerIndices() const { return _markerIndex; }
  /**
   * Returns the index of the specified marker in the target data, or
   * returns -1 if the marker is not present in the target data.
   * @param marker index of a marker in the reference data
   */
  // int targetMarkerIndex(int marker) const { return _targetMarkerIndex[marker]; }
  /**
   * Returns an array of length {@code this.nMarkers()} whose {@code k}-th
   * element is the index of the {@code k}-th marker in the list of target
   * markers or is -1 if the marker is not present in the target data.
   */
  // QList<int> targetMarkerIndices() const { return _targetMarkerIndex; }
  /**
   * Add the reference haplotype pairs that are restricted
   * to the target data markers to the specified list.
   * @param list a list of haplotype pairs for target data markers
   */
  void addRestrictedRefHapPairs(QList<HapPair> &list) const { list.append(_restRefHapPairs); }
  /**
   * Returns a list of reference haplotype pairs that are restricted
   * to the target data markers, or returns {@code null}
   * if there are no reference samples.
   */
  const SampleHapPairs &restrictedRefSampleHapPairs() const { return _restrictedRefSampleHapPairs; }
  /**
   * Returns a list of reference haplotype pairs, or returns {@code null}
   * if there are no reference samples.
   */
  const RefHapPairs &refSampleHapPairs() const { return _refSampleHapPairs; }
  /**
   * Returns the genotype likelihoods for the
   * target samples at the target data markers.
   */
  SplicedGL targetGL() const { return _targetGL; }
  /**
   * Returns an array whose initial element is {@code 0} and whose
   * {@code j}-th element for {@code j > 0} is the recombination rate
   * between the target markers with indices {@code (j - 1)} and {@code j}.
   *
   * @return inter-marker recombination rates for the target markers
   */
  // QList<float> recombRate() {
  //     return _recombRate==null ? null : _recombRate;
  // }

private:
  // QList<float> recombRate(Markers markers, GeneticMap map,
  //   float mapScale);

  int _window;
  SampleHapPairs _initHaps;
  int _prevSpliceStart;
  int _nextOverlapStart;
  int _nextSpliceStart;
  int _nextTargetSpliceStart;
  int _nextTargetOverlapStart;

  SplicedGL _targetGL;
  // NuclearFamilies families;
  // Weights weights;

  int _nRefSamples;
  Samples _refSamples;
  Samples _targetSamples;
  Samples _allSamples;

  Markers _markers;
  Markers _targetMarkers;
  // QList<int> _targetMarkerIndex;
  QList<int> _markerIndex;

  QList<HapPair> _restRefHapPairs;
  RefHapPairs _refSampleHapPairs;
  SampleHapPairs _restrictedRefSampleHapPairs;

  // QList<float> _recombRate;
};

class R2Estimator
{
public:
  void clear();
  void addSampleData(QVector<double> doseProbs);

  /**
   * Returns the estimated squared correlation between the most probable
   * ALT allele dose and the true ALT allele dose for the current
   * genotype data.
   * Returns 0 if the marker is monomorphic or if the most probable ALT
   * allele dose is monomorphic.
   */
  double allelicR2();

  /**
   * Returns the estimated squared correlation between the expected
   * ALT allele dose and the true ALT allele dose for the current
   * genotype data. Returns 0 if the marker is monomorphic.
   */
  double doseR2();

  /**
   * Returns the current number of genotypes with allele dose data.
   */
  int nGenotypes() { return _nGenotypes; }

private:
  int _nGenotypes;
  double _sumCall;
  double _sumSquareCall;
  double _sumExpected;
  double _sumExpectedSquare;
  double _sumSquareExpected;
  double _sumCallExpected;
};

class ConstrainedAlleleProbs;
class ConstrainedGLAlleleProbs;

class ImputeDataWriter
{
public:
  ImputeDataWriter(const Samples &samples) : _samples(samples) {}

  void printWindowOutput(const CurrentData &cd,
                         const SampleHapPairs &targetHapPairs,
                         const ConstrainedAlleleProbs &alProbs,
                         const Par &par);

  Samples samples() const { return _samples; }
  // Markers markers() const { return _markers; }

  virtual void writeHeader() = 0;
  virtual void writeEOF() = 0;

protected:
  virtual void initializeWindowBuffering(const int initSize) = 0;
  virtual void appendPhasedVariantData() = 0;
  virtual void finishAndWriteRec() = 0;
  /////////////////// virtual void debugWrite() = 0; ///////////////////////

  Samples _samples;
  // Markers _markers;

  QVector<bool> _isImputed;
  int _start;
  int _end;
  bool _printDS;
  bool _printGP;

  QVector<double> _gt3Probs;

  R2Estimator _r2Est;

  int _mNum;
  Marker _marker;
  int _nAlleles;
  int _nGenotypes;
  int _allele1;
  int _allele2;
  QVector<double> _gtProbs;
  QVector<double> _dose;
  QVector<double> _cumAlleleProbs;

  QVector<double> _alProbs1;
  QVector<double> _alProbs2;

private:
  void setIsImputed(const CurrentData &cd);
  void printWindowData(const ConstrainedAlleleProbs &alProbs);
  void initializeForWindow(const int initSize);
  void resetRec(const Marker &marker);
  void constructSampleDataForMarker();
  int maxIndex(QVector<double> &da, int expLength);
};

#endif
