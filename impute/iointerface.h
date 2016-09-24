/* Copyright 2016 Golden Helix, Inc. */
#ifndef IOUTILITIES_H
#define IOUTILITIES_H

#include "impute/haplotypepair.h"
#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"

class Par
{
public:
  virtual int window() const { return 50000; }
  virtual int overlap() const { return 3000; }
};

class GenericDataReader
{
public:
  virtual bool canAdvanceWindow() const = 0;
  virtual void makeNewWindow(int overlap) = 0;
  virtual void addNewDataToNewWindow(int windowSize) = 0;
  virtual int windowSize() const = 0;
  virtual int maxWindowSize() const = 0;
  virtual bool lastWindowOnChrom() const = 0;
};

class RefDataReader : public GenericDataReader
{
public:
  virtual bool canAdvanceWindow() const { return false; }
  void makeNewWindow(int overlap);
  void addNewDataToNewWindow(int windowSize);
  int windowSize() const { return _vcfRefRecs.length(); }
  int maxWindowSize() const { return 15000; }
  bool lastWindowOnChrom() const;

  Samples samples() const { return _samples; }
  Markers markers() const { return _markers; }
  QList<BitSetRefGT> refRecs() const { return _vcfRefRecs; }
protected:
  virtual bool hasNextRec() const { return false; }
  virtual BitSetRefGT nextRec() const { return BitSetRefGT(); }
  virtual void advanceRec() {}
  bool sameChrom(BitSetRefGT a, BitSetRefGT b) const;
  bool samePosition(BitSetRefGT a, BitSetRefGT b) const;
  int currentChromIndex() const;

  Samples _samples;
  Markers _markers;
  QList<BitSetRefGT> _oldVcfRefRecs;
  QList<BitSetRefGT> _vcfRefRecs;
};

class TargDataReader : public GenericDataReader
{
public:
  virtual bool canAdvanceWindow() const = 0;

  void makeNewWindow(int overlap);
  void addNewDataToNewWindow(int windowSize);
  int windowSize() const { return _vcfEmissions.length(); }
  int maxWindowSize() const { return 10000; }
  bool lastWindowOnChrom() const;

  Samples samples() const { return _samples; }
  Markers markers() const { return _markers; }
  QList<BitSetGT> vcfRecs() const { return _vcfEmissions; }
protected:
  virtual bool hasNextRec() const = 0;
  virtual BitSetGT nextRec() const = 0;
  virtual void advanceRec() = 0;

  bool sameChrom(BitSetGT a, BitSetGT b) const;
  bool samePosition(BitSetGT a, BitSetGT b) const;
  int currentChromIndex() const;

  Samples _samples;
  Markers _markers;
  QList<BitSetGT> _oldVcfEmissions;
  QList<BitSetGT> _vcfEmissions;
};

class VcfWindow
{
public:
  VcfWindow() : _overlap(0), _cumMarkerCnt(0) {}
  /**
   * Advances the sliding window of VCF records. The advanced window
   * (a list of {@code BitSetRefGT} or of {@code BitSetGT} objects) is kept in the
   * DataReader object.  The size of the advanced window and the
   * number of markers of overlap between the marker window
   * immediately before method invocation and the marker window
   * immediately after method invocation may differ from the
   * requested values.  If the advanced window size or overlap is
   * less than the requested value, the actual value will be as
   * large as possible. If {@code this.lastWindowOnChrom() == true}
   * before method invocation, then there will be no overlap between
   * the advanced window and the previous window.
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

class RestrictedVcfWindow
{
public:
  RestrictedVcfWindow() {}
};

class CurrentData;

class InputData
{
public:
  virtual bool canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const = 0;
  virtual void advanceWindow(int overlap, int windowSize, TargDataReader &tr,
                             RefDataReader &rr) = 0;
  virtual void setCdData(CurrentData &cd, const Par &par, const SampleHapPairs &overlapHaps,
                         const TargDataReader &tr, const RefDataReader &rr) = 0;
};

class TargetData : public InputData
{
public:
  TargetData() : _window(0) {}
  bool canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const;
  void advanceWindow(int overlap, int windowSize, TargDataReader &tr, RefDataReader &rr);
  void setCdData(CurrentData &cd, const Par &par, const SampleHapPairs &overlapHaps,
                 const TargDataReader &tr, const RefDataReader &rr);

private:
  int nextOverlapStart(int targetOverlap, const TargDataReader &tr);

  VcfWindow _vcfWindow;
  int _window;
  // Markers _markers;
  SplicedGL _gl;
};

class AllData : public InputData
{
public:
  AllData() : _window(0) {}
  bool canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const;
  void advanceWindow(int overlap, int windowSize, TargDataReader &tr, RefDataReader &rr);
  void setCdData(CurrentData &cd, const Par &par, const SampleHapPairs &overlapHaps,
                 const TargDataReader &tr, const RefDataReader &rr);

private:
  VcfWindow _refWindow;
  RestrictedVcfWindow _targetWindow;

  int _window;
  QList<int> refIndices;
  QList<int> targetIndices;
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
  int window() { return _window; }
  /**
   * Returns the first marker index in the overlap between this
   * marker window and the next marker window, or
   * returns {@code this.nMarkers()} there is no overlap.
   */
  int nextOverlapStart() { return _nextOverlapStart; }
  /**
   * Returns the first target marker index in the overlap between this
   * marker window and the next marker window, or
   * returns {@code this.nMarkers()} if there is no overlap or if there are
   * no target markers in the overlap.
   */
  int nextTargetOverlapStart() { return _nextTargetOverlapStart; }
  /**
   * Returns the first marker index after the splice point with
   * the previous marker window. Returns 0 if the current marker window
   * is the first marker window.
   */
  int prevSpliceStart() { return _prevSpliceStart; }
  /**
   * Returns the first marker index after the splice point between this
   * marker window and the next marker window, or returns
   * {@code this.nMarkers()} if there is no overlap or if there are
   * no markers after the splice point.
   */
  int nextSpliceStart() { return _nextSpliceStart; }
  /**
   * Returns the first target marker index after the splice point with
   * the previous marker window. Returns 0 if the current marker window
   * is the first marker window.
   */
  int prevTargetSpliceStart() { return (_initHaps.nHaps() == 0) ? 0 : _initHaps.nMarkers(); }
  /**
   * Returns the first target marker index after the splice point between this
   * marker window and the next marker window, or returns
   * {@code this.nTargetMarkers()} if there is no overlap or if there are
   * no target markers after the splice point
   */
  int nextTargetSpliceStart() { return _nextTargetSpliceStart; }
  /**
   * Returns the target data haplotype pairs in the segment of the current
   * marker window preceding the splice point with the previous marker window:
   * {@code this.targetMarkers().restrict(0, this.prevTargetSplice())}
   */
  SampleHapPairs initHaps() { return _initHaps; }
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
  int nRefSamples() { return _nRefSamples; }
  /**
   * Returns the list of reference samples.
   */
  Samples refSamples() { return _refSamples; }
  /**
   * Returns the number of target samples.
   */
  int nTargetSamples() { return _targetSamples.nSamples(); }
  /**
   * Returns the list of target samples.
   */
  Samples targetSamples() { return _targetSamples; }
  /**
   * Returns the number of reference and target samples.
   */
  int nAllSamples() { return _allSamples.nSamples(); }
  /**
   * Returns a list of all target and reference samples.
   * Target samples are listed first in the same order as the list returned
   * by {@code this.targetSamples()}. Reference samples are listed last
   * in the same order as the list returned by {@code this.refSamples()}.
   */
  Samples allSamples() { return _allSamples; }
  /**
   * Returns the number of target data markers.
   */
  int nTargetMarkers() { return _targetMarkers.nMarkers(); }
  /**
   * Returns the list of target data markers.
   */
  Markers targetMarkers() { return _targetMarkers; }
  /**
   * Returns the number of reference data markers.
   */
  int nMarkers() { return _markers.nMarkers(); }
  /**
   * Returns the list of reference data markers.
   */
  Markers markers() { return _markers; }
  /**
   * Returns the index of the specified marker in the reference data markers.
   * @param targetMarker index of a marker in the list of target data markers
   */
  int markerIndex(int targetMarker) { return _markerIndex[targetMarker]; }
  /**
   * Returns an array of length {@code this.nTargetMarkers()} which maps
   * the {@code k}-th marker in the list of target data markers to the
   * index of the marker in the list of reference data markers.
   */
  QList<int> markerIndices() { return _markerIndex; }
  /**
   * Returns the index of the specified marker in the target data, or
   * returns -1 if the marker is not present in the target data.
   * @param marker index of a marker in the reference data
   */
  int targetMarkerIndex(int marker) { return _targetMarkerIndex[marker]; }
  /**
   * Returns an array of length {@code this.nMarkers()} whose {@code k}-th
   * element is the index of the {@code k}-th marker in the list of target
   * markers or is -1 if the marker is not present in the target data.
   */
  QList<int> targetMarkerIndices() { return _targetMarkerIndex; }
  /**
   * Add the reference haplotype pairs that are restricted
   * to the target data markers to the specified list.
   * @param list a list of haplotype pairs for target data markers
   */
  void addRestrictedRefHapPairs(QList<HapPair> list) { list.append(_restRefHapPairs); }
  /**
   * Returns a list of reference haplotype pairs that are restricted
   * to the target data markers, or returns {@code null}
   * if there are no reference samples.
   */
  SampleHapPairs restrictedRefSampleHapPairs() { return _restrictedRefSampleHapPairs; }
  /**
   * Returns a list of reference haplotype pairs, or returns {@code null}
   * if there are no reference samples.
   */
  SampleHapPairs refSampleHapPairs() { return _refSampleHapPairs; }
  /**
   * Returns the genotype likelihoods for the
   * target samples at the target data markers.
   */
  SplicedGL targetGL() { return _targetGL; }
  /**
   * Returns an array whose initial element is {@code 0} and whose
   * {@code j}-th element for {@code j > 0} is the recombination rate
   * between the target markers with indices {@code (j - 1)} and {@code j}.
   *
   * @return inter-marker recombination rates for the target markers
   */
  // QList<double> recombRate() {
  //     return _recombRate==null ? null : _recombRate;
  // }

private:
  /* Returns the index of the first marker in the overlap */
  int nextOverlapStart(int targetOverlap, const TargDataReader &tr);

  // QList<double> recombRate(Markers markers, GeneticMap map,
  //   double mapScale);

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
  QList<int> _targetMarkerIndex;
  QList<int> _markerIndex;

  QList<HapPair> _restRefHapPairs;
  SampleHapPairs _refSampleHapPairs;
  SampleHapPairs _restrictedRefSampleHapPairs;

  QList<float> _recombRate;
};

#endif
