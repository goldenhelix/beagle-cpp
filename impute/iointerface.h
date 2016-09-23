/* Copyright 2016 Golden Helix, Inc. */
#ifndef IOUTILITIES_H
#define IOUTILITIES_H

#include "impute/haplotypepair.h"
#include "impute/markers.h"
#include "impute/samples.h"
#include "impute/vcfemission.h"

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
  QList<BitSetRefGT> refRecs() const { return _vcfRefRecs; }
protected:
  virtual bool hasNextRec() const { return false; }
  virtual BitSetRefGT nextRec() const { return BitSetRefGT(); }
  virtual void advanceRec() {}
  bool sameChrom(BitSetRefGT a, BitSetRefGT b) const;
  bool samePosition(BitSetRefGT a, BitSetRefGT b) const;
  int currentChromIndex() const;

  Samples _samples;
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
  QList<BitSetGT> vcfRecs() const { return _vcfEmissions; }
protected:
  virtual bool hasNextRec() const = 0;
  virtual BitSetGT nextRec() const = 0;
  virtual void advanceRec() = 0;

  bool sameChrom(BitSetGT a, BitSetGT b) const;
  bool samePosition(BitSetGT a, BitSetGT b) const;
  int currentChromIndex() const;

  Samples _samples;
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
  // virtual void setCdData(CurrentData &cd, const SampleHapPairs &overlapHaps,
  //		 const TargDataReader &tr, const RefDataReader &rr);
};

class TargetData : public InputData
{
public:
  TargetData() : _window(0) {}
  bool canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const;
  void advanceWindow(int overlap, int windowSize, TargDataReader &tr, RefDataReader &rr);
  // void setCdData(CurrentData &cd, const SampleHapPairs &overlapHaps,
  //	 const TargDataReader &tr, const RefDataReader &rr);

private:
  VcfWindow _vcfWindow;

  int _window;
  Markers _markers;
  SplicedGL _gl;
};

class AllData : public InputData
{
public:
  AllData() : _window(0) {}
  bool canAdvanceWindow(const TargDataReader &tr, const RefDataReader &rr) const;
  void advanceWindow(int overlap, int windowSize, TargDataReader &tr, RefDataReader &rr);
  // void setCdData(CurrentData &cd, const SampleHapPairs &overlapHaps,
  //	 const TargDataReader &tr, const RefDataReader &rr);

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
  friend class TargetData;  //// (If we go with private data, then we need "friends".)
  friend class AllData;

public:
  CurrentData();
};

#endif
