#ifndef WORKERRUNNER_H
#define WORKERRUNNER_H

#include "impute/imputedriver.h"
#include "impute/imputationdata.h"

#include <QHash>
#include <QObject>
#include <QList>

class ComputeWorker;

class ComputeRunner : public QObject
{
  Q_OBJECT;

public:
  void start(int sampleCount);
  bool hasNextWorkerResult();
  int nextResultIdx() const;
  void sendWorkerNextSample(int workerId);

  virtual bool hasResultsForSample(int sampleIdx) const = 0;

signals:
  void computeSample(int sampleIdx, int workerId);
  void shutdownWorkers();
  void nextResultArrived();

protected:
  void setupNewWorker(ComputeWorker* worker, QThread* thread);

  int _sampleCount;
  int _sentSampleIdx;
  int _retrievedSampleIdx;
  QList<ComputeWorker*> _workers;
};

class ComputeWorker : public QObject
{
  Q_OBJECT;

public:
  void setWorkerId(int workerId) { _workerId = workerId; }

public slots:
  virtual void computeSample(int single, int workerId) = 0;
  void throwThyselfOnSword();

protected:
  int _workerId;
};


class SingleBaumWorker;
class SingleBaumRunner : public ComputeRunner
{
  Q_OBJECT;

public:
  SingleBaumRunner();
  void moveWorkerToThread(SingleBaumWorker* worker, QThread* thread);
  QList<HapPair> nextWorkerResult();
  bool hasResultsForSample(int sampleIdx) const;

public slots:
  void resultsReady(int sampleIdx, int workerId, QList<HapPair> result);

private:
  QHash< int, QList<HapPair> > _results; // to be retrieved results
};

class SingleBaumWorker : public ComputeWorker
{
  Q_OBJECT;

public:
  SingleBaumWorker(Dag& dag, SplicedGL& gl, int seed, int nSamplingsPerIndividual, bool lowmem);

public slots:
  void computeSample(int single, int workerId);

signals:
  void resultsReady(int sampleIdx, int workerId, QList<HapPair>);

private:
  SingleBaum _baum;
};

class LSImputeWorker;
class LSImputeRunner : public ComputeRunner
{
  Q_OBJECT;

public:
  LSImputeRunner();
  void moveWorkerToThread(LSImputeWorker* worker, QThread* thread);
  QList<HapAlleleProbs> nextWorkerResult();
  bool hasResultsForSample(int sampleIdx) const;

public slots:
  void resultsReady(int sampleIdx, int workerId, QList<HapAlleleProbs> result);

private:
  QHash< int, QList<HapAlleleProbs> > _results; // to be retrieved results
};

class LSImputeWorker : public ComputeWorker
{
  Q_OBJECT;

public:
  LSImputeWorker(const Par &par, const CurrentData &cd,
                 const SampleHapPairs &targetHapPairs, double scaleFactor);

public slots:
  void computeSample(int sample, int workerId);

signals:
  void resultsReady(int sampleIdx, int workerId, QList<HapAlleleProbs>);

private:
  PositionMap _imputationMap;
  ImputationData _impData;
  LSHapBaum _hb;
};

#endif
