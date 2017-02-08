#include "impute/workerrunner.h"

#include <QEventLoop>
#include <QThread>

// ---------------
// ComputeRunner
// ----------------

void ComputeRunner::start(int sampleCount)
{
  _sampleCount = sampleCount;
  _sentSampleIdx = 0;
  _retrievedSampleIdx = 0;
  for (int workerId = 0; workerId < _workers.size(); workerId++) {
    sendWorkerNextSample(workerId);
  }
}

bool ComputeRunner::hasNextWorkerResult()
{
  if (_retrievedSampleIdx == _sampleCount) {
    emit shutdownWorkers();
    return false;
  }

  // If we already have the result for the next idx, return now
  if (hasResultsForSample(_retrievedSampleIdx))
    return true;

  // Otherwise wait for it
  QEventLoop loop;
  connect(this, SIGNAL(nextResultArrived()), &loop, SLOT(quit()));
  loop.exec();

  return true;
}

int ComputeRunner::nextResultIdx() const
{
  return _retrievedSampleIdx;
}

void ComputeRunner::sendWorkerNextSample(int workerId)
{
  if (_sentSampleIdx < _sampleCount) {
    int sampleIdx = _sentSampleIdx++;
    emit computeSample(sampleIdx, workerId);
  }
}

void ComputeRunner::setupNewWorker(ComputeWorker* worker, QThread* thread)
{
  int workerId = _workers.size();
  worker->setWorkerId(workerId);
  _workers << worker;
  worker->moveToThread(thread);

  // Because 'this' and 'thread' now have different thread affinities,
  // these connections are "queued" and cross thread boundaries safely
  // using serialized events pass using each threads event loop. We
  // won't receive any messages unless we are running an event loop in
  // this thread (which we run manually in hasNextWorkerResult)
  connect(this, SIGNAL(computeSample(int, int)), worker, SLOT(computeSample(int, int)));
  connect(this, SIGNAL(shutdownWorkers()), worker, SLOT(throwThyselfOnSword()));
}

// -------------
// ComputeWorker
// -------------
void ComputeWorker::throwThyselfOnSword()
{
  _workerId = -1;
  this->deleteLater();
}


// ---------------
// SingleBaumRunner
// ----------------

Q_DECLARE_METATYPE( QList<HapPair> );
SingleBaumRunner::SingleBaumRunner()
{
  qRegisterMetaType< QList<HapPair> >();
}

void SingleBaumRunner::moveWorkerToThread(SingleBaumWorker* worker, QThread* thread)
{
  ComputeRunner::setupNewWorker(worker, thread);

  // Our results signal
  connect(worker, SIGNAL(resultsReady(int, int, QList<HapPair>)), this,
          SLOT(resultsReady(int, int, QList<HapPair>)));
}

QList<HapPair> SingleBaumRunner::nextWorkerResult()
{
  int idx = _retrievedSampleIdx++;
  return _results.take(idx);
}

bool SingleBaumRunner::hasResultsForSample(int sampleIdx) const
{
  return _results.contains(sampleIdx);
}

void SingleBaumRunner::resultsReady(int sampleIdx, int workerId, QList<HapPair> result)
{
  _results[sampleIdx] = result;
  sendWorkerNextSample(workerId);
  if (sampleIdx == _retrievedSampleIdx)
    emit nextResultArrived();
}


// ---------------
// SingleBaumWorker
// ----------------

SingleBaumWorker::SingleBaumWorker(const Dag& dag, const SplicedGL& gl, int seed, int nSamplingsPerIndividual,
                                   bool lowmem)
  : _baum(dag, gl, seed, nSamplingsPerIndividual, lowmem)
{
}

void SingleBaumWorker::computeSample(int single, int workerId)
{
  if (workerId != _workerId)
    return;
  QList<HapPair> result = _baum.randomSample(single);
  emit resultsReady(single, workerId, result);
}

// ---------------
// RecombSingleBaumRunner
// ----------------

// NOTE: QList<HapPair> has already been declared as a metatype.
RecombSingleBaumRunner::RecombSingleBaumRunner()
{
  qRegisterMetaType< QList<HapPair> >();
}

void RecombSingleBaumRunner::moveWorkerToThread(RecombSingleBaumWorker* worker, QThread* thread)
{
  ComputeRunner::setupNewWorker(worker, thread);

  // Our results signal
  connect(worker, SIGNAL(resultsReady(int, int, QList<HapPair>)), this,
          SLOT(resultsReady(int, int, QList<HapPair>)));
}

QList<HapPair> RecombSingleBaumRunner::nextWorkerResult()
{
  int idx = _retrievedSampleIdx++;
  return _results.take(idx);
}

bool RecombSingleBaumRunner::hasResultsForSample(int sampleIdx) const
{
  return _results.contains(sampleIdx);
}

void RecombSingleBaumRunner::resultsReady(int sampleIdx, int workerId, QList<HapPair> result)
{
  _results[sampleIdx] = result;
  sendWorkerNextSample(workerId);
  if (sampleIdx == _retrievedSampleIdx)
    emit nextResultArrived();
}

// ---------------
// RecombSingleBaumWorker
// ----------------

RecombSingleBaumWorker::RecombSingleBaumWorker(const SamplerData &samplerData, int seed,
                                               int nSamplingsPerIndividual, bool lowmem)
  : _rbaum(samplerData, seed, nSamplingsPerIndividual, lowmem)
{
}

void RecombSingleBaumWorker::computeSample(int single, int workerId)
{
  if (workerId != _workerId)
    return;
  QList<HapPair> result = _rbaum.randomSample(single);
  emit resultsReady(single, workerId, result);
}

// ---------------
// LSImputeRunner
// ----------------

Q_DECLARE_METATYPE( QList<HapAlleleProbs> );
LSImputeRunner::LSImputeRunner()
{
  qRegisterMetaType< QList<HapAlleleProbs> >();
}

void LSImputeRunner::moveWorkerToThread(LSImputeWorker* worker, QThread* thread)
{
  ComputeRunner::setupNewWorker(worker, thread);

  // Our results signal
  connect(worker, SIGNAL(resultsReady(int, int, QList<HapAlleleProbs>)), this,
          SLOT(resultsReady(int, int, QList<HapAlleleProbs>)));
}

QList<HapAlleleProbs> LSImputeRunner::nextWorkerResult()
{
  int idx = _retrievedSampleIdx++;
  return _results.take(idx);
}

bool LSImputeRunner::hasResultsForSample(int sampleIdx) const
{
  return _results.contains(sampleIdx);
}

void LSImputeRunner::resultsReady(int sampleIdx, int workerId, QList<HapAlleleProbs> result)
{
  _results[sampleIdx] = result;
  sendWorkerNextSample(workerId);
  if (sampleIdx == _retrievedSampleIdx)
    emit nextResultArrived();
}


// ---------------
// LSImputeWorker
// ----------------

LSImputeWorker::LSImputeWorker(const Par &par, const CurrentData &cd,
                 const SampleHapPairs &targetHapPairs, double scaleFactor)
    : _imputationMap(scaleFactor), _impData(par, cd, targetHapPairs, _imputationMap), _hb(_impData, par.lowMem())
{
}

void LSImputeWorker::computeSample(int sample, int workerId)
{
  if (workerId != _workerId)
    return;
  QList<HapAlleleProbs> result;
  int hap1 = 2 * sample;
  int hap2 = 2 * sample + 1;

  result << _hb.randomHapSample(hap1);
  result << _hb.randomHapSample(hap2);
  emit resultsReady(sample, workerId, result);
}
