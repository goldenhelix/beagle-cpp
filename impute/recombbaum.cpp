#include "impute/recombbaum.h"

#include "impute/samplerdata.h"
#include "impute/singlepermittedstates.h"

#include <QMap>

#include <float.h> 
#include <math.h>

#define MIN_VALUE_FOR_BAUM      100.0 * FLT_MIN

#define LOAD_FACTOR             0.75

static qint64 hash1(int node1, int node2)
{
  qint64 hash = 5;
  hash = 71 * hash + node1;
  hash = 71 * hash + node2;
  return hash;
}

static qint64 hash2(int node1, int node2)
{
  qint64 hash = 7;
  hash = 97 * hash + node1;
  hash = 97 * hash + node2;
  return hash;
}

RecombSingleNodes::RecombSingleNodes() : _size(0)
{
  _capacity = (1<<10);
  _rehashThreshold = (int) (LOAD_FACTOR * _capacity);
  _index.fill(0, _capacity);
  _node1.fill(0, _capacity);
  _node2.fill(0, _capacity);
  _value.fill(0, _capacity);
}

void RecombSingleNodes::initialize(int nNodes)
{
  Q_ASSERT_X(nNodes >= 1,
             "RecombSingleNodes::initialize",
             "nNodes < 1");

  _nNodes = nNodes;
  _sumNode1Value.fill(0, nNodes);
  _sumNode2Value.fill(0, nNodes);
  _sumValue = 0.0f;
}

int RecombSingleNodes::index(int node1, int node2) const
{
  qint64 h1 = hash1(node1, node2);
  qint64 h2 = hash2(node1, node2);
  if ((h2 & 1)==0) {
    // h2 must be relatively prime to maxSize, which is a power of 2
    ++h2;
  }
  qint64 l = h1;
  for (int k=0; k<_capacity; ++k) {
    int i = (int) (l % _capacity);

    if (_value[i]==0.0 || (_node1[i]==node1 && _node2[i]==node2))
      return i;

    l += h2;
  }
  Q_ASSERT_X(false,
             "SingleNodes::index",
             "Can't find hash index!");
  return -1;
}

void RecombSingleNodes::rehash()
{
  Q_ASSERT_X(_size>=_rehashThreshold,
             "SingleNodes::rehash",
             "_size < _rehashThreshold");

  int newMaxSize = 2*_capacity;

  Q_ASSERT_X(newMaxSize >= 0,
             "SingleNodes::rehash",
             "newMaxSize < 0--hash table overflow");

  QVector<int> oldIndex = _index;
  QVector<int> oldNode1 = _node1;
  QVector<int> oldNode2 = _node2;
  QVector<float> oldValue = _value;

  _capacity = newMaxSize;
  _index.fill(0, newMaxSize);
  _node1.fill(0, newMaxSize);
  _node2.fill(0, newMaxSize);
  _value.fill(0, newMaxSize);

  for (int j=0; j<_size; ++j) {
    int oldInd = oldIndex[j];
    int newIndex = index(oldNode1[oldInd], oldNode2[oldInd]);
    _index[j] = newIndex;
    _node1[newIndex] = oldNode1[oldInd];
    _node2[newIndex] = oldNode2[oldInd];
    _value[newIndex] = oldValue[oldInd];
  }
  _rehashThreshold = (int) (LOAD_FACTOR * _capacity);
}

void RecombSingleNodes::sumUpdate(int node1, int node2, float value)
{
  Q_ASSERT_X(value > 0.0, "SingleNodes::sumUpdate", "value <= 0.0");

  int i = index(node1, node2);
  bool addNode = (_value[i]==0.0);
  _value[i] += value;
  _sumNode1Value[node1] += value;
  _sumNode2Value[node2] += value;
  _sumValue += value;
  if (addNode) {
    _index[_size++] = i;
    _node1[i] = node1;
    _node2[i] = node2;
    if (_size>=_rehashThreshold) {
      rehash();
    }
  }
}

void RecombSingleNodes::checkSize(int index) const
{
  Q_ASSERT_X(index < _size,
             "SingleNodes::checkSize",
             "index>=_size");
}

float RecombSingleNodes::value(int node1, int node2) const
{
  Q_ASSERT_X(node1 >= 0  &&  node1 < _nNodes,
             "SingleNodes::value",
             "node1 < 0 || node1 >= _nNodes");
  Q_ASSERT_X(node2 >= 0  &&  node2 < _nNodes,
             "SingleNodes::value",
             "node2 < 0 || node2 >= _nNodes");

  return _value[index(node1, node2)];
}

void RecombSingleNodes::clear()
{
  for (int j=0; j<_size; ++j) {
    _value[_index[j]] = 0.0f;
    _sumNode1Value[_node1[_index[j]]] = 0.0f;
    _sumNode2Value[_node2[_index[j]]] = 0.0f;
  }
  _sumValue = 0.0f;
  _size = 0;
}


RecombSingleBaumLevel::RecombSingleBaumLevel(const SamplerData &samplerData)
  : _marker(-1), _sample(-1), _size(0), _samplerData(samplerData),
    _dag(samplerData.rdag().dag()), _gl(samplerData.gl()),
    _fwdValueSum(0.0), _bwdValueSum(0.0) {}

void RecombSingleBaumLevel::setForwardValues(RecombSingleNodes &nodes,
                                             SinglePermittedStates &permittedStates,
                                             int marker, int sample)
{
  _marker = marker;
  _sample = sample;
  // _nGenotypes = gl.marker(marker).nGenotypes();
  _size = 0;
  _fwdValueSum = 0.0f;
  _bwdValueSum = 0.0f;
  // initializeGtProbs(); // initialized here due to gtProbs() contract
  setStates(nodes, permittedStates);
  setChildNodes(nodes);
}

void RecombSingleBaumLevel::reset()
{
  _size = 0;
  _edges1.clear();
  _edges2.clear();
  _fwdValues.clear();
}

void RecombSingleBaumLevel::setStates(const RecombSingleNodes &nodes,
                                      SinglePermittedStates &permittedStates)
{
  float valueSum = 0.0f;
  reset();
  permittedStates.setMarker(_marker);
  while (permittedStates.hasNext()) {
    permittedStates.next();
    int edge1 = permittedStates.edge1();
    int edge2 = permittedStates.edge2();
    float thisFwdValue = fwdValue(edge1, edge2, nodes);
    if (thisFwdValue > 0.0) {
      _edges1.append(edge1);
      _edges2.append(edge2);
      _fwdValues.append(thisFwdValue);
      valueSum += thisFwdValue;
    }
  }
  _size = _fwdValues.length();

  Q_ASSERT_X(valueSum > 0.0f,
             "RecombSingleBaumLevel::setStates",
             "valueSum <= 0f");

  for (int k=0; k < _size; ++k)
    _fwdValues[k] /= valueSum;

  _fwdValueSum = valueSum;
}

float RecombSingleBaumLevel::fwdValue(int edge1, int edge2, const RecombSingleNodes &nodes)
{
  float fwdValue = 0.0f;
  int symbol1 = _dag.symbol(_marker, edge1);
  int symbol2 = _dag.symbol(_marker, edge2);
  float emProb = _gl.gl(_marker, _sample, symbol1, symbol2);

  if (emProb > 0.0) {
    float pRecom = _samplerData.pRecomb(_marker);
    float rec0 = (1-pRecom)*(1-pRecom);
    float rec1 = pRecom*(1-pRecom);
    float rec2 = pRecom*pRecom;
    float ep1 = _dag.condEdgeProb(_marker, edge1);
    float ep2 = _dag.condEdgeProb(_marker, edge2);
    float ep = ep1*ep2;

    int pn1 = _dag.parentNode(_marker, edge1);
    int pn2 = _dag.parentNode(_marker, edge2);
    float pnp1 = _dag.parentProb(_marker, pn1);
    float pnp2 = _dag.parentProb(_marker, pn2);

    fwdValue = rec0*emProb*ep*nodes.value(pn1, pn2);
    fwdValue += rec1*emProb*ep*pnp2*nodes.sumNode1Value(pn1);
    fwdValue += rec1*emProb*ep*pnp1*nodes.sumNode2Value(pn2);
    fwdValue += rec2*emProb*ep*pnp1*pnp2*nodes.sumValue();
    if (fwdValue < MIN_VALUE_FOR_BAUM) {
      fwdValue = MIN_VALUE_FOR_BAUM;
    }
  }
  return fwdValue;
}

void RecombSingleBaumLevel::setChildNodes(RecombSingleNodes &nodes)
{
  nodes.clear();
  for (int k=0; k<_size; ++k) {
    int node1 = _dag.childNode(_marker, _edges1[k]);
    int node2 = _dag.childNode(_marker, _edges2[k]);
    nodes.sumUpdate(node1, node2, _fwdValues[k]);
  }
}

void RecombSingleBaumLevel::checkIndex(int state) const
{
  Q_ASSERT_X(state < _size, "SingleBaumLevel::checkIndex", "state >= _size");
}


RecombSingleBaum::RecombSingleBaum(const SamplerData &samplerData, int seed,
                                   int nSamplingsPerIndividual, bool lowMem)
  : _samplerData(samplerData), _rdag(samplerData.rdag()), _dag(samplerData.rdag().dag()),
  _gl(samplerData.gl()), _windowIndex(-9999), _arrayIndex(-9999),
  _seed(seed), _nSamplingsPerIndividual(nSamplingsPerIndividual)
{
  Q_ASSERT_X(nSamplingsPerIndividual >= 1,
             "RecombSingleBaum::RecombSingleBaum",
             "nSamplingsPerIndividual < 1");

  _nMarkers = samplerData.nMarkers();
  // _random = new Random(seed);

  QList<int> zeroList;
  for(int j=0; j < _nMarkers; j++)
    zeroList.append(0);

  for(int i=0; i < nSamplingsPerIndividual; i++)
  {
    _node1.append(0);
    _node2.append(0);
    _baseTrProb.append(0.0);
    _maxSum.append(0.0);
    _alleles1.append(zeroList);
    _alleles2.append(zeroList);
  }

  int size = _dag.nLevels();
  if (lowMem)
    size = (int) (sqrt(1.0 + 8.0*_dag.nLevels())/2.0) + 2;

  for (int j=0; j < size; ++j)
    _levels.append(RecombSingleBaumLevel(samplerData));

  _fwdNodes.initialize(_dag.maxNodes());
}

QList<HapPair> RecombSingleBaum::randomSample(int sample)
{
  SinglePermittedStates permittedStates(_rdag, sample);
  forwardAlgorithm(sample, permittedStates);
  initSampleAlleles(currentLevel(), sample);
  for (int j=_nMarkers-2; j >= 0; --j)
  {
    RecombSingleBaumLevel &level = previousLevel(sample, permittedStates);
    sampleAlleles(level, sample);
  }
  pruneLevels();
  return hapList(sample);
}

void RecombSingleBaum::pruneLevels()
{
  // Might as well clean out all the QLists now....
  for (int j=0; j < _levels.length(); ++j)
    _levels[j].reset();
}

/***
void RecombSingleBaum::pruneLevels()
{
  int meanSize = estMeanSize();
  int capacityThreshold = 3*meanSize;
  int newCapacity = 3*meanSize/2 + 1;
  for (int j=0; j < _levels.length(); ++j) {
    if (_levels[j].capacity() > capacityThreshold) {
      _levels[j].reset(newCapacity);
    }
  }
}

int RecombSingleBaum::estMeanSize()
{
  int nLevelsToSample = 20;
  int sizeSum = 0;
  for (int j=0; j<nLevelsToSample; ++j) {
    /// sizeSum += levels[random.nextInt(levels.length())].size();   /// %%%%%%%%
    sizeSum += levels[j % _levels.length()].size();                  /// %%%%%%%%
  }
  return (int) (sizeSum / nLevelsToSample);
}
***/

QList<HapPair> RecombSingleBaum::hapList(int sample) const
{
  QList<HapPair> hList;
  for (int copy=0; copy < _nSamplingsPerIndividual; ++copy) {
    HapPair hpair(_gl.markers(), _gl.samples(), sample,
                  _alleles1[copy], _alleles2[copy]);
    hList.append(hpair);
  }
  return hList;
}

void RecombSingleBaum::initSampleAlleles(const RecombSingleBaumLevel &level, int sample)
{
  for (int j=0; j < _nSamplingsPerIndividual; ++j)
    saveCurrentData(level, sample, j, initialRandomState(level, j));  /// %%%%%%%%
}

int RecombSingleBaum::initialRandomState(const RecombSingleBaumLevel &level, int copy)
{
  //////////////////  float d = random.nextFloat();
  //////////////////  float d = 0.5;
  float d = (float) copy / (float)(_nSamplingsPerIndividual + 1);
  float sum = 0.0f;
  for (int j=0, n=level.size(); j<n; ++j) {
    sum += level.forwardValue(j);
    if (d <= sum) {
      return j;
    }
  }
  // GR: Added floor of 0 here
  return qMax(level.size()-1, 0); // error in finite bit arithmetic encountered
}

void RecombSingleBaum::saveCurrentData(const RecombSingleBaumLevel &level, int sample,
                                       int copy, int stateIndex)
{
  int m = level.marker();
  int e1 = level.edge1(stateIndex);
  int e2 = level.edge2(stateIndex);
  int s1 = level.symbol1(stateIndex);
  int s2 = level.symbol2(stateIndex);
  _node1[copy] = level.parentNode1(stateIndex);
  _node2[copy] = level.parentNode2(stateIndex);
  float p1 = _dag.edgeProb(m, e1);
  float p2 = _dag.edgeProb(m, e2);
  _baseTrProb[copy] = p1*p2;

  _maxSum[copy] = level.forwardValue(stateIndex) * level.forwardValuesSum()
    / _gl.gl(m, sample, s1, s2);
  _alleles1[copy][m] = s1;
  _alleles2[copy][m] = s2;
}

void RecombSingleBaum::sampleAlleles(const RecombSingleBaumLevel &level, int sample)
{
  for (int j=0; j < _nSamplingsPerIndividual; ++j)
    saveCurrentData(level, sample, j, randomPreviousState(level, j));
}

int RecombSingleBaum::randomPreviousState(const RecombSingleBaumLevel &level, int copy)
{
  int m = level.marker();
  float np1 = _dag.parentProb(m+1, _node1[copy]);
  float np2 = _dag.parentProb(m+1, _node2[copy]);
  float pRecomb = _samplerData.pRecomb(m+1);
  /// float d = random.nextFloat() * _maxSum[copy];                                   /// %%%%%%%%
  float d = ((float) copy / (float)(_nSamplingsPerIndividual + 1)) * _maxSum[copy];   /// %%%%%%%%
  float sum = 0.0f;
  for (int j=0, n=level.size(); j<n; ++j) {
    float tp = 0.0f;
    bool noJump1 = level.childNode1(j)==_node1[copy];
    bool noJump2 = level.childNode2(j)==_node2[copy];
    if (noJump1 && noJump2) {
      tp += (1-pRecomb)*(1-pRecomb)*_baseTrProb[copy]/ (np1*np2);
    }
    if (noJump1) {
      tp += (1-pRecomb)*pRecomb*_baseTrProb[copy] / np1;
    }
    if (noJump2) {
      tp += pRecomb*(1-pRecomb)*_baseTrProb[copy] / np2;
    }
    tp += pRecomb*pRecomb*_baseTrProb[copy];

    sum += (level.forwardValue(j)*tp);
    if (d <= sum) {
      return j;
    }
  }
  return qMax(level.size()-1, 0); // if reached due to rounding
}

RecombSingleBaumLevel& RecombSingleBaum::nextLevel()
{
  _arrayIndex++;
  if (_arrayIndex == _levels.length())
  {
    _windowIndex++;
    _arrayIndex = _windowIndex;
  }
  return _levels[_arrayIndex];
}

RecombSingleBaumLevel &RecombSingleBaum::previousLevel(int sample,
                                                       SinglePermittedStates &permittedStates)
{
  if (_arrayIndex == _windowIndex) {
    _windowIndex--;
    _arrayIndex = _windowIndex;
    _levels[_arrayIndex].setChildNodes(_fwdNodes);
    int startLevel = _levels[_windowIndex].marker() + 1;
    int endLevel = startLevel + (_levels.length() - (_windowIndex + 1) );
    for (int marker=startLevel; marker<endLevel; marker++)
      nextLevel().setForwardValues(_fwdNodes, permittedStates, marker, sample);
    return currentLevel();
  }
  else
    return _levels[--_arrayIndex];
}

void RecombSingleBaum::forwardAlgorithm(int sample,
                                        SinglePermittedStates &permittedStates)
{
  _fwdNodes.clear();
  _fwdNodes.sumUpdate(0, 0, 1.0f);
  _windowIndex = -1;
  _arrayIndex = _levels.length() - 1;
  for (int marker = 0; marker < _nMarkers; marker++)
    nextLevel().setForwardValues(_fwdNodes, permittedStates, marker, sample);
}

