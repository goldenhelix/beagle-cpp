#include "impute/baumhmm.h"

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

SingleNodes::SingleNodes()
{
  _size = 0;
  _capacity = (1<<10);
  _rehashThreshold = (int) (LOAD_FACTOR * _capacity);
  _index.fill(0, _capacity);
  _node1.fill(0, _capacity);
  _node2.fill(0, _capacity);
  _value.fill(0, _capacity);
}

int SingleNodes::index(int node1, int node2) const
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

void SingleNodes::rehash()
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

void SingleNodes::sumUpdate(int node1, int node2, float value)
{
  Q_ASSERT_X(node1 >= 0, "SingleNodes::sumUpdate", "node1 < 0");
  Q_ASSERT_X(node2 >= 0, "SingleNodes::sumUpdate", "node2 < 0");
  Q_ASSERT_X(value > 0.0, "SingleNodes::sumUpdate", "value <= 0.0");

  int i = index(node1, node2);
  bool addNode = (_value[i]==0.0);
  _value[i] += value;
  if (addNode) {
    _index[_size++] = i;
    _node1[i] = node1;
    _node2[i] = node2;
    if (_size>=_rehashThreshold) {
      rehash();
    }
  }
}

void SingleNodes::checkSize(int index) const
{
  Q_ASSERT_X(index < _size,
             "SingleNodes::checkSize",
             "index>=_size");
}

float SingleNodes::value(int node1, int node2) const
{
  Q_ASSERT_X(node1 >= 0,
             "SingleNodes::value",
             "node1 < 0");
  Q_ASSERT_X(node2 >= 0,
             "SingleNodes::value",
             "node2 < 0");

  return _value[index(node1, node2)];
}

void SingleNodes::clear()
{
  for (int j=0; j<_size; ++j) {
    _value[_index[j]] = 0.0;
  }
  _size = 0;
}

SingleBaumLevel::SingleBaumLevel(const Dag &dag, const SplicedGL &gl)
  : _marker(-1), _sample(-1), _size(0), _dag(dag), _gl(gl),
    _fwdValueSum(0.0), _bwdValueSum(0.0)
{
  Q_ASSERT_X(dag.markers() == gl.markers(),
             "SingleBaumLevel::SingleBaumLevel",
             "inconsistent markers");
}

SingleBaumLevel& SingleBaumLevel::operator=(const SingleBaumLevel& other)
{
  //default assignment operator wont generate const casts
  const_cast<Dag&>(_dag) = other._dag;
  const_cast<SplicedGL&>(_gl) = other._gl;

  _marker = other._marker;
  _sample = other._sample;
  _size = other._size;

  _edges1 = other._edges1;
  _edges2 = other._edges2;
  _fwdValues = other._fwdValues;
  _bwdValues = other._bwdValues;
   
  _fwdValueSum = other._fwdValueSum;
  _bwdValueSum = other._bwdValueSum;
  return *this;
}

void SingleBaumLevel::setForwardValues(SingleNodes &nodes, int marker, int sample)
{
  _marker = marker;
  _sample = sample;
  // _nGenotypes = gl.marker(marker).nGenotypes();
  _size = 0;
  _fwdValueSum = 0.0;
  _bwdValueSum = 0.0;
  // initializeGtProbs(); // initialized here due to gtProbs() contract
  setStates(nodes);
  setChildNodes(nodes);
}

void SingleBaumLevel::setStates(const SingleNodes &nodes)
{
  float valueSum = 0.0;
  _edges1.clear();
  _edges2.clear();
  _fwdValues.clear();

  for (int j=0, n=nodes.size(); j<n; ++j) {
    int node1 = nodes.enumNode1(j);
    int node2 = nodes.enumNode2(j);

    for (int i1=0, nI1=_dag.nOutEdges(_marker, node1); i1<nI1; ++i1)
    {
      int edge1 = _dag.outEdge(_marker, node1, i1);
      int symbol1 = _dag.symbol(_marker, edge1);
      for (int i2=0, nI2=_dag.nOutEdges(_marker, node2); i2<nI2; ++i2)
      {
        int edge2 = _dag.outEdge(_marker, node2, i2);
        int symbol2 = _dag.symbol(_marker, edge2);
        float ep = _gl.gl(_marker, _sample, symbol1, symbol2);
        if (ep > 0.0)
        {
          _edges1.append(edge1);
          _edges2.append(edge2);
          float tp1 = _dag.condEdgeProb(_marker, edge1);
          float tp2 = _dag.condEdgeProb(_marker, edge2);

          float nodeValue = nodes.enumValue(j);

          float fwdValue = ep * nodeValue * (tp1 * tp2);

          if (fwdValue < MIN_VALUE_FOR_BAUM  &&  nodeValue > 0.0)
            fwdValue = MIN_VALUE_FOR_BAUM;

          _fwdValues.append(fwdValue);
          valueSum += fwdValue;
        }
      }
    }
  }
  _size = _fwdValues.length();
        
  Q_ASSERT_X(valueSum > 0.0 || _size==0, "SingleBaumLevel::setStates",
             "both valueSum and _size are zero");
  for (int k=0; k < _size; ++k) {
    _fwdValues[k] /= valueSum;
  }
  _fwdValueSum = valueSum;
}

void SingleBaumLevel::setChildNodes(SingleNodes &nodes)
{
  nodes.clear();
  for (int k=0; k<_size; ++k) {
    int node1 = _dag.childNode(_marker, _edges1[k]);
    int node2 = _dag.childNode(_marker, _edges2[k]);
    nodes.sumUpdate(node1, node2, _fwdValues[k]);
  }
}

void SingleBaumLevel::checkIndex(int state) const
{
  Q_ASSERT_X(state < _size, "SingleBaumLevel::checkIndex", "state >= _size");
}


SingleBaum::SingleBaum(const Dag &dag, const SplicedGL &gl, int seed, int nSamplingsPerIndividual,
                       bool lowMem) : _dag(dag), _gl(gl), _seed(seed), _windowIndex(-9999), _arrayIndex(-9999)
{
  Q_ASSERT_X(dag.markers() == gl.markers(),
             "SingleBaum::SingleBaum",
             "inconsistent markers");
  Q_ASSERT_X(nSamplingsPerIndividual >= 1,
             "SingleBaum::SingleBaum",
             "nSamplingsPerIndividual < 1");

  _nMarkers = dag.nLevels();
  _nSamplingsPerIndividual = nSamplingsPerIndividual;

  qsrand(seed);

  QList<int> zeroList;
  for(int j=0; j < _nMarkers; j++)
    zeroList.append(0);

  for(int i=0; i < nSamplingsPerIndividual; i++)
  {
    _node1.append(0);
    _node2.append(0);
    _nodeValue.append(0.0);
    _alleles1.append(zeroList);
    _alleles2.append(zeroList);
  }

  int size = dag.nLevels();
  if (lowMem)
    size = (int) (sqrt(1.0 + 8.0*dag.nLevels())/2.0) + 2;

  for (int j=0; j < size; ++j)
    _levels.append(SingleBaumLevel(dag, gl));
}

QList<HapPair> SingleBaum::randomSample(int sample)
{
  forwardAlgorithm(sample);
  initSampleAlleles(currentLevel(), sample);
  for (int j=_nMarkers-2; j >= 0; --j)
  {
    SingleBaumLevel &level = previousLevel(sample);
    sampleAlleles(level, sample);
  }
  return hapList(sample);
}

QList<HapPair> SingleBaum::hapList(int sample) const
{
  QList<HapPair> hList;
  for (int copy=0; copy < _nSamplingsPerIndividual; ++copy) {
    HapPair hpair(_gl.markers(), _gl.samples(), sample,
                  _alleles1[copy], _alleles2[copy]);
    hList.append(hpair);
  }
  return hList;
}

void SingleBaum::initSampleAlleles(const SingleBaumLevel &level, int sample)
{
  int m = level.marker();
  for (int copy=0; copy < _nSamplingsPerIndividual; ++copy)
  {
    int state = initialRandomState(level, copy);
    _node1[copy] = level.parentNode1(state);
    _node2[copy] = level.parentNode2(state);
    _nodeValue[copy] =  parentSum(level, sample, state);
    _alleles1[copy][m] = level.symbol1(state);
    _alleles2[copy][m] = level.symbol2(state);
  }
}

int SingleBaum::initialRandomState(const SingleBaumLevel &level, int copy)
{

#ifdef SIMULATE_RANDOM
  double d = (double) copy / (double)(_nSamplingsPerIndividual + 1);
#else
  double d = (double) qrand() / (double) RAND_MAX;
#endif

  double sum = 0.0;
  for (int j=0, n=level.size(); j<n; ++j) {
    sum += level.forwardValue(j);
    if (d <= sum) {
      return j;
    }
  }
  // GR: Added floor of 0 here
  return qMax(level.size()-1, 0); // if reached due to rounding
}

double SingleBaum::parentSum(const SingleBaumLevel &level, int sample, int state) const
{
  int marker = level.marker();
  double fwdValue = level.forwardValuesSum()*level.forwardValue(state);
  int edge1 = level.edge1(state);
  int edge2 = level.edge2(state);
  double tp1 = _dag.condEdgeProb(marker, edge1);
  double tp2 = _dag.condEdgeProb(marker, edge2);
  int symbol1 = _dag.symbol(marker, edge1);
  int symbol2 = _dag.symbol(marker, edge2);
  double ep = _gl.gl(marker, sample, symbol1, symbol2);
  return fwdValue / ( ep*tp1*tp2 );
}

void SingleBaum::sampleAlleles(const SingleBaumLevel &level, int sample)
{
  int m = level.marker();
  for (int copy=0; copy < _nSamplingsPerIndividual; ++copy)
  {
    int state = randomPreviousState(level, _node1[copy], _node2[copy],
                                    _nodeValue[copy], copy);
    _node1[copy] = level.parentNode1(state);
    _node2[copy] = level.parentNode2(state);
    _nodeValue[copy] =  parentSum(level, sample, state);
    _alleles1[copy][m] = level.symbol1(state);
    _alleles2[copy][m] = level.symbol2(state);
  }
}

int SingleBaum::randomPreviousState(const SingleBaumLevel &level, int node1,
                                    int node2, double nodeValue, int copy)
{

#ifdef SIMULATE_RANDOM
  double d = ((double)copy / (double)(_nSamplingsPerIndividual + 1)) * nodeValue;
#else
  double d = ((double) qrand() / (double) RAND_MAX) * nodeValue;
#endif

  double sum = 0.0;
  for (int j=0, n=level.size(); j<n; ++j) {
    if ( node1==level.childNode1(j)
         && node2==level.childNode2(j) ) {
      sum += level.forwardValue(j);
      if (d <= sum) {
        return j;
      }
    }
  }
  return level.size()-1; // error in finite bit arithmetic encountered
}

SingleBaumLevel& SingleBaum::nextLevel()
{
  _arrayIndex++;
  if (_arrayIndex == _levels.length())
  {
    _windowIndex++;
    _arrayIndex = _windowIndex;
  }
  return _levels[_arrayIndex];
}

SingleBaumLevel& SingleBaum::previousLevel(int sample)
{
  if (_arrayIndex == _windowIndex) {
    _windowIndex--;
    _arrayIndex = _windowIndex;
    _levels[_arrayIndex].setChildNodes(_fwdNodes);
    int startLevel = _levels[_windowIndex].marker() + 1;
    int endLevel = startLevel + (_levels.length() - (_windowIndex + 1) );
    for (int marker=startLevel; marker<endLevel; marker++)
      nextLevel().setForwardValues(_fwdNodes, marker, sample);
    return currentLevel();
  }
  else
    return _levels[--_arrayIndex];
}

void SingleBaum::forwardAlgorithm(int sample)
{
  _fwdNodes.clear();
  _fwdNodes.sumUpdate(0, 0, 1.0);
  _windowIndex = -1;
  _arrayIndex = _levels.length() - 1;
  for (int marker = 0; marker < _nMarkers; marker++)
    nextLevel().setForwardValues(_fwdNodes, marker, sample);
}

