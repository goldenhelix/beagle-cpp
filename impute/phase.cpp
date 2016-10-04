#include "impute/phase.h"

#include <QMap>

#include <math.h>

#define MIN_VALUE_FOR_BAUM   (100.*FLT_MIN)


void SingleNodes::sumUpdate(int node1, int node2, double value)
{
  Q_ASSERT_X(node1 >= 0, "SingleNodes::sumUpdate", "node1 < 0");
  Q_ASSERT_X(node2 >= 0, "SingleNodes::sumUpdate", "node2 < 0");
  Q_ASSERT_X(value > 0.0, "SingleNodes::sumUpdate", "value <= 0.0");

  IntPair p(node1, node2);
  _nodes.insert(p, _nodes.value(p, 0.0) + value);
}

QMapIterator<IntPair, double> SingleNodes::nodeIterator()
{
  QMapIterator<IntPair, double> i(_nodes);
  return i;
}

double SingleNodes::value(int node1, int node2) const
{
  Q_ASSERT_X(node1 >= 0,
	     "SingleNodes::value",
	     "node1 < 0");
  Q_ASSERT_X(node2 >= 0,
	     "SingleNodes::value",
	     "node2 < 0");

  IntPair p(node1, node2);
  return _nodes[p];
}

void SingleNodes::clear()
{
  _nodes.clear();
}


SingleBaumLevel::SingleBaumLevel(Dag *dag, SplicedGL *gl)
  : _marker(-1), _sample(-1), _size(0), _dag(dag), _gl(gl),
    _fwdValueSum(0.0), _bwdValueSum(0.0)
{
  Q_ASSERT_X(dag->markers() == gl->markers(),
             "SingleBaumLevel::SingleBaumLevel",
             "inconsistent markers");
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

void SingleBaumLevel::setStates(SingleNodes &nodes)
{
  double valueSum = 0.0;
  _edges1.clear();
  _edges2.clear();
  _fwdValues.clear();
  QMapIterator<IntPair, double> nit = nodes.nodeIterator();
  while (nit.hasNext())
  {
    nit.next();
    IntPair ip = nit.key();
    int node1 = ip.firstInt();
    int node2 = ip.secondInt();
    double nodeValue = nit.value();
    for (int i1=0, nI1=_dag->nOutEdges(_marker, node1); i1<nI1; ++i1)
    {
      int edge1 = _dag->outEdge(_marker, node1, i1);
      int symbol1 = _dag->symbol(_marker, edge1);
      for (int i2=0, nI2=_dag->nOutEdges(_marker, node2); i2<nI2; ++i2)
      {
        int edge2 = _dag->outEdge(_marker, node2, i2);
        int symbol2 = _dag->symbol(_marker, edge2);
        double ep = _gl->gl(_marker, _sample, symbol1, symbol2);
        if (ep > 0.0)
        {
          _edges1.append(edge1);
          _edges2.append(edge2);
          double tp1 = _dag->condEdgeProb(_marker, edge1);
          double tp2 = _dag->condEdgeProb(_marker, edge2);
          double fwdValue = ep * nodeValue * (tp1 * tp2);
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
    int node1 = _dag->childNode(_marker, _edges1[k]);
    int node2 = _dag->childNode(_marker, _edges2[k]);
    nodes.sumUpdate(node1, node2, _fwdValues[k]);
  }
}

void SingleBaumLevel::setBackwardValues(SingleNodes &nodes)
{
  _bwdValues.clear();
  for (int j=0; j < _size; ++j)
  {
    int node1 = _dag->childNode(_marker, _edges1[j]);
    int node2 = _dag->childNode(_marker, _edges2[j]);
    double backwardValue = nodes.value(node1, node2);
    _bwdValues.append(backwardValue);
    _bwdValueSum += backwardValue;
  }
  nodes.clear();

  // float gtProbsSum = 0f;

  for (int j=0; j < _size; ++j)
  {
    _bwdValues[j] /= _bwdValueSum;
    int edge1 = _edges1[j];
    int edge2 = _edges2[j];
    int symb1 = symbol1(j);
    int symb2 = symbol2(j);
    double tp1 = _dag->condEdgeProb(_marker, edge1);
    double tp2 = _dag->condEdgeProb(_marker, edge2);

    // float stateProb = fwdValues[j] * bwdValues[j];
    // int gtIndex = BasicGL.genotype(symb1, symb2);
    // // gtProbs assumed to be initialized in setForwardValues() method
    // gtProbs[gtIndex] += stateProb;
    // gtProbsSum += stateProb;

    double ep = _gl->gl(_marker, _sample, symb1, symb2);
    double bwdValue = _bwdValues[j] * (tp1 * tp2) * ep;
    if (bwdValue < MIN_VALUE_FOR_BAUM && _bwdValues[j] > 0.0)
      bwdValue = MIN_VALUE_FOR_BAUM;

    if (bwdValue > 0.0)
    {
      int pn1 = _dag->parentNode(_marker, edge1);
      int pn2 = _dag->parentNode(_marker, edge2);
      nodes.sumUpdate(pn1, pn2, bwdValue);
    }
  }

  // for (int j=0; j<nGenotypes; ++j) {
  //     gtProbs[j] /= gtProbsSum;
  // }
}

void SingleBaumLevel::checkIndex(int state) const
{
  Q_ASSERT_X(state < _size, "SingleBaumLevel::checkIndex", "state >= _size");
}


SingleBaum::SingleBaum(Dag &dag, SplicedGL &gl, int seed, int nSamplesPerIndividual,
                       bool lowMem) : _windowIndex(-9999), _arrayIndex(-9999)
{
  Q_ASSERT_X(dag.markers() == gl.markers(),
	     "SingleBaum::SingleBaum",
             "inconsistent markers");
  Q_ASSERT_X(nSamplesPerIndividual >= 1,
             "SingleBaum::SingleBaum",
             "nSamplesPerIndividual < 1");

  _dag = &dag;
  _gl = &gl;

  _nMarkers = dag.nLevels();
  _nSamplesPerIndividual = nSamplesPerIndividual;
  _seed = seed;
  // _random = new Random(seed);

  QList<int> zeroList;
  for(int j=0; j < _nMarkers; j++)
    zeroList.append(0);

  for(int i=0; i < nSamplesPerIndividual; i++)
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
    _levels.append(SingleBaumLevel(&dag, &gl));
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
  for (int copy=0; copy < _nSamplesPerIndividual; ++copy) {
    HapPair hpair(_gl->markers(), _gl->samples(), sample,
                  _alleles1[copy], _alleles2[copy]);
    hList.append(hpair);
  }
  return hList;
}

void SingleBaum::initSampleAlleles(const SingleBaumLevel &level, int sample)
{
  int m = level.marker();
  for (int copy=0; copy < _nSamplesPerIndividual; ++copy)
  {
    int state = initialRandomState(level);
    _node1[copy] = level.parentNode1(state);
    _node2[copy] = level.parentNode2(state);
    _nodeValue[copy] =  parentSum(level, sample, state);
    _alleles1[copy][m] = level.symbol1(state);
    _alleles2[copy][m] = level.symbol2(state);
  }
}

int SingleBaum::initialRandomState(const SingleBaumLevel &level)
{
  //////////////////  double d = random.nextDouble();
  double d = 0.5;
  double sum = 0.0;
  for (int j=0, n=level.size(); j<n; ++j) {
    sum += level.forwardValue(j);
    if (d <= sum) {
      return j;
    }
  }
  return level.size()-1; // if reached due to rounding
}

double SingleBaum::parentSum(const SingleBaumLevel &level, int sample, int state) const
{
  int marker = level.marker();
  double fwdValue = level.forwardValuesSum()*level.forwardValue(state);
  int edge1 = level.edge1(state);
  int edge2 = level.edge2(state);
  double tp1 = _dag->condEdgeProb(marker, edge1);
  double tp2 = _dag->condEdgeProb(marker, edge2);
  int symbol1 = _dag->symbol(marker, edge1);
  int symbol2 = _dag->symbol(marker, edge2);
  double ep = _gl->gl(marker, sample, symbol1, symbol2);
  return fwdValue / ( ep*tp1*tp2 );
}

void SingleBaum::sampleAlleles(const SingleBaumLevel &level, int sample)
{
  int m = level.marker();
  for (int copy=0; copy < _nSamplesPerIndividual; ++copy)
  {
    int state = randomPreviousState(level, _node1[copy], _node2[copy],
                                    _nodeValue[copy]);
    _node1[copy] = level.parentNode1(state);
    _node2[copy] = level.parentNode2(state);
    _nodeValue[copy] =  parentSum(level, sample, state);
    _alleles1[copy][m] = level.symbol1(state);
    _alleles2[copy][m] = level.symbol2(state);
  }
}

int SingleBaum::randomPreviousState(const SingleBaumLevel &level, int node1,
                                    int node2, double nodeValue)
{
  //////////////////  double d = random.nextDouble() * nodeValue;
  double d = 0.5 * nodeValue;
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
  for (int marker=0; marker < _nMarkers; marker++)
    nextLevel().setForwardValues(_fwdNodes, marker, sample);
}

