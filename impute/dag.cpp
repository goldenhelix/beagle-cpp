#include "impute/dag.h"

#include <math.h>
#include <stdio.h>

/* should be defined in math.h
static double log10(double x)
{
  return log(x) / log(10.);
}
*/

LinkageEquilibriumDag::LinkageEquilibriumDag(const SplicedGL &gl, double minFreq)
{
  Q_ASSERT_X(minFreq > 0.0 || minFreq < 0.5,
             "LinkageEquilibriumDag::LinkageEquilibriumDag",
             "bad value for minFreq");

  int nMarkers = gl.nMarkers();
  int localMaxAlleles = 0;
  _markers = gl.markers();

  for (int marker = 0; marker < nMarkers; ++marker) {
    _alleleFreq.append(alleleFrequencies(gl, marker, minFreq));
    if (_alleleFreq[marker].length() > localMaxAlleles) {
      localMaxAlleles = _alleleFreq[marker].length();
    }
  }
  _maxAlleles = localMaxAlleles;
  _sumAlleles = gl.markers().sumAlleles();
}

QList<double> LinkageEquilibriumDag::alleleFrequencies(const SplicedGL &gl,
                                                       int marker,
                                                       double minFreq)
{
  int nSamples = gl.nSamples();
  int nAlleles = gl.marker(marker).nAlleles();
  QList<double> alleleFreq;
  QList<double> scaledFreq;
  for (int a1 = 0; a1 < nAlleles; a1++) {
    alleleFreq.append(0.0);
    scaledFreq.append(0.0);
  }
  for (int sample = 0; sample < nSamples; ++sample) {
    for (int a1 = 0; a1 < nAlleles; ++a1) {
      for (int a2 = 0; a2 < nAlleles; ++a2) {
        double likelihood = gl.gl(marker, sample, a1, a2);
        scaledFreq[a1] += likelihood;
        scaledFreq[a2] += likelihood;
      }
    }
    divideEntriesBySum(scaledFreq);
    for (int j = 0; j < scaledFreq.length(); ++j) {
      alleleFreq[j] += scaledFreq[j];
      scaledFreq[j] = 0.0f;
    }
  }
  divideEntriesBySum(alleleFreq);
  enforceMinFrequency(alleleFreq, minFreq);
  return alleleFreq;
}

void LinkageEquilibriumDag::divideEntriesBySum(QList<double> &fa)
{
  double sum = 0.0f;
  for (int j = 0; j < fa.length(); ++j)
    sum += fa[j];

  for (int j = 0; j < fa.length(); ++j)
    fa[j] /= sum;
}

void LinkageEquilibriumDag::enforceMinFrequency(QList<double> &alleleFreq, double minFreq)
{
  bool changedFreq = false;
  for (int j = 0; j < alleleFreq.length(); ++j) {
    if (alleleFreq[j] < minFreq) {
      alleleFreq[j] = minFreq;
      changedFreq = true;
    }
  }
  if (changedFreq)
    divideEntriesBySum(alleleFreq);
}

int LinkageEquilibriumDag::outEdge(int level, int parentNode, int outEdge) const
{
  checkParentNode(level, parentNode);
  checkEdge(level, outEdge);
  return outEdge;
}

int LinkageEquilibriumDag::nInEdges(int level, int childNode) const
{
  checkLevel(level);
  return _alleleFreq[level].length();
}

int LinkageEquilibriumDag::inEdge(int level, int childNode, int inEdge) const
{
  checkEdge(level, inEdge);
  Q_ASSERT_X(childNode == 0, "LinkageEquilibriumDag::inEdge", "childNode != 0");
  return inEdge;
}

bool LinkageEquilibriumDag::isChildOf(int parentLevel, int parentEdge, int childEdge) const
{
  checkEdge(parentLevel, parentEdge);
  checkEdge(parentLevel + 1, childEdge);
  return true;
}

QList<double> LinkageEquilibriumDag::posArray() const
{
  QList<double> pos;
  for (int j = 0; j < _alleleFreq.length(); ++j) {
    double condEdgeProb = 0.0;
    for (int a = 0; a < _alleleFreq[j].length(); ++a)
      condEdgeProb += _alleleFreq[j][a] * _alleleFreq[j][a];

    if (j == 0)
      pos.append(-log10(condEdgeProb));
    else
      pos.append(pos[j - 1] - log10(condEdgeProb));
  }
  return pos;
}

void LinkageEquilibriumDag::checkLevel(int level) const
{
  Q_ASSERT_X(level >= 0 && level < _alleleFreq.length(),
             "LinkageEquilibriumDag::checkLevel",
             "invalid level");
}

void LinkageEquilibriumDag::checkEdge(int level, int edge) const
{
  Q_ASSERT_X(edge >= 0 && edge < _alleleFreq[level].length(),
             "LinkageEquilibriumDag::checkEdge",
             "invalid edge for this level");
}

void LinkageEquilibriumDag::checkParentNode(int level, int node) const
{
  checkLevel(level);
  Q_ASSERT_X(node == 0, "LinkageEquilibriumDag::checkParentNode", "invalid parent node (!= 0)");
}
