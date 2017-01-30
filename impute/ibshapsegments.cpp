#include "impute/ibshapsegments.h"

#define INDEXMAP_NIL_VALUE   INT_MIN
#define INDEXMAP_MAX_VALUE   INT_MAX
#define INDEXMAP_DEL_VALUE       -67

#include <QMap>

/**
 * Class {@code IndexMap} is a map whose keys are a bounded set of
 * non-negative integers and whose values are integers.
 *
 * Class {@code IndexMap} supports a {@code clear()} method, but it does not
 * support a {@code remove()} method.
 *
 * "Class {@code IndexMap} is not thread-safe."
 */
class IndexMap
{
public:

  /**
   * Creates a new instance of {@code IndexMap} whose {@code nil()} method
   * will return the specified {@code nil} value.
   * @param maxKey the maximum key
   * @param nil the value that will be returned by the instance's
   * {@code get()} method if a key has no assigned value
   */
  IndexMap(int maxKey, int nil);

  /**
   * Copy constructor.
   */
  IndexMap(const IndexMap &other) : _nil(other._nil), _values(other._values),
                                    _keys(other._keys), _size(other._size) {}

  /**
   * Returns the value that is returned by {@code this.get()} if
   * a key has no assigned value.
   */
  int nil() const { return _nil; }

  /**
   * Adds the specified key and value to the map. If the map
   * contains a value for the specified key when the method is invoked,
   * the old value is replaced by the specified value.
   *
   * @param key the key
   * @param value the value
   * @return the previous value associated with {@code key}, or
   * {@code this.nil()} if no such previous value exists
   */
  int put(int key, int value);

  /**
   * Returns the value associated with the specified key, or
   * {@code this.nil()} if the specified key is not contained in the map.
   * @param key the key
   * @return the value associated with the specified key, or
   * {@code this.nil()} if the specified key is not contained in the map.
   */
  int operator[](int key) const { return _values[key]; }

  /**
   * Returns the number of key-value pairs in the map.
   */
  int size() const { return _size; }

  /**
   * Returns the maximum key.
   */
  int maxKey() const { return _keys.length() - 1; }

  /**
   * Removes all key-value pairs from the map.
   */
  void clear();

  /**
   * Returns the specified key in an enumeration of the keys in the map.
   * @param index an index of an element in the enumeration
   */
  int enumeratedKey(int index) const { return _keys[index]; }

  /**
   * Returns the value associated with the specified key
   * in an enumeration of the keys in the map.
   * If {@code (index >= 0 && index < this.size())}, then the returned value
   * will satisfy:
   * {@code (*this)[this.enumeratedKey(index)]=(*this).enumeratedValue(index)}.
   * @param index an index of an element in the enumeration
   * @return the value associated with the specified key
   * in an enumeration of the keys in the map
   */
  int enumeratedValue(int index) const { return _values[_keys[index]]; }

private:
  int _nil;
  QVector<int> _values;
  QVector<int> _keys;
  int _size;
};

class Haplotype
{
public:
  Haplotype(const SampleHapPairs &haps, int hapIndex, int start, int end)
    : _haps(haps), _hapIndex(hapIndex), _start(start), _end(end) {}

  bool operator<(const Haplotype &other) const;

private:
  const SampleHapPairs& _haps;
  int _hapIndex;
  int _start;
  int _end;
};

IndexMap::IndexMap(int maxKey, int nil)
  : _size(0)
{
  Q_ASSERT_X(maxKey >= 0,
             "IndexMap::IndexMap",
             "maxKey < 0");

  _nil = nil;
  _values.fill(nil, maxKey+1);
  _keys.fill(0, maxKey+1);
}

int IndexMap::put(int key, int value)
{
  Q_ASSERT_X(value != _nil,
             "IndexMap::put",
             "value == nil()");

  int prevValue = _values[key];
  if (prevValue == _nil) {
    _keys[_size++] = key;
  }
  _values[key] = value;
  return prevValue;
}

void IndexMap::clear()
{
  for (int j=0; j < _size; ++j)
    _values[_keys[j]] = _nil;

  _size = 0;
}

bool Haplotype::operator< (const Haplotype &other) const
{
  // Optimizations made possible by our actual usage--we assume:
  //
  // (1) We're working with references to the same SampleHapPairs object;
  // (2) Both starts are the same.
  // (3) Both ends are the same.
  //
  // Our test is to compare the (possibly different) haplotypes at the
  // (usually or always) different haplotype indices over the same
  // marker range using the same haplotype list.

  int otherHapIndex = other._hapIndex;

  for(int m=_start; m < _end; m++)
  {
    if(_haps.allele(m, _hapIndex) != _haps.allele(m, otherHapIndex))
      return (_haps.allele(m, _hapIndex) < _haps.allele(m, otherHapIndex));
  }

  return false;    // All equal--therefore, neither one is less than the other.
}

int IbsHapSegUtility::lowerBoundIndex(const QList<double> &pos, int first, double value)
{
  int it, count, step;
  count = pos.length() - first;
 
  while (count > 0)
  {
    it = first; 
    step = count / 2; 
    it += step;

    if (pos[it] < value)
    {
      first = ++it; 
      count -= step + 1; 
    }
    else
      count = step;
  }

  return first;
}

HapSegment::HapSegment(int hap, int start, int end)
  :  _hap(hap), _start(start), _end(end)
{
  Q_ASSERT_X(start <= end,
             "HapSegment::HapSegment",
             "start > end");
}

bool HapSegment::operator< (const HapSegment &other) const
{
  if (_start != other._start)
    return (_start < other._start);
  else if (_end != other._end)
    return (_end < other._end);
  else
    return (_hap < other._hap);
}

HapSegmentES::HapSegmentES(const HapSegment &hs)
  : HapSegment(hs.hap(), hs.start(), hs.end()) {}

bool HapSegmentES::operator< (const HapSegmentES &other) const
{
  if (_end != other._end)
    return (_end > other._end);     // Means "less than" here in this context....
  else if (_start != other._start)
    return (_start < other._start);
  else
    return (_hap < other._hap);
}


void IbsHapSegments::initialize(const SampleHapPairs &haps, const QList<double> &pos, double minIbsLength)
{
  checkArguments(haps, pos, minIbsLength);
  _haps = haps;
  _pos = pos;
  _minIbsLength = minIbsLength;
  findWindowStarts();
  findIdSets();
  /// globalRsb.dumpWs(_windowStarts);
  /// globalRsb.dumpIdSets(_idSets);
}

void IbsHapSegments::checkArguments(const SampleHapPairs &haps, const QList<double> &pos,
                                    double minIbsLength) const
{
  Q_ASSERT_X(minIbsLength > 0.0f,
             "IbsHapSegments::checkArguments",
             "minIbsLength <= 0.0f");

  Q_ASSERT_X(haps.nMarkers() == pos.length(),
             "IbsHapSegments::checkArguments",
             "haps.nMarkers()!= pos.length");

  Q_ASSERT_X(pos[0] >= 0,
             "IbsHapSegments::checkArguments",
             "pos[0]<0");

  for (int j=1; j<pos.length(); ++j)
  {
    Q_ASSERT_X(pos[j] >= pos[j-1],
               "IbsHapSegments::checkArguments",
               "Positions are not non-decreasing.");
  }
}

void IbsHapSegments::findWindowStarts()
{
  double step = _minIbsLength/2.0f;
  int index = 0;
  do {
    _windowStarts.append(index);
    double nextPos = _pos[index] + step;
    index = IbsHapSegUtility::lowerBoundIndex(_pos, index, nextPos);
  } while (index < _pos.length());
}


void IbsHapSegments::findIdSets()
{
  int wsLength = _windowStarts.length();
  int nMarkers = _haps.nMarkers();
  for(int index = 0; index < wsLength; index++)
    fillOneIdSet(_windowStarts[index],  (index+1 < _windowStarts.length()) ? _windowStarts[index+1] : nMarkers );
}

void IbsHapSegments::fillOneIdSet(int ipFirst, int ipSecond)
{
  int nHaps = _haps.nHaps();
  QMap< Haplotype, QList<int> > idSetMap;
  for(int hapIndex=0; hapIndex < nHaps; hapIndex++)
  {
    Haplotype newHap(_haps, hapIndex, ipFirst, ipSecond);
    idSetMap[newHap].append(hapIndex);
  }   

  QVector< QList<int> > idSet(nHaps);
  foreach(QList<int> hapIndexList, idSetMap.values())
  {
    foreach(int i, hapIndexList)
      idSet[i] = hapIndexList;
  }
  _idSets.append(idSet);
}

void IbsHapSegments::find(QList<HapSegment> &outHapSegs, int hap) const
{
  IndexMap prev(_haps.nHaps()-1, INDEXMAP_NIL_VALUE);
  IndexMap next(_haps.nHaps()-1, INDEXMAP_NIL_VALUE);
  int window = 0;
  matches(hap, window, prev);
  while (++window < _idSets.length())
  {
    matches(hap, window, next);
    extend(prev, next);
    save(hap, prev, _windowStarts[window], outHapSegs);

    prev = next;
    next.clear();
  }

  save(hap, prev, _haps.nMarkers(), outHapSegs);
}

void IbsHapSegments::filteredFind(QList<HapSegment> &outHapSegs, int hap) const
{
  IndexMap prev(_haps.nHaps()-1, INDEXMAP_NIL_VALUE);
  IndexMap next(_haps.nHaps()-1, INDEXMAP_NIL_VALUE);
  int window = 0;
  matches(hap, window, prev);
  while (++window < _idSets.length())
  {
    matches(hap, window, next);
    int minExtendedStartWindow = extend(prev, next);
    filteredSave(hap, prev, minExtendedStartWindow,
                 _windowStarts[window], outHapSegs);
    prev = next;
    next.clear();
  }

  filteredSave(hap, prev, window, _haps.nMarkers(), outHapSegs);
}

void IbsHapSegments::matches(int hap, int window, IndexMap &map) const
{
  Q_ASSERT_X(map.size() == 0,
             "IbsHapSegments::matches",
             "map.size() != 0");

  foreach (int h, _idSets[window][hap])
  {
    if (h!=hap)
      map.put(h, window);
  }
}

int IbsHapSegments::extend(IndexMap &prev, IndexMap &next) const
{
  int nil = next.nil();
  int minStart = INDEXMAP_MAX_VALUE;
  for (int i=0, n=next.size(); i<n; ++i)
  {
    int hap = next.enumeratedKey(i);
    int prevStart = prev[hap];
    if (prevStart != nil)
    {
      next.put(hap, prevStart);
      prev.put(hap, INDEXMAP_DEL_VALUE);
      if (prevStart < minStart)
      minStart = prevStart;
    }
  }
  return minStart;
}

void IbsHapSegments::save(int hap1, const IndexMap &prev,
                          int prevExclEnd, QList<HapSegment> &segments) const
{
  for (int i=0, n=prev.size(); i<n; ++i)
  {
    int hap2 = prev.enumeratedKey(i);
    int startWindow = prev.enumeratedValue(i);
    if (startWindow != INDEXMAP_DEL_VALUE)
    {
      int start = findStart(hap1, hap2, _windowStarts[startWindow]);
      int inclEnd = findInclusiveEnd(hap1, hap2, prevExclEnd);
      if ( (_pos[inclEnd] - _pos[start]) >= _minIbsLength)
        segments.append( HapSegment(hap2, start, inclEnd) );
    }
  }
}

void IbsHapSegments::filteredSave(int hap1, const IndexMap &prev,
                                  int minExtendedStartWindow,
                                  int prevExclEnd, QList<HapSegment> &segments) const
{
  for (int i=0, n=prev.size(); i<n; ++i)
  {
    int hap2 = prev.enumeratedKey(i);
    int startWindow = prev.enumeratedValue(i);
    if (startWindow != INDEXMAP_DEL_VALUE  &&  startWindow <= minExtendedStartWindow)
    {
      int start = findStart(hap1, hap2, _windowStarts[startWindow]);
      int inclEnd = findInclusiveEnd(hap1, hap2, prevExclEnd);
      if ( (_pos[inclEnd] - _pos[start]) >= _minIbsLength)
        segments.append( HapSegment(hap2, start, inclEnd) );
    }
  }
}


int IbsHapSegments::findStart(int hap1, int hap2, int start) const
{
  while (start > 0  &&  _haps.allele(start-1, hap1) == _haps.allele(start-1, hap2))
    --start;

  return start;
}

int IbsHapSegments::findInclusiveEnd(int hap1, int hap2, int end) const
{
  while (end < _haps.nMarkers()  &&  _haps.allele(end, hap1) == _haps.allele(end, hap2))
    ++end;

  return end-1;
}

RSBDump globalRsb;

void RSBDump::setOkToDump(bool otd)
{
  /*
  _okToDump = otd;
  if(otd)
  {
    if (!_outHasBeenOpened)
    {
      _out.setFileName("rsbdump9k4.txt");
      if (!_out.open(QIODevice::WriteOnly)) {
        printf("Unable to write to the text file.\n");
        return;
      }
      _outHasBeenOpened = true;
    }
  }
  */
  _okToDump = false;
}

void RSBDump::setNewSample(int sampNum)
{
  if(_okToDump)
  {
    _out.flush();
    _sampleNum = sampNum;
    _out.write("SAMPLE ");
    _out.write(QByteArray::number(_sampleNum));
    _out.write(":\n");
    if (sampNum > 0)
      _okToDump = false;
  }
}

void RSBDump::dumpNext(int rv, int e1, int e2, int i1, int i2)
{
  if(_okToDump)
  {
    _out.write(" rv");
    _out.write(QByteArray::number(rv));
    _out.write(" e1");
    _out.write(QByteArray::number(e1));
    _out.write(" e2");
    _out.write(QByteArray::number(e2));
    _out.write(" i1");
    _out.write(QByteArray::number(i1));
    _out.write(" i2");
    _out.write(QByteArray::number(i2));
    cr();
  }
}

void RSBDump::dumpFwd1(int m, int sy1, int sy2, float emp)
{
  if(_okToDump)
  {
    _out.write(" m");
    _out.write(QByteArray::number(m));
    _out.write(" sy1_");
    _out.write(QByteArray::number(sy1));
    _out.write(" sy2_");
    _out.write(QByteArray::number(sy2));
    _out.write(" emp: ");
    _out.write(QByteArray::number(emp, 'f', 3));
  }
}

void RSBDump::dumpFwd2(float pr, float r0, float r1, float r2, float ep1, float ep2)
{
  if(_okToDump)
  {
    _out.write(" pr: ");
    _out.write(QByteArray::number(pr, 'f', 3));
    _out.write(" r0: ");
    _out.write(QByteArray::number(r0, 'f', 3));
    _out.write(" r1: ");
    _out.write(QByteArray::number(r1, 'f', 3));
    _out.write(" r2: ");
    _out.write(QByteArray::number(r2, 'f', 3));
    _out.write(" ep1: ");
    _out.write(QByteArray::number(ep1, 'f', 3));
    _out.write(" ep2: ");
    _out.write(QByteArray::number(ep2, 'f', 3));
  }
}

void RSBDump::dumpFwd3(int pn1, int pn2, float pnp1, float pnp2)
{
  if(_okToDump)
  {
    _out.write(" pn1: ");
    _out.write(QByteArray::number(pn1));
    _out.write(" pn2: ");
    _out.write(QByteArray::number(pn2));
    _out.write(" pnp1: ");
    _out.write(QByteArray::number(pnp1, 'f', 3));
  }
}

void RSBDump::cr()
{
  if(_okToDump)
    _out.write("\n");
}

void RSBDump::dumpSCD(int m, int e1, int e2, int s1, int s2, int node1, int node2, float p1, float p2,
                      float maxSum, float fv, float fs, float gl)
{
  _out.write(" m: ");
  _out.write(QByteArray::number(m));
  _out.write(" e1: ");
  _out.write(QByteArray::number(e1));
  _out.write(" e2: ");
  _out.write(QByteArray::number(e2));
  _out.write(" s1: ");
  _out.write(QByteArray::number(s1));
  _out.write(" s2: ");
  _out.write(QByteArray::number(s2));
  _out.write(" node1: ");
  _out.write(QByteArray::number(node1));
  _out.write(" node2: ");
  _out.write(QByteArray::number(node2));
  _out.write(" p1: ");
  _out.write(QByteArray::number(p1, 'f', 3));
  _out.write(" p2: ");
  _out.write(QByteArray::number(p2, 'f', 3));
  _out.write(" maxSum: ");
  _out.write(QByteArray::number(maxSum, 'f', 3));
  _out.write(" fv: ");
  _out.write(QByteArray::number(fv, 'f', 3));
  _out.write(" fs: ");
  _out.write(QByteArray::number(fs, 'f', 3));
  _out.write(" gl: ");
  _out.write(QByteArray::number(gl, 'f', 3));
  _out.write("\n");
  _out.flush();
}

void RSBDump::dumpWs(const QList<int> &ws)
{
  if(_okToDump)
  {
    _out.write("Window Starts:");

    for(int i=0; i<ws.length(); i++)
    {
      if(i % 20 == 0)
        cr();
      else
        _out.write("/");
      _out.write(QByteArray::number(ws[i]));
    }
    cr();
  }
}

void RSBDump::dumpIdSets(const QList< QVector< QList<int> > > &idSets)
{
  if(_okToDump)
  {
    for(int o=0; o < idSets.length(); o++)
    {
      _out.write("ID Set ");
      _out.write(QByteArray::number(o));
      _out.write(":");
      cr();
      QVector< QList<int> > middle = idSets[o];
      for(int m=0; m < middle.length(); m++)
      {
        QList<int> inner = middle[m];
        for(int i=0; i < inner.length(); i++)
        {
          _out.write(" ");
          _out.write(QByteArray::number(inner[i]));
        }
        cr();
      }
    }
  }
}

void RSBDump::dumpHSegs(int hap, const QList<HapSegment> &ibsSegs)
{
  if (_okToDump)
  {
    _out.write("Hap ");
    _out.write(QByteArray::number(hap));
    _out.write(":");
    cr();
    foreach(HapSegment seg, ibsSegs)
    {
      _out.write(" ");
      _out.write(QByteArray::number(seg.hap()));
      _out.write("/");
      _out.write(QByteArray::number(seg.start()));
      _out.write("/");
      _out.write(QByteArray::number(seg.end()));
      cr();
    }
    _out.flush();
  }
}

void RSBDump::treeHeader(int hap, int start, int end, int size)
{
  _out.write("Tree from Hap ");
  _out.write(QByteArray::number(hap));
  _out.write("--Start ");
  _out.write(QByteArray::number(start));
  _out.write(" End ");
  _out.write(QByteArray::number(end));
  _out.write(" Size ");
  _out.write(QByteArray::number(size));
  _out.write(":");
  cr();
}

void RSBDump::dumpTSegs(int center, char *type, int depth, const QList<HapSegment> &ibsSegs)
{
  for (int margin = 0; margin < depth - 1; margin++)
    _out.write(" ");
  _out.write(type);
  _out.write(" Center: ");
  _out.write(QByteArray::number(center));
  cr();
  foreach(HapSegment seg, ibsSegs)
  {
    for (int margin = 0; margin < depth - 1; margin++)
      _out.write(" ");
    _out.write(type);
    _out.write(" ");
    _out.write(QByteArray::number(seg.hap()));
    _out.write("/");
    _out.write(QByteArray::number(seg.start()));
    _out.write("/");
    _out.write(QByteArray::number(seg.end()));
    cr();
  }
  _out.flush();
}

void RSBDump::dumpMSegs(int marker, const QList<HapSegment> &ibsSegs)
{
  if (_okToDump)
  {
    _out.write("Marker ");
    _out.write(QByteArray::number(marker));
    _out.write(":");
    cr();
    foreach(HapSegment seg, ibsSegs)
    {
      _out.write(" ");
      _out.write(QByteArray::number(seg.hap()));
      _out.write("/");
      _out.write(QByteArray::number(seg.start()));
      _out.write("/");
      _out.write(QByteArray::number(seg.end()));
      cr();
    }
    _out.flush();
  }
}

void RSBDump::dumpModifier(char *svse, int maxStart, int minEnd, const QList<HapSegment> &ibsSegs)
{
  if (_okToDump)
  {
    _out.write(svse);
    _out.write(" maxStart ");
    _out.write(QByteArray::number(maxStart));
    _out.write(" minEnd ");
    _out.write(QByteArray::number(minEnd));
    _out.write(":  ");
    foreach(HapSegment seg, ibsSegs)
    {
      _out.write(" ");
      _out.write(QByteArray::number(seg.hap()));
      _out.write("/");
      _out.write(QByteArray::number(seg.start()));
      _out.write("/");
      _out.write(QByteArray::number(seg.end()));
    }
    cr();
    _out.flush();
  }
}

