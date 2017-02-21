/* Copyright notice. */
#ifndef SINGLEPERMITTEDSTATES_H
#define SINGLEPERMITTEDSTATES_H

#include "impute/ibshapsegments.h"

template <class E>
class SortedSet
{
 public:
  SortedSet() {}

  bool insert(const E &element)
  {
    if(_sortMap.contains(element))
      return false;
    else
    {
      _sortMap.insert(element, element);
      return true;
    }
  }

  bool remove(const E &element)
  {
    if(_sortMap.contains(element))
    {
      _sortMap.remove(element);
      return true;
    }
    else
      return false;
  }

  void clear() { _sortMap.clear(); }

  bool contains(const E &element) const { return _sortMap.contains(element); }
  QMapIterator<E, E> iterator() const { return QMapIterator<E, E>(_sortMap); }
  QList<E> wholeList() const { return _sortMap.values(); }

  bool hasIntersectAll(int start, int end) const
  {
    QMap<E, E>::const_iterator i = _sortMap.constBegin();

    while(i != _sortMap.constEnd()  &&  i.value().start() <= start)
    {
      if(i.value().end() >= end)
        return true;

      i = _sortMap.lowerBound(E(0, i.value().start(), end));  // Hap, start, end
    }

    return false;
  }

private:
  QMap<E, E> _sortMap;
};

template <class E, class F>
class Node
{
public:

  Node(int center) : center(center), leftChild(0), rightChild(0) {}

  bool add(const E &element)
  {
    Q_ASSERT_X(element.start() <= center  &&  element.end() >= center,
               "Node<E>::add",
               "element does not overlap center.");

    bool addedStart = sortedStart.insert(element);
    F elementF(element);
    bool addedEnd = sortedEnd.insert(elementF);

    Q_ASSERT_X(addedStart == addedEnd,
               "Node<E>::add",
               "addedStart != addedEnd");

    return addedStart;
  }

  bool contains(const E &element) const
  {
    bool startContains = sortedStart.contains(element);
    F elementF(element);
    Q_ASSERT_X(startContains == sortedEnd.contains(elementF),
               "Node<E>::contains",
               "startContains != sortedEnd.contains(elementF)");

    return startContains;
  }

  bool remove(const E &element)
  {
    bool removedStart = sortedStart.remove(element);
    F elementF(element);
    bool removedEnd = sortedEnd.remove(elementF);

    Q_ASSERT_X(removedStart == removedEnd,
               "Node<E>::remove",
               "removedStart != removedEnd");

    return removedStart;
  }

  void intersect(int point, QList<E> &collection) const
  {
    if (point <= center)
    {
      bool finished = false;
      QMapIterator<E, E> it = sortedStart.iterator();
      while (it.hasNext() && finished==false) {
        E e = it.next().value();
        if (e.start() <= point) {
          collection.append(e);
        }
        else {
          finished = true;
        }
      }
    }
    else {
      bool finished = false;
      QMapIterator<F, F> it = sortedEnd.iterator();
      while (it.hasNext() && finished==false) {
        F f = it.next().value();
        if (f.end() >= point) {
          const F &fp = f;
          collection.append(fp);   // This will downcast from an F to an E.
        }
        else {
          finished = true;
        }
      }
    }
  }

  // Optimized version:
  void intersect(int point, QVector<int> &collection) const
  {
    if (point <= center)
    {
      bool finished = false;
      QMapIterator<E, E> it = sortedStart.iterator();
      while (it.hasNext() && finished==false) {
        E e = it.next().value();
        if (e.start() <= point) {
          collection.append(e.hap());
        }
        else {
          finished = true;
        }
      }
    }
    else {
      bool finished = false;
      QMapIterator<F, F> it = sortedEnd.iterator();
      while (it.hasNext() && finished==false) {
        F f = it.next().value();
        if (f.end() >= point) {
          collection.append(f.hap());
        }
        else {
          finished = true;
        }
      }
    }
  }

  void intersectPart(int start, int end, QList<E> &collection) const
  {
    if (end < center) {
      bool finished = false;
      QMapIterator<E, E> it = sortedStart.iterator();
      while (it.hasNext() && finished==false) {
        E e = it.next().value();
        if (e.start() <= end) {
          collection.append(e);
        }
        else {
          finished = true;
        }
      }
    }
    else if (start > center) {
      bool finished = false;
      QMapIterator<F, F> it = sortedEnd.iterator();
      while (it.hasNext() && finished==false) {
        F f = it.next().value();
        if (start <= f.end()) {
          const F &fp = f;
          collection.append(fp);   // This will downcast from an F to an E.
        }
        else {
          finished = true;
        }
      }
    }
    else {
      collection.append(sortedStart.wholeList());
    }
  }

  void intersectAll(int start, int end, QList<E> &collection) const
  {
    bool finished = false;
    QMapIterator<E, E> it = sortedStart.iterator();
    while (it.hasNext() && finished==false) {
      E e = it.next().value();
      if (e.start() <= start) {
        if (e.end() >= end) {
          collection.append(e);
        }
      }
      else {
        finished = true;
      }
    }
  }

  bool hasIntersectAll(int start, int end) const
  {
    return sortedStart.hasIntersectAll(start, end);
  }

  void clear() {
    sortedStart.clear();
    sortedEnd.clear();
  }

  int center;
  SortedSet<E> sortedStart;
  SortedSet<F> sortedEnd;
  Node<E, F> *parent;
  Node<E, F> *leftChild;
  Node<E, F> *rightChild;
};

/**
 * (Templated) class {@code CenteredIntIntervalTree} implements a
 * centered interval tree that stores {@code IntInterval}-like
 * objects--specifically, {@code HapSegment} objects.
 *
 * "Instances of class {@code CenteredIntIntervalTree} are not thread-safe."

 * @param <E> the class of objects ({@code IntInterval} or {@code
 * HapSegment}) stored by {@code this}.
 * @param <F> the end-sorted counterpart of class {@code HapSegment},
 * namely, {@code HapSegmentES}.
 */
template <class E, class F>
class CenteredIntIntervalTree
{
public:

  /**
   * "Creates a new {@code CenteredIntIntervalTree} instance for the
   * specified range.
   * @param start the minimum start value (inclusive) for intervals stored in
   * this interval tree
   * @param end the maximum end value (inclusive) for intervals stored
   * in this interval tree"
   *
   * CenteredIntIntervalTree(int start, int end) : _start(start), _end(end), _size(0)
   * {
   *   initializeRootNode();
   * }
   */

  /**
   * Default constructor.
   */
  CenteredIntIntervalTree() : _start(0), _end(-1), _size(0), _root(0) {}

  /**
   * Use this constructor (for an automatic variable) instead of a
   * "SinglePermittedStates::getTree" method.
   */
  CenteredIntIntervalTree(int nMarkers, const QList<E> &ibsSegs) : _start(0), _size(0)
  {
    initialize(nMarkers, ibsSegs);
  }

  ~CenteredIntIntervalTree()
  {
    clear();
  }

  /**
   * Use this method (on a class variable) instead of a
   * "SinglePermittedStates::getTree" method.
   */
  void initialize(int nMarkers, const QList<E> &ibsSegs)
  {
    _end = nMarkers - 1;
    initializeRootNode();

    foreach(E hs, ibsSegs)
      add(hs);
  }

  int start() const { return _start; }

  int end() const { return _end; }

  void clear()
  {
    clear(_root);
    _size = 0;
  }

  bool add(const E &element)
  {
    Q_ASSERT_X(element.start() >= _start  &&  element.end() <= _end,
               "CenteredIntIntervalTree<E, F>::add",
               "element out of range.");

    bool added = add(_root, element);
    if (added)
      _size++;

    return added;
  }

  bool contains(const E &element) const { return contains(_root, element); }

  bool remove(const E &element)
  {
    bool removed = remove(_root, element);
    if (removed)
      _size--;

    return removed;
  }

  void intersect(int point, QList<E> &collection) const { intersect(_root, point, collection); }
  // Optimized overload:
  void intersect(int point, QVector<int> &collection) const { intersect(_root, point, collection); }

  void intersectPart(int start, int end, QList<E> &collection) const { intersectPart(_root, start, end, collection); }

  void intersectAll(int start, int end, QList<E> &collection) const { intersectAll(_root, start, end, collection); }

  bool hasIntersectAll(int start, int end) const { return hasIntersectAll(_root, start, end); }

  bool isEmpty() const { return _size==0; }

  int size() const { return _size; }

  void getDumpInfo(int &start, int &end, int &size, Node<E, F> *&root) const
  {
    start = _start;
    end = _end;
    size = _size;
    root = _root;
  }

  /// QList<E> toQList() const
  /// {
  ///   QList<E> list;
  ///   toQList(_root, list);
  ///   return list;
  /// }

private:
  void initializeRootNode()
  {
    Q_ASSERT_X(_end >= _start,
               "CenteredIntIntervalTree const.",
               "_end < _start");

    int length = (_end - _start + 1);
    int center = _start + (length/2);
    _root = new Node<E, F>(center);
    _root->parent = 0;
  }

  void clear(Node<E, F> *&tree)
  {
    if (!tree) {
      return;
    }
    tree->clear();
    clear(tree->leftChild);
    clear(tree->rightChild);
    delete tree;
    tree = 0;
  }

  bool add(Node<E, F> *tree, const E &element)
  {
    if (element.end() < tree->center) {
      if (!tree->leftChild) {
        int nextOffset = findNextOffset(tree);
        tree->leftChild = new Node<E, F>(tree->center - nextOffset);
        tree->leftChild->parent = tree;
      }
      return add(tree->leftChild, element);
    }
    else if (element.start() > tree->center) {
      if (!tree->rightChild) {
        int nextOffset = findNextOffset(tree);
        tree->rightChild = new Node<E, F>(tree->center + nextOffset);
        tree->rightChild->parent = tree;
      }
      return add(tree->rightChild, element);
    }
    else {
      return tree->add(element);
    }
  }

  int findNextOffset(Node<E, F> *node) const {
    int lastOffset;
    if (!node->parent) {
      lastOffset = (_end - _start + 1)/2;
    }
    else {
      lastOffset = abs(node->center - node->parent->center);
    }
    Q_ASSERT_X(lastOffset > 0,
               "CenteredIntIntervalTree::findNextOffset",
               "lastOffset <= 0");
    int offset = (lastOffset+1)/2;
    return offset;
  }

  bool contains(Node<E, F> *tree, const E &element) const {
    if (!tree) {
      return false;
    }
    else if (element.end() < tree->center) {
      return contains(tree->leftChild, element);
    }
    else if (element.start() > tree->center) {
      return contains(tree->rightChild, element);
    }
    else {
      return tree->contains(element);
    }
  }

  bool remove(Node<E, F> *tree, const E &element) {
    if (!tree) {
      return false;
    }
    else if (element.end() < tree->center) {
      return remove(tree->leftChild, element);
    }
    else if (element.start() > tree->center) {
      return remove(tree->rightChild, element);
    }
    else {
      return tree->remove(element);
    }
  }

  void intersect(Node<E, F> *tree, int point, QList<E> &collection) const {
    if (!tree) {
      return;
    }
    tree->intersect(point, collection);
    if (point < tree->center) {
      intersect(tree->leftChild, point, collection);
    }
    else if (point > tree->center) {
      intersect(tree->rightChild, point, collection);
    }
  }

  // Optimized overload:
  void intersect(Node<E, F> *tree, int point, QVector<int> &collection) const {
    if (!tree) {
      return;
    }
    tree->intersect(point, collection);
    if (point < tree->center) {
      intersect(tree->leftChild, point, collection);
    }
    else if (point > tree->center) {
      intersect(tree->rightChild, point, collection);
    }
  }

  void intersectPart(Node<E, F> *tree, int start, int end, QList<E> &collection) const {
    if (!tree) {
      return;
    }
    tree->intersectPart(start, end, collection);
    if (start < tree->center) {
      intersectPart(tree->leftChild, start, end, collection);
    }
    if (end > tree->center) {
      intersectPart(tree->rightChild, start, end, collection);
    }
  }

  void intersectAll(Node<E, F> *tree, int start, int end, QList<E> &collection) const {
    if (!tree) {
      return;
    }
    tree->intersectAll(start, end, collection);
    if (end < tree->center) {
      intersectAll(tree->leftChild, start, end, collection);
    }
    if (start > tree->center) {
      intersectAll(tree->rightChild, start, end, collection);
    }
  }

  bool hasIntersectAll(Node<E, F> *tree, int start, int end) const {

    if (!tree)
      return false;

    if(tree->hasIntersectAll(start, end))
      return true;

    if (end < tree->center)
      return hasIntersectAll(tree->leftChild, start, end);

    if (start > tree->center)
      return hasIntersectAll(tree->rightChild, start, end);

    return false;
  }

  /// void toQList(Node<E, F> *tree, QList<E> &list) const {
  ///   if (!tree) {
  ///     return;
  ///   }
  ///   toQList(tree->leftChild, list);
  ///   list.append(tree->sortedStart.wholeList());
  ///   toQList(tree->rightChild, list);
  /// }

  int _start;
  int _end;
  int _size;
  Node<E, F> *_root;
};


/**
 * Class {@code IndexSet} is a set that stores non-negative indices that are
 * less than or equal to a specified maximum value.
 *
 * Class {@code IndexSet} supports a {@code clear()} method, but it does not
 * support a {@code remove()} method.
 * 
 * "Class {@code IndexSet} is not thread-safe."
 */
class IndexSet
{
public:

  /* Default constructor. */
  IndexSet() : _size(0) {}

  /**
   * Creates a new instance of {@code IndexSet} that can contain
   * non-negative integer indices that are less than or equal to the specified
   * maximum value.
   */
  void initialize(int max);

  /**
   * Adds the specified element to the set.
   *
   * @param element an element to add to this set.
   * @return {@code true} if the set was changed by the operation, and
   * {@code false} otherwise.
   */
  bool add(int element);

  /**
   * Returns {@code true} if the set contains the specified element,
   * and returns {@code false} otherwise.
   * @param element an element
   * @return {@code true} if the set contains the specified element
   */
  bool contains(int element) const { return _inSet[element]; }

  /**
   * Returns the number of elements in this set.
   */
  int size() const { return _size; }

  /**
   * Returns the maximum permitted element in the set.
   */
  int maxPermittedElement() const { return _indices.length()-1; }

  /**
   * Removes all elements from the set.
   */
  void clear();

  /**
   * Returns the specified element in an enumeration of the elements in the
   * set.
   * @param enumIndex an index of an element in the enumeration
   * @return the specified element in an enumeration of the elements in the
   * set
   */
  int enumeratedValue(int enumIndex) const;

  /**
   * Returns a QVector containing the (active) elements in this set.
   */
  QVector<int> toQVector() const { return _indices.mid(0, _size); }

private:
  QVector<bool> _inSet;
  QVector<int> _indices;
  int _size;
};

class RestrictedDag;

class SinglePermittedStates
{
public:
  SinglePermittedStates(const RestrictedDag &rdag, int sample);

  int nMarkers() const { return _nMarkers; }
  int marker() const { return _marker; }
  void setMarker(int marker);
  bool hasNext() const { return _i1 < _indices1.size(); }
  void next();
  int edge1() const { return _edge1; }
  int edge2() const { return _edge2; }

private:
  void extendSegment(QList<HapSegment> &extendedSegs, int hap, const QList<HapSegment> &ibsSegs) const;
  void convertToIndices(int marker, const CenteredIntIntervalTree<HapSegment, HapSegmentES> &tree,
                        IndexSet &set);

  int modifyStart(const HapSegment &targetHS, CenteredIntIntervalTree<HapSegment, HapSegmentES> &tree) const;
  int modifyEnd(const HapSegment &targetHS, CenteredIntIntervalTree<HapSegment, HapSegmentES> &tree) const;

  const RestrictedDag &_rdag;

  const QList<double>& _pos;
  double _ibdExtend;

  int _nMarkers;
  IndexSet _indices1;
  IndexSet _indices2;
  CenteredIntIntervalTree<HapSegment, HapSegmentES> _tree1;
  CenteredIntIntervalTree<HapSegment, HapSegmentES> _tree2;

  int _marker;
  int _i1;
  int _i2;
  int _edge1;
  int _edge2;
  bool _rev;

  QVector<int> _hsegHaps;
  int _maxHaps;
};

#endif
