/* Copyright notice. */
#ifndef SINGLEPERMITTEDSTATES_H
#define SINGLEPERMITTEDSTATES_H

#include "impute/ibshapsegments.h"

template <class E>
class SortedStartSet
{
 public:
  SortedStartSet() {}

  bool insert(const E &element)
  {
    if(_startMap.contains(element.start()))
      return false;
    else
    {
      _startMap.insert(element.start(), element);
      return true;
    }
  }

  bool contains(const E &element) const { return _startMap.contains(element.start()); }

  bool remove(const E &element)
  {
    if(_startMap.contains(element.start()))
    {
      _startMap.remove(element.start());
      return true;
    }
    else
      return false;
  }

  QMapIterator<int, E> iterator() const { return QMapIterator<int, E>(_startMap); }

  void clear() { _startMap.clear(); }

  QList<E> wholeList() const { return _startMap.values(); }  // (Needed only for SortedStartSet)

private:
  QMap<int, E> _startMap;
};

template <class E>
class SortedEndSet
{
public:
  SortedEndSet() {}

  bool insert(const E &element)
  {
    if(_endMap.contains(element.end()))
      return false;
    else
    {
      _endMap.insert(element.end(), element);
      return true;
    }
  }

  bool contains(const E &element) const { return _endMap.contains(element.end()); }

  bool remove(const E &element)
  {
    if(_endMap.contains(element.end()))
    {
      _endMap.remove(element.end());
      return true;
    }
    else
      return false;
  }

  QMapIterator<int, E> iterator() const { return QMapIterator<int, E>(_endMap); }

  void clear() { _endMap.clear(); }

private:
  QMap<int, E> _endMap;
};

template <class E>
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
    bool addedEnd = sortedEnd.insert(element);

    Q_ASSERT_X(addedStart == addedEnd,
               "Node<E>::add",
               "addedStart != addedEnd");

    return addedStart;
  }

  bool contains(const E &element) const
  {
    bool startContains = sortedStart.contains(element);

    Q_ASSERT_X(startContains == sortedEnd.contains(element),
               "Node<E>::contains",
               "startContains != sortedEnd.contains(element)");

    return startContains;
  }

  bool remove(const E &element)
  {
    bool removedStart = sortedStart.remove(element);
    bool removedEnd = sortedEnd.remove(element);

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
      QMapIterator<int, E> it = sortedStart.iterator();
      while (it.hasNext() && finished==false) {
        const E &e = it.next().value();
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
      QMapIterator<int, E> it = sortedEnd.iterator();
      while (it.hasNext() && finished==false) {
        const E &e = it.next().value();
        if (e.end() >= point) {
          collection.append(e);
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
      QMapIterator<int, E> it = sortedStart.iterator();
      while (it.hasNext() && finished==false) {
        const E &e = it.next().value();
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
      QMapIterator<int, E> it = sortedEnd.iterator();
      while (it.hasNext() && finished==false) {
        const E &e = it.next().value();
        if (start <= e.end()) {
          collection.append(e);
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
    QMapIterator<int, E> it = sortedStart.iterator();
    while (it.hasNext() && finished==false) {
      const E &e = it.next().value();
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

  void clear() {
    sortedStart.clear();
    sortedEnd.clear();
  }

  int center;
  SortedStartSet<E> sortedStart;
  SortedEndSet<E> sortedEnd;
  Node<E> *parent;
  Node<E> *leftChild;
  Node<E> *rightChild;
};

/**
 * (Templated) class {@code CenteredIntIntervalTree} implements a
 * centered interval tree that stores {@code IntInterval}-like
 * objects--specifically, {@code HapSegment} objects.
 *
 * "Instances of class {@code CenteredIntIntervalTree} are not thread-safe."

 * @param <E> the class of objects ({@code IntInterval} or {@code
 * HapSegment}) stored by {@code this}.
 */
template <class E>
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
               "CenteredIntIntervalTree<E>::add",
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

  void intersectPart(int start, int end, QList<E> &collection) const { intersectPart(_root, start, end, collection); }

  void intersectAll(int start, int end, QList<E> &collection) const { intersectAll(_root, start, end, collection); }

  bool isEmpty() const { return _size==0; }

  int size() const { return _size; }

  QList<E> toQList() const
  {
    QList<E> list;
    toQList(_root, list);
    return list;
  }

private:
  void initializeRootNode()
  {
    Q_ASSERT_X(_end >= _start,
               "CenteredIntIntervalTree const.",
               "_end < _start");

    int length = (_end - _start + 1);
    int center = _start + (length/2);
    _root = new Node<E>(center);
  }

  void clear(Node<E> *&tree)
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

  bool add(Node<E> *tree, const E &element)
  {
    if (element.end() < tree->center) {
      if (!tree->leftChild) {
        int nextOffset = findNextOffset(tree);
        tree->leftChild = new Node<E>(tree->center - nextOffset);
        tree->leftChild->parent = tree;
      }
      return add(tree->leftChild, element);
    }
    else if (element.start() > tree->center) {
      if (!tree->rightChild) {
        int nextOffset = findNextOffset(tree);
        tree->rightChild = new Node<E>(tree->center + nextOffset);
        tree->rightChild->parent = tree;
      }
      return add(tree->rightChild, element);
    }
    else {
      return tree->add(element);
    }
  }

  int findNextOffset(Node<E> *node) const {
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

  bool contains(Node<E> *tree, const E &element) const {
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

  bool remove(Node<E> *tree, const E &element) {
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

  void intersect(Node<E> *tree, int point, QList<E> &collection) const {
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

  void intersectPart(Node<E> *tree, int start, int end, QList<E> &collection) const {
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

  void intersectAll(Node<E> *tree, int start, int end, QList<E> &collection) const {
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

  void toQList(Node<E> *tree, QList<E> &list) const {
    if (!tree) {
      return;
    }
    toQList(tree->leftChild, list);
    list.append(tree->sortedStart.wholeList());
    toQList(tree->rightChild, list);
  }

  int _start;
  int _end;
  int _size;
  Node<E> *_root;
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
  void convertToIndices(int marker, const CenteredIntIntervalTree<HapSegment> &tree,
                        IndexSet &set) const;

  const RestrictedDag &_rdag;

  int _nMarkers;
  IndexSet _indices1;
  IndexSet _indices2;
  CenteredIntIntervalTree<HapSegment> _tree1;
  CenteredIntIntervalTree<HapSegment> _tree2;

  int _marker;
  int _i1;
  int _i2;
  int _edge1;
  int _edge2;
  bool _rev;
};

#endif
