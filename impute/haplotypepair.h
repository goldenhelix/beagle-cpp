/* Copyright 2016 Golden Helix, Inc. */
#ifndef HAPLOTYPEPAIR_H
#define HAPLOTYPEPAIR_H

#include <QString>
#include <QPair>

class HaplotypePair
{
public:
  void setData(QPair<QString,QString> data);

  QString a() const { return _a; }
  QString b() const { return _b; }
private:
  QString _a;
  QString _b;
};

#endif
