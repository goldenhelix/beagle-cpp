#include "haplotypepair.h"

void HaplotypePair::setData(QPair<QString,QString> data)
{
  _a = data.first;
  _b = data.second;
}
