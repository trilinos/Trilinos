#ifndef ML_SMOOTHER_H
#define ML_SMOOTHER_H
#include "ml_config.h"
#include <iostream>

using namespace std;

namespace MLAPI {

class Vector;
class Space;

class Smoother {

public:

  virtual ~Smoother()
  {}

  virtual int ApplyInverse(const Vector& lhs, Vector& rhs) const = 0;

  virtual Vector& ApplyInverse(const Vector& lhs) const = 0;

  virtual Space& RangeSpace() const = 0;

  virtual Space& DomainSpace() const = 0;

}; // Smoother

} // namespace MLAPI
#endif // ML_SMOOTHER_H
