#ifndef MLAPI_PRECONDITIONER_H
#define MLAPI_PRECONDITIONER_H

#include "ml_config.h"

namespace MLAPI {

class Vector;
class Space;

class Preconditioner {

public:

  virtual int Solve(const Vector& LHS, Vector& RHS) const = 0;

  virtual const Space& DomainSpace() const = 0;

  virtual const Space& RangeSpace() const = 0;
};
} // namespace MLAPI

#endif
