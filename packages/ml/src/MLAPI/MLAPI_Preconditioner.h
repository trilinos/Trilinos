#ifndef MLAPI_PRECONDITIONER_H
#define MLAPI_PRECONDITIONER_H

namespace MLAPI {

class DoubleVector;
class Space;

class Preconditioner {

public:

  virtual int Solve(const DoubleVector& LHS, DoubleVector& RHS) const = 0;

  virtual const Space& DomainSpace() const = 0;

  virtual const Space& RangeSpace() const = 0;
};
} // namespace MLAPI

#endif
