#ifndef MLAPI_TWOLEVELDDADDITIVE_H
#define MLAPI_TWOLEVELDDADDITIVE_H

#include "ml_include.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Preconditioner.h"

namespace MLAPI {

class DoubleVector;
class Operator;
class InverseOperator;

class TwoLevelDDAdditive : public Preconditioner {

public:

  TwoLevelDDAdditive(const Operator* A_f,
                     const InverseOperator* S_f, const InverseOperator* S_c,
                     const Operator* R, const Operator *P,
                     int MaxIters = 1) :
    A_f_(*A_f),
    S_f_(*S_f),
    S_c_(*S_c),
    R_(*R),
    P_(*P),
    MaxIters_(MaxIters)
  {}
      
  int Solve(const DoubleVector& r_f, DoubleVector& x_f) const
  {
    
    // FIXME: here I suppose that the starting solution is zero

    DoubleVector r_c(S_c_.DomainSpace());

    // apply one-level preconditioner
    x_f = S_f_ / r_f;
    // restrict to coarse
    r_c = R_ * r_f;
    // solve coarse problem (with smoother)
    r_c = S_c_ / r_c;
    // prolongate back and add to solution
    x_f = x_f + P_ * r_c;

    return(0);
  }

  const Space& DomainSpace() const {
    return(A_f_.DomainSpace());
  }

  const Space& RangeSpace() const {
    return(A_f_.RangeSpace());
  }

private:
  const Operator& A_f_;
  const Operator& R_;
  const Operator& P_;
  const InverseOperator& S_f_;
  const InverseOperator& S_c_;
  int MaxIters_;

};
} // namespace MLAPI

#endif
