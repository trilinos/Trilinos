#ifndef MLAPI_TWOLEVELDDADDITIVE_H
#define MLAPI_TWOLEVELDDADDITIVE_H

#include "ml_config.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Smoother.h"
#include "MLAPI_Vector.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Preconditioner.h"
#include "Epetra_MultiVector.h"

namespace MLAPI {

class TwoLevelDDAdditive : public Preconditioner {

public:

  TwoLevelDDAdditive(const Operator* A_f,
                     const Smoother* S_f, const Smoother* S_c,
                     const Operator* R, const Operator *P,
                     int MaxIters = 1) :
    A_f_(*A_f),
    S_f_(*S_f),
    S_c_(*S_c),
    R_(*R),
    P_(*P),
    MaxIters_(MaxIters)
  {}
      
  int Solve(const Vector& r_f, Vector& x_f) const
  {
    
    // FIXME: here I suppose that the starting solution is zero

    Vector r_c(S_c_.DomainSpace());

    for (int i = 0 ; i < MaxIters_ ; ++i) {
      // apply one-level preconditioner
      x_f = S_f_ / r_f;
      // restrict to coarse
      r_c = R_ * r_f;
      // solve coarse problem (with smoother)
      r_c = S_c_ / r_c;
      // prolongate back and add to solution
      x_f = x_f + P_ * r_c;
    }

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
  const Smoother& S_f_;
  const Smoother& S_c_;
  int MaxIters_;

};
} // namespace MLAPI

#endif
