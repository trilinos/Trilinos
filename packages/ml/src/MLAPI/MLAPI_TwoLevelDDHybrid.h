#ifndef MLAPI_TWOLEVELCYCLE_H
#define MLAPI_TWOLEVELCYCLE_H

#include "ml_config.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Vector.h"
#include "MLAPI_Smoother.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Preconditioner.h"
#include "Epetra_MultiVector.h"

namespace MLAPI {

class TwoLevelDDHybrid : public Preconditioner {

public:

  TwoLevelDDHybrid(const Operator* A_f, 
                   const Smoother* Spre_f, const Smoother* Spost_f,
                   const Smoother* S_c,
                   const Operator* R, const Operator *P,
                   int MaxIters = 1) :
    A_f_(*A_f),
    R_(*R),
    P_(*P),
    Spre_f_(*Spre_f),
    Spost_f_(*Spost_f),
    S_c_(*S_c),
    MaxIters_(MaxIters)
  {}
      
  int Solve(const Vector& b_f, Vector& x_f) const
  {
    
    x_f = 0.0;

    Vector r_f(A_f_.DomainSpace());
    Vector r_c(S_c_.DomainSpace());

    for (int i = 0 ; i < MaxIters_ ; ++i) {
      // compute residual
      r_f = b_f; //- A_f_ * x_f;
      // apply smoother with zero initial guess
      x_f = Spre_f_ / r_f;
      // new residual
      r_f = b_f - A_f_ * x_f;
      // restrict to coarse
      r_c = R_ * r_f;
      // solve coarse problem (with smoother)
      r_c = S_c_ / r_c;
      // prolongate back and add to solution
      x_f = x_f + P_ * r_c;
      // apply post-smoother
      Spost_f_.ApplyInverse(b_f,x_f);;
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
  const Smoother& Spre_f_;
  const Smoother& Spost_f_;
  const Smoother& S_c_;
  const Operator& R_;
  const Operator& P_;
  int MaxIters_;

};
} // namespace MLAPI

#endif
