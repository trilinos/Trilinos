#ifndef MLAPI_TWOLEVELCYCLE_H
#define MLAPI_TWOLEVELCYCLE_H

#include "ml_include.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_Preconditioner.h"

namespace MLAPI {

class DoubleVector;
class Operator;
class Smoother;

class TwoLevelDDHybrid : public Preconditioner {

public:

  TwoLevelDDHybrid(const Operator* A_f, 
                   const Smoother* Spre_f, const Smoother* Spost_f,
                   const Smoother* S_c,
                   const Operator* R, const Operator *P) :
    A_f_(*A_f),
    R_(*R),
    P_(*P),
    Spre_f_(Spre_f),
    Spost_f_(Spost_f),
    S_c_(*S_c)
  {}
      
  int Solve(const DoubleVector& b_f, DoubleVector& x_f) const
  {
    
    DoubleVector r_f(A_f_.DomainSpace());
    DoubleVector r_c(S_c_.DomainSpace());

    // apply smoother with zero initial guess
    if (Spre_f_ != 0)
      x_f = (*Spre_f_) / b_f;
    // new residual
    r_f = b_f - A_f_ * x_f;
    // restrict to coarse
    r_c = R_ * r_f;
    // solve coarse problem
    r_c = S_c_ / r_c;
    // prolongate back and add to solution
    x_f = x_f + P_ * r_c;
    // apply post-smoother if available
    if (Spost_f_)
      Spost_f_->ApplyInverse(b_f,x_f);;

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
  const Smoother* Spre_f_;
  const Smoother* Spost_f_;
  const Smoother& S_c_;
  const Operator& R_;
  const Operator& P_;

};
} // namespace MLAPI

#endif
