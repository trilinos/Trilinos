#ifndef _TSF_PRECLINEARPROBLEM_H_
#define _TSF_PRECLINEARPROBLEM_H_

//! TSF_LinearProblem:  The Trilinos Virtual Linear Problem Class.
/*! The TSF_RDP_PrecLinearProblem class is a wrapper that encapsulates the 
  general information needed for solving a linear system of equations.  
*/

#include "TSF_RDP_LinearProblem.h"
#include "TSF_RDP_MultiVector.h"

class TSF_RDP_PrecLinearProblem: public virtual TSF_RDP_LinearProblem {
    
  public:

  //! TSF_PrecLinearProblem Destructor.
  virtual ~TSF_PrecLinearProblem(void){};

  //! Set Preconditioner \f$M\f$.
  virtual int SetPreconditioner(TSF_RDP_Preconditioner & M) = 0;

  //! Get Preconditioner \f$M\f$.
  virtual int GetPreconditioner(TSF_RDP_Preconditioner * & M) = 0;

  //! Returns the preconditioned residual(\f$ \|r\| = \|b-Ax\|\f$) for this problem using the current solution.
  virtual double PrecResidual() = 0;

  //! Returns the preconditioned scaled residual (\f$ \|r\|/\|b\|\f$) for this problem using the current solution.
  virtual double ScaledPrecResidual()= 0;
};

#endif /* _TSF_PRECLINEARPROBLEM_H_ */
