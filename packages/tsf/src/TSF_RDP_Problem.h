#ifndef _TSF_RDP_PROBLEM_H_
#define _TSF_RDP_PROBLEM_H_

//! TSF_RDP_Problem:  The Trilinos Virtual  Problem Class.
/*! The TSF_Problem class is a wrapper that encapsulates the 
  general information needed for solving a  system of equations.  
*/

#include "TSF_RDP_Operator.h"
#include "TSF_RDP_MultiVector.h"

class TSF_RDP_Problem {
    
  public:

  //! TSF_RDP_Problem Destructor.
  virtual ~TSF_RDP_Problem(void){};

  //! Set Operator \f$A\f$.
  int SetOperator(TSF_RDP_Operator & A) = 0;

  //! Get Operator \f$A\f$.
  int GetOperator(TSF_RDP_Operator * & A) = 0;

  //! Set Right Hand Side \f$B\f$.
  int SetRHS(TSF_RDP_MultiVector & B) = 0;

  //! Get Right Hand Side \f$B\f$.
  int GetRHS(TSF_RDP_MultiVector * & B) = 0;

  //! Set Initial guess/solution  \f$X\f$.
  int SetLHS(TSF_RDP_MultiVector & X) = 0;

  //! Get Initial guess/solution \f$B\f$.
  int GetLHS(TSF_RDP_MultiVector * & X) = 0;

  //! Check input parameters, print summary if output level > 1.
    /*! Returns 0 if all input parameters are valid */
    int CheckInput() = 0;

  //! Returns the true unscaled residual(\f$ \|r\| = \|b-Ax\|\f$) for this problem using the current solution.
  double TrueResidual() = 0;

  //! Returns the true scaled residual (\f$ \|r\|/\|b\|\f$) for this problem using the current solution.
  double ScaledResidual()= 0;

  //! Returns the scaled residual tolerance \f$\gamma\f$ such that we want (\f$ \|r\|/\|b\|<\gamma\f$).
  double ScaledResidualTolerance()= 0;

  //! Determines if the problem has converged.
  /*! This function is for the purposes of allowing the solver to query whether or not the problem is converged.
      In should be implemented to test whatever stopping criteria the TSF_Problem wants to use to decide
      if convergence was met.  The solver should provide its current residual estimate as an input to this method.
      If no estimate is available, SolverResidual should be set to 0.0.

      This method removes the responsibility for determining convergence from the solver and places it on the 
       problem.
  */
  bool Converged(double SolverResidual) = 0;
};

#endif /* _TSF_PROBLEM_H_ */
