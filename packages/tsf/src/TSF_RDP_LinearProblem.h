#ifndef _TSF_RDP_LINEARPROBLEM_H_
#define _TSF_RDP_LINEARPROBLEM_H_

//! TSF_LinearProblem:  The Trilinos Virtual Linear Problem Class.
/*! The TSF_LinearProblem class is a wrapper that encapsulates the 
  general information needed for solving a linear system of equations.  
*/

#include "TSF_RDP_LinearOperator.h"
#include "TSF_RDP_Vector.h"
#include "TSF_RDP_Problem.h"
class TSF_RDP_LinearProblem: public virtual TSF_RDP_Problem {
    
  public:

  //! TSF_RDP_LinearProblem Destructor.
  virtual ~TSF_RDP_LinearProblem(void){};

  //! Perform left scaling of a linear problem.
  /*! Applies the scaling vector D to the left side of the operator A defined by SetOperator() and
      to the right hand side B defined by SetRHS(). If A is not a matrix,
      D must be applied implicitly.
    \param In
           D - Vector containing scaling values.  D[i] will be applied 
               to the ith row of the operator A and RHS B.
    \return Integer error code, set to 0 if successful.
  */
  virtual int LeftScale(const TSF_RDP_Vector & D) = 0;

  //! Perform right scaling of a linear problem.
  /*! Applies the scaling vector D to the right side of the operator A defined by SetOperator().
      Apply the inverse of D to the left hand side X defined by SetLHS().  If A is not a matrix,
      D must be applied implicitly.
    \param In
           D - Vector containing scaling values.  D[j] will be applied 
               to the jth column of A.  1/D[j] will be applied to the
               jth row of B.

    \return Integer error code, set to 0 if successful.
  */
  virtual int RightScale(const TSF_RDP_Vector & D) = 0;

  //! Apply operator that was specified in SetOperator() to a given TSF_RDP_Multivector X, put results in Y.
  /*! Assuming OperatorSetup() is called, applies the operator.
    \param In 
           X - User supplied input multivector.
    \param Out 
           Y - Output multivector.

    \return Integer error code, set to 0 if successful.

  */
  virtual int ApplyOperator(TSF_RDP_MultiVector & X, TSF_RDP_MultiVector & Y) = 0;

};

#endif /* _TSF_RDP_LINEARPROBLEM_H_ */
