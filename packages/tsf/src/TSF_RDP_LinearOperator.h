#ifndef _TSF_RDP_LINEAROPERATOR_H_
#define _TSF_RDP_LINEAROPERATOR_H_

//! TSF_RDP_LinearOperator:  The Trilinos Real Double Precision Linear Operator Class.
/*! The TSF_RDP_LinearOperator class is a pure virtual class that specifies
    the required interfaces that any Trilinos-compliant real double precision linear
    operator must implement.  It extends the TSF_RDP_Operator class.
*/
#include "TSF_Parameter.h"
#include "TSF_RDP_Operator.h"
#include "TSF_RDP_Vector.h"

class TSF_RDP_LinearOperator: public virtual TSF_RDP_Operator {
    
  public:

  //! TSF_RDP_LinearOperator Destructor.
  virtual ~TSF_RDP_LinearOperator(void){};

  //! Perform left scaling of a linear operator.
  /*! Applies the scaling vector D to the left side of the Linear Operator.  If the operator is not available
      in matrix form, then this vector must be applied each time the Apply method is called for this operator.
    \param In
           D - Vector containing scaling values.  D[i] will be applied 
               to the ith row of the operator.
    \return Integer error code, set to 0 if successful.
  */
  virtual int LeftScale(const TSF_RDP_Vector & D) = 0;

  //! Perform right scaling of a linear operator.
  /*! Applies the scaling vector D to the right side of the Linear Operator.  If the operator is not available
      in matrix form, then this vector must be applied each time the Apply method is called for this operator.
    \param In
           D - Vector containing scaling values.  D[i] will be applied 
               to the jth column of operator.
    \return Integer error code, set to 0 if successful.
  */
  virtual int RightScale(const TSF_RDP_Vector & D) = 0;


};

#endif /* _TSF_RDP_LINEAROPERATOR_H_ */
