#ifndef _TSF_RDP_MULTIVECTOR_H_
#define _TSF_RDP_MULTIVECTOR_H_

//! TSF_RDP_Multivector:  The Trilinos Virtual Linear Problem Class.
/*! The TSF_RDP_Multivector class is a wrapper that encapsulates the 
  general information needed for solving a linear system of equations.  
*/

#include "TSF_RDP_LinearOperator.h"
#include "TSF_RDP_Vector.h"

class TSF_RDP_MultiVector {
    
  public:

  //! TSF_RDP_MultiVector Destructor.
  virtual ~TSF_RDP_Multivector(void){};

  //! Perform scaling of a multivector.
  /*! Applies the scaling vector D to the \e this multivector.
    \param In
           D - Vector containing scaling values.  D[i] will be applied 
               to the ith row of multivector.
    \return Integer error code, set to 0 if successful.
  */
  int Scale(const TSF_RDP_Vector & D) = 0;

};

#endif /* _TSF_RDP_MULTIVECTOR_H_ */
