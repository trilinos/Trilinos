#ifndef _TSF_RDP_PRECONDITIONER_H_
#define _TSF_RDP_PRECONDITIONER_H_

//! TSF_RDP_Preconditioner:  The Trilinos Real Double Precision Operator Class.
/*! The TSF_RDP_Preconditioner class is a pure virtual class that specifies
    the required interfaces that any Trilinos-compliant real double precision 
    operator must implement.
*/
#include "TSF_RDP_Parameter.h"
#include "TSF_RDP_Operator.h"

class TSF_RDP_Preconditioner: public virtual TSF_RDP_Operator {
    
  public:
  //! TSF_RDP_Preconditioner Destructor.
  virtual ~TSF_RDP_Preconditioner(void){};

  //! Performs any required setup that is needed by Apply(), must be called after SetParameters().
  virtual int PreconditionerSetup(void) = 0;

  //! Returns true if there is a non-trivial left preconditioner.
  virtual bool LeftPreconditioner() = 0;

  //! Returns true if there is a non-trivial right preconditioner.
  virtual bool RightPreconditioner() = 0;

  //! Apply left side of two-sided preconditioning to a given TSF_Multivector X, put results in Y.
  /*! Assuming PreconditionerSetup is called, applies the left side member of a two-sided preconditioner.
    \param In 
           X - User supplied input multivector.
    \param Out 
           Y - Output multivector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int SolveLeft(TSF_RDP_MultiVector & X, TSF_RDP_MultiVector & Y) = 0;

  //! Apply right side of two-sided preconditioning to a given TSF_Multivector X, put results in Y.
  /*! Assuming PreconditionerSetup is called, applies the right side member of a two-sided preconditioner.
    \param In 
           X - User supplied input multivector.
    \param Out 
           Y - Output multivector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int SolveRight(TSF_RDP_MultiVector & X, TSF_RDP_MultiVector & Y) = 0;

  //! Apply preconditioning to a given TSF_Multivector X, put results in Y.
  /*! Assuming PreconditionerSetup is called, applies the complete preconditioner to X.
    \param In 
           X - User supplied input multivector.
    \param Out 
           Y - Output multivector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Solve(TSF_RDP_MultiVector & X, TSF_RDP_MultiVector & Y) = 0;

};

#endif /* _TSF_RDP_PRECONDITIONER_H_ */
