#ifndef _SPOOLESOO_H_
#define _SPOOLESOO_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
//  #include "Epetra_LinearProblem.h"
#include "Epetra_Operator.h"
#include "Epetra_Object.h"


//! SpoolesOO:  An object-oriented wrapper for Spooles.
/*! 
  SpoolesOO will solve a linear systems of equations: \f$ AX=B \f$, using Epetra
  objects and the Spooles solver library, where \f$A\f$ is an Epetra_Operator or Epetra_RowMatrix (note
  that the Epetra_Operator class is a base class for Epetra_RowMatrix so that Epetra_RowMatrix \e isa
  Epetra_Operator.) \f$X\f$ and \f$B\f$ are Epetra_MultiVector objects.

  \warning SpoolesOO does not presently support solution of more than one simultaneous right-hand-side.
*/

class SpoolesOO {
    
  public:
  SpoolesOO(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);

  SpoolesOO();

  virtual ~SpoolesOO(void);
  
  int SetUserMatrix(Epetra_RowMatrix * UserMatrix);

  int SetLHS(Epetra_MultiVector * X);

  int SetRHS(Epetra_MultiVector * B);

  Epetra_RowMatrix * GetUserMatrix() const {return(UserMatrix_);};

  Epetra_MultiVector * GetLHS() const {return(X_);};

  Epetra_MultiVector * GetRHS() const {return(B_);};

  bool GetTrans( ) const { return Transpose_ ;} ;

  void SetTrans( bool trans ) { Transpose_ = trans ;} ; 

  int SetSpoolesDefaults();

  int Solve() ;

 protected:

  Epetra_Operator * UserOperator_;
  Epetra_RowMatrix * UserMatrix_;
  Epetra_Operator * PrecOperator_;
  Epetra_RowMatrix * PrecMatrix_;
  Epetra_MultiVector * X_;
  Epetra_MultiVector * B_;

  
  bool Transpose_ ;

  int x_LDA_;
  double *x_;
  int b_LDA_;
  double *b_;
  bool inConstructor_;
};


#endif /* _SPOOLESOO_H_ */

