// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef _KUNDERTOO_H_
#define _KUNDERTOO_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
//! KundertOO:  An object-oriented wrapper for Kundert.
/*! 
  KundertOO will solve a linear systems of equations: \f$ AX=B \f$, using Epetra
  objects and the Kundert solver library, where \f$A\f$ is an Epetra_Operator or Epetra_RowMatrix (note
  that the Epetra_Operator class is a base class for Epetra_RowMatrix so that Epetra_RowMatrix \e isa
  Epetra_Operator.) \f$X\f$ and \f$B\f$ are Epetra_MultiVector objects.

  \warning KundertOO does not presently support solution of more than one simultaneous right-hand-side.
*/

class KundertOO {
    
  public:
  //@{ \name Constructor methods
  //! Kundert Constructor.
  /*! Creates an Kundert instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object.
  */
  KundertOO(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);

  //! Kundert Destructor.
  /*! Completely deletes a Kundert object.  
  */
  KundertOO();
  //@}

  virtual ~KundertOO(void);
  
  bool GetTrans( ) const { return Transpose_ ;} ;

  void SetTrans( bool trans ) { Transpose_ = trans ;} ; 

  //@{ \name Solve method
  //!  All computation is performed during the call to Solve() 
  /*!  Factor controls whether or not the matrix should be factored prior to the solve.
       Default is true.
   */
  int Solve() ;

 protected:
  int SetUserMatrix(Epetra_RowMatrix * UserMatrix);

  int SetLHS(Epetra_MultiVector * X);

  int SetRHS(Epetra_MultiVector * B);

  Epetra_RowMatrix * GetUserMatrix() const {return(UserMatrix_);};

  Epetra_MultiVector * GetLHS() const {return(X_);};

  Epetra_MultiVector * GetRHS() const {return(B_);};

  int SetKundertDefaults();


  Epetra_Operator * UserOperator_;
  Epetra_RowMatrix * UserMatrix_;
  //  Epetra_Operator * PrecOperator_;
  //  Epetra_RowMatrix * PrecMatrix_;
  Epetra_MultiVector * X_;
  Epetra_MultiVector * B_;

  bool Transpose_ ;

  int x_LDA_;
  double *x_;
  int b_LDA_;
  double *b_;
  bool inConstructor_;
};


#endif /* _KUNDERTOO_H_ */

