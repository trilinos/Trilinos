/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#ifndef IFPACK_CRSILUT_H
#define IFPACK_CRSILUT_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_ScalingType.h"
#include "Ifpack_OverlapGraph.h"
#include "Ifpack_OverlapFactorObject.h"
#include "Ifpack_OverlapSolveObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Object.h"
class Epetra_Comm;
class Epetra_Map;
class Epetra_RowMatrix;
class Epetra_Vector;

#ifdef HAVE_IFPACK_TEUCHOS
namespace Teuchos {
  class ParameterList;
}
#endif

//! Ifpack_Jacobi: Jacobi preconditioner of a given Epetra_RowMatrix.
/*! This class supports the construction and use of Jacobi's basic iterative method as a 
    preconditioner for a Krylov iterative method.  It is also possible to use this class to solve a 
    problem using Jacobi's method only.  Formally Jacobi's method is an iteration of the form:

    \f[ x_{k+1} = D^{-1}(E+F)x_k + D^{-1}b \f]

    where \f$E\f$, \f$D\f$ and \f$F\f$ are the strictly lower triangle, diagonal, and upper triangular
    parts, resp. of the user matrix \f$A\f$. (See Saad \em {Iterative Methods for Sparse Linear Systems}, Ch. 4).
    To start the Jacobi iteration, we use an initial guess of \f$ x_0 = 0\f$, so a single step is equivalent to scaling
    the input vector by the inverse of the diagonal of \f$A\f$.

    Use of more than one step is often not beneficial.

    There are two parameters for this class:
<ol> 
<li> UseReciprocal - Specifies whether or not the explicit reciprocal of the diagonal elements should
     computed and used.  Jacobi's method requires division using the diagonal of the user matrix.  On
     many processors, multiplication by the reciprocal is much faster than division.  Thus, it is
     advantageous to explicitly compute the inverse of the diagonal terms of the user matrix and multiply
     by the reciprocal on each iteration.  At the same time, multiplying by the reciprocal is numerically
     less accurate than division, so in some situations it may be advantageous to set this parameter to false.
<li> NumSteps - Jacobi's method is an iterative method in its own right.  Typically, when used to precondition a
     Krylov method, only a single step (iteration) of Jacobi is appropriate.  More than one iteration can significantly
     increase cost and even cause the Krylov method to fail.  

*/

class Ifpack_Jacobi: public Epetra_Object, public Epetra_CompObject, public Ifpack_OverlapFactorObject, public Ifpack_OverlapSolveObject {
  
 public:
  //@{ \name Constructors/Destructor

  //! Constructor using Ifpack_OverlapGraph.
  /*! Creates an object from the overlap graph. 
    \param OverlapGraph (In) - Graph describing the graph that should be used for the factors.
    \param UseReciprocal (In/Default) - Explicitly compute and use the reciprocal of diagonal elements if set to true.
    \param NumSteps (In/Default) - Specify the number of Jacobi steps to perform on each call to the Solve() method.

  */
  Ifpack_Jacobi(const Ifpack_OverlapGraph * OverlapGraph, bool UseReciprocal = true, int NumSteps = 1);

  //! Constructor using Epetra_RowMatrix.
  /*! Creates an Ifpack_Graph object from the user graph implicitly defined by the
	 Epetra_RowMatrix interface. 
    \param RowMatrix (In) - An object that has implemented the Epetra_RowMatrix interface.
    \param UseReciprocal (In/Default) - Explicitly compute and use the reciprocal of diagonal elements if set to true.
    \param NumSteps (In/Default) - Specify the number of Jacobi steps to perform on each call to the Solve() method.

  */
  Ifpack_Jacobi(const Epetra_RowMatrix * UserMatrix, bool UseReciprocal = true, int NumSteps = 1);
  
  //! Copy constructor.
  Ifpack_Jacobi(const Ifpack_Jacobi & Source);

  //! Ifpack_Jacobi Destructor
  virtual ~Ifpack_Jacobi();
  //@}

  //@{ \name Attribute access methods.

#ifdef HAVE_IFPACK_TEUCHOS
  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the configure argument
     '--enable-ifpack-teuchos' was used.
     This method recognizes two parameter names: use_reciprocal and num_steps.
     These names are case insensitive. For each, the ParameterEntry must have
     type int.
  */
  int SetParameters(const Teuchos::ParameterList& parameterlist,
                    bool cerr_warning_if_unused=false);
#endif

  //! Returns current value of UseReciprocal.
  bool UseReciprocal() const {return(UseReciprocal_);};

  //! Returns current value of NumSteps.
  int NumSteps() const {return(NumSteps_);};

  //! Returns current vector of diagonal values.
  const Epetra_Vector & DiagValues() const {return(*DiagValues_);};
  //@}
  
 protected:
  //@{ \name Methods needed to implement OverlapFactorObject.

  //! Processes the overlapped user matrix for computing the ILUT preconditioner: WARNING: THIS ROUTINE IS NOT USER CALLABLE, CALL InitValues().
  int ProcessOverlapMatrix(const Epetra_RowMatrix &A);
  //! Compute ILUT factors L and U: WARNING: THIS ROUTINE IS NOT USER CALLABLE, CALL Factor().
  int DerivedFactor();
  //@}

 private:

 bool UseReciprocal_;
 int NumSteps_;
 Epetra_Vector * DiagValues_;

};

#endif /* IFPACK_CRSILUT_H */
