/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * October 20, 2002, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

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
