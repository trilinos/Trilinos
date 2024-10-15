/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_SOLVEOBJECT_H
#define IFPACK_SOLVEOBJECT_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_CombineMode.h"
class Epetra_Flops;

//! Ifpack_OverlapSolveObject: Provides Overlapped Forward/back solve services for Ifpack.

class Ifpack_OverlapSolveObject: public virtual Epetra_Operator {

 public:
  //@{ \name Constructors/Destructor
  //! Constructor.
  Ifpack_OverlapSolveObject(char * Label, const Epetra_Comm & Comm);

  //! Copy constructor.
  Ifpack_OverlapSolveObject(const Ifpack_OverlapSolveObject & Source);

  //! Ifpack_OverlapSolveObject Destructor
  virtual ~Ifpack_OverlapSolveObject();
  //@}

  //@{ \name Initialization methods.

  //! Generate Ifpack_OverlapGraph object using current settings.
  void SetOverlapMode( Epetra_CombineMode OverlapMode) {OverlapMode_ = OverlapMode; return;}

  //! Define the operator to be used for the lower triangle
  int SetLowerOperator (Epetra_CrsMatrix * L, bool UseLTrans) {L_ = L; UseLTrans_ = UseLTrans; return(0);};

  //! Define the vector to be used for the diagonal
  int SetDiagonal (Epetra_Vector * D, bool UseDInv) {D_ = D; UseDInv_ = UseDInv; return(0);};

  //! Define the operator to be used for the upper triangle
  int SetUpperOperator (Epetra_CrsMatrix * U, bool UseUTrans) {U_ = U; UseUTrans_ = UseUTrans; return(0);};
  //@}

  //@{ \name Mathematical functions.
  
  
  //! Returns the result of a Ifpack_CrsIlut forward/back solve on a Epetra_MultiVector X in Y (works for Epetra_Vectors also).
  /*! 
    \param In
    Trans -If true, solve transpose problem.
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of multiplying U, D and L in that order on an Epetra_MultiVector X in Y.
  /*! 
    \param In
    Trans -If true, multiply by L^T, D and U^T in that order.
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Multiply(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the maximum over all the condition number estimate for each local ILU set of factors.
  /*! This functions computes a local condition number estimate on each processor and return the
      maximum over all processor of the estimate.
   \param In
    Trans -If true, solve transpose problem.
    \param Out
    ConditionNumberEstimate - The maximum across all processors of 
    the infinity-norm estimate of the condition number of the inverse of LDU.
  */
  int Condest(bool Trans, double & ConditionNumberEstimate) const;
  //@}

  //@{ \name Attribute access methods.

  //! Returns the overlap mode used to combine terms that are redundantly computed.
  /*! Since rows of the graph, and any related matrices are multiply owned, some values
      in the subdomain solves will be computed on multiple processors.  The overlap mode
      is used to determine how the redundant values that come in from other processors will
      be handled.
  */
  Epetra_CombineMode OverlapMode() const {return(OverlapMode_);}
  
  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {return(L().NumGlobalNonzeros()+U().NumGlobalNonzeros());};
  
  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {return(L().NumMyNonzeros()+U().NumMyNonzeros());};
  
  //! Returns the address of the L factor associated with this factored matrix.
  const Epetra_CrsMatrix & L() const {return(*L_);};
    
  //! Returns the address of the D factor associated with this factored matrix.
  const Epetra_Vector & D() const {return(*D_);};
    
  //! Returns the address of the L factor associated with this factored matrix.
  const Epetra_CrsMatrix & U() const {return(*U_);};
  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

    //! Returns a character string describing the operator
    const char * Label() const {return(Label_);};
    
    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! Note that this implementation of Apply does NOT perform a forward back solve with
        the LDU factorization.  Instead it applies these operators via multiplication with 
	U, D and L respectively.  The ApplyInverse() method performs a solve.

    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Multiply(Ifpack_OverlapSolveObject::UseTranspose(), X, Y));};

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! In this implementation, we use several existing attributes to determine how virtual
        method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
	the Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
	if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
	be accounted for when doing a triangular solve.

    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Solve(Ifpack_OverlapSolveObject::UseTranspose(), X, Y));};

    //! Returns 0.0 because this class cannot compute Inf-norm.
    double NormInf() const {return(0.0);};

    //! Returns false because this class cannot compute an Inf-norm.
    bool HasNormInf() const {return(false);};

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const {return(UseTranspose_);};

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {return(U_->OperatorDomainMap());};

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map & OperatorRangeMap() const{return(L_->OperatorRangeMap());};

    //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
    const Epetra_Comm & Comm() const{return(Comm_);};
  //@}

 protected:
    virtual int SetupXY(bool Trans, 
			const Epetra_MultiVector& Xin, const Epetra_MultiVector& Yin,
			Epetra_MultiVector * & Xout, Epetra_MultiVector * & Yout) const=0;
 private:
  char * Label_;
  Epetra_CrsMatrix * L_;
  bool UseLTrans_;
  Epetra_Vector * D_;
  bool UseDInv_;
  Epetra_CrsMatrix * U_;
  bool UseUTrans_;
  bool UseTranspose_;
  const Epetra_Comm & Comm_;
  mutable double Condest_;
  Epetra_Flops * Counter_;
  Epetra_CombineMode OverlapMode_;

};
#endif // IFPACK_SOLVEOBJECT_H
