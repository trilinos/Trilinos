// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
#ifndef LOCA_EPETRA_COMPACTWYOP_H
#define LOCA_EPETRA_COMPACTWYOP_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"     // class data element

#include "Epetra_Operator.h"	// base class
#include "Epetra_LocalMap.h"    // class data element
#include "Epetra_MultiVector.h" // class data element
#include "NOX_Common.H"         // for std::string data member
#include "NOX_Abstract_MultiVector.H" // for DenseMatrix

// forward declarations
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

  namespace Epetra {
    
    /*! 
     * \brief An Epetra operator for solving extended sets of equations using
     * Householder transformations.
     */
    /*!
     * This class implements the \f$P\f$ operator as described in the
     * LOCA::BorderedSolver::EpetraHouseholder documentation for solving
     * an extended set of equations.  It uses the \f$Q\f$ factor from 
     * a QR factorization using the compact WY representation.
     */
    class CompactWYOp : public Epetra_Operator {
      
    public:

      //! Constructor 
      /*!
       * \param global_data [in] The global data object
       * \param jacOperator [in] Jacobian operator J
       * \param A_multiVec [in] Multivector representing A
       * \param Y_x_multiVec [in] Multivector representing the solution 
       *  component of the Y matrix in the compact WY representation
       * \param Y_p_matrix [in] Matrix representing the parameter 
       *  component of the Y matrix in the compact WY representation
       * \param T_matrix [in] Matrix representing the T matrix in the 
       *  compact WY representation.
       */
      CompactWYOp(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<const Epetra_Operator>& jacOperator, 
	const Teuchos::RCP<const Epetra_MultiVector>& A_multiVec, 
	const Teuchos::RCP<const Epetra_MultiVector>& Y_x_multiVec,
	const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& Y_p_matrix,
	const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& T_matrix);

      //! Destructor
      virtual ~CompactWYOp();

      //! The operator currently does not support a transpose
      /*! 
       * Setting this to true throws an error.
       */
      virtual int SetUseTranspose(bool UseTranspose);

      /*! 
       * \brief Returns the result of a Epetra_Operator applied to a 
       * Epetra_MultiVector Input in Result as described above.
       */
      virtual int Apply(const Epetra_MultiVector& Input, 
			Epetra_MultiVector& Result) const;

      /*!
       * \brief This method does nothing.
       */
      virtual int ApplyInverse(const Epetra_MultiVector& X, 
			       Epetra_MultiVector& Y) const;
  
      //! Returns an approximate infinity norm of the operator matrix.
      /*!
       * This is defined only if %NormInf() of the underlying operator \f$J\f$
       * is defined and is given by 
       * \f$\|J\|_\infty+\|A\|_\infty\f$.
       */
      virtual double NormInf() const;
    
      //! Returns a character std::string describing the operator
      virtual const char* Label () const;
  
      //! Returns the current UseTranspose setting.  Always returns false.
      virtual bool UseTranspose() const;

      /*! 
       * \brief Returns true if the \e this object can provide an 
       * approximate Inf-norm, false otherwise.
       */
      virtual bool HasNormInf() const;

      /*! 
       * \brief Returns a reference to the Epetra_Comm communicator 
       * associated with this operator.
       */
      virtual const Epetra_Comm & Comm() const;

      /*!
       * \brief Returns the Epetra_Map object associated with the 
       * domain of this matrix operator.
       */
      virtual const Epetra_Map& OperatorDomainMap () const;

      /*! 
       * \brief Returns the Epetra_Map object associated with the 
       * range of this matrix operator.
       */
      virtual const Epetra_Map& OperatorRangeMap () const;

      /*!
       * \brief Initialize operator.  Call this before starting a linear
       * solve.  The Epetra_MultiVector argument \em x must be of the
       * same size and distribution as arguments to Apply().
       */
      void init(const Epetra_MultiVector& x);

      /*!
       * \brief Finish up solve.  Call this after a linear solve is
       * finished to inform the operator that the solve is completed.
       */
      void finish();

      /*!
       * \brief Applies the operator Q with a zero parameter component
       * on input.
       */
      void applyCompactWY(const Epetra_MultiVector& x, 
			  Epetra_MultiVector& result_x, 
			  Epetra_MultiVector& result_p) const;

    private:

      //! Private to prohibit copying
      CompactWYOp(const CompactWYOp&);

      //! Private to prohibit copying
      CompactWYOp& operator=(const CompactWYOp&);

    protected:

      //! Global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Label for operator
      std::string label;

      //! Local map for generating temporary matrices
      Epetra_LocalMap localMap;
    
      //! Stores operator representing J
      Teuchos::RCP<const Epetra_Operator> J;
  
      //! Stores multivector representing A
      Teuchos::RCP<const Epetra_MultiVector> A;

      //! Stores multivector representing solution component of Y
      Teuchos::RCP<const Epetra_MultiVector> Y_x;

      //! Stores multivector representing parameter component of Y
      Epetra_MultiVector Y_p;

      //! Stores multivector representing T
      Epetra_MultiVector T;

      //! Temporary matrix
      Epetra_MultiVector* tmpMat1;

      //! Temporary multivec
      Epetra_MultiVector* tmpMV;

      //! BLAS wrappers, used for triangular matrix-matrix products
      Teuchos::BLAS<int,double> dblas;

    };
  } // namespace Epetra
} // namespace LOCA

#endif 
