// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef LOCA_BORDEREDSYSTEM_LAPACKDIRECTSOLVE_H
#define LOCA_BORDEREDSYSTEM_LAPACKDIRECTSOLVE_H

#include "LOCA_BorderedSystem_AbstractStrategy.H"  // base class
#include "NOX_LAPACK_Matrix.H"                     // class data element

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace MultiContinuation {
    class ConstraintInterfaceMVDX;
  }
  namespace LAPACK {
    class Group;
  }
}

namespace LOCA {

  namespace BorderedSystem {
    
    //! Bordered system solver strategy based on direct factorization
    /*! 
     * This class solves the extended system of equations
     * \f[
     *     \begin{bmatrix}
     *          J & A    \\
     *        B^T & C
     *     \end{bmatrix}
     *     \begin{bmatrix}
     *        X \\
     *        Y 
     *     \end{bmatrix} = 
     *     \begin{bmatrix}
     *        F \\
     *        G
     *     \end{bmatrix}
     * \f]
     * using a direct solve by augmenting the \f$J\f$ with \f$A\f$, \f$B\f$,
     * and \f$C\f$.  To support this, the group representing \f$J\f$ must
     * be a LOCA::LAPACK::Group and the contraint object representing
     * \f$B\f$ must be a LOCA::MultiContinuation::ConstraintInterfaceMVDX.
     */
    class LAPACKDirectSolve : public LOCA::BorderedSystem::AbstractStrategy {

    public:

      //! Constructor. 
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list 
       * \param solverParams [in] Bordered solver parameters.  Currently none
       * are referenced.
       */
      LAPACKDirectSolve(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& solverParams);

      //! Destructor
      virtual ~LAPACKDirectSolve();

      /*! 
       * \brief Set flag indicating whether F and A are stored in a 
       * continguous multivector.
       */
      virtual void setIsContiguous(bool flag);

      //! Set blocks
      /*!
       * The \c blockA or \c blockC pointer may be null if either is zero.
       * Whether block B is zero will be determined by querying \c blockB
       * via ConstraintInterface::isConstraintDerivativesXZero.
       */
      virtual void setMatrixBlocks(
	  const Teuchos::RefCountPtr<const NOX::Abstract::Group>& group,
	  const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>& blockA,
	  const Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
	  const Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix>& blockC);
  
      /*! 
       * \brief Computed extended matrix-multivector product
       */
      /*!
       * Computes
       * \f[
       *     \begin{bmatrix}
       *        U \\
       *        V
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *          J & A    \\
       *        B^T & C
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y 
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *          J*X + A*Y \\
       *        B^T*X + C*Y
       *     \end{bmatrix}.
       * \f]
       */
      virtual NOX::Abstract::Group::ReturnType 
      apply(const NOX::Abstract::MultiVector& X,
	    const NOX::Abstract::MultiVector::DenseMatrix& Y,
	    NOX::Abstract::MultiVector& U,
	    NOX::Abstract::MultiVector::DenseMatrix& V) const;

      /*! 
       * \brief Computed extended matrix transpose-multivector product
       */
      /*!
       * Computes
       * \f[
       *     \begin{bmatrix}
       *        U \\
       *        V
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *        J^T & B    \\
       *        A^T & C
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y 
       *     \end{bmatrix} = 
       *     \begin{bmatrix}
       *        J^T*X + B*Y \\
       *        A^T*X + C^T*Y
       *     \end{bmatrix}.
       * \f]
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyTranspose(const NOX::Abstract::MultiVector& X,
		     const NOX::Abstract::MultiVector::DenseMatrix& Y,
		     NOX::Abstract::MultiVector& U,
		     NOX::Abstract::MultiVector::DenseMatrix& V) const;

      /*! 
       * \brief Solves the extended system as defined above using bordering
       */
      /*!
       * The \em params argument is the linear solver parameters. If
       * \em isZeroF or \em isZeroG is true, than the corresponding
       * \em F or \em G pointers may be NULL.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyInverse(NOX::Parameter::List& params,
		   const NOX::Abstract::MultiVector* F,
		   const NOX::Abstract::MultiVector::DenseMatrix* G,
		   NOX::Abstract::MultiVector& X,
		   NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    private:

      //! Private to prohibit copying
      LAPACKDirectSolve(const LAPACKDirectSolve&);

      //! Private to prohibit copying
      LAPACKDirectSolve& operator = (const LAPACKDirectSolve&);

    protected:

      //! Global data object
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Solver parameters
      Teuchos::RefCountPtr<NOX::Parameter::List> solverParams;

      //! Pointer to group storing J
      Teuchos::RefCountPtr<const LOCA::LAPACK::Group> grp;

      //! Pointer to A block
      Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> A;

      //! Pointer to B block
      Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> B;

      //! Pointer to C block
      Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> C;

      //! The augmented Jacobian matrix
      Teuchos::RefCountPtr<NOX::LAPACK::Matrix> augmentedJ;

      //! Pivots used when factoring augmented Jacobian
      vector<int> pivots;

      //! Jacobian dimension
      int n;

      //! Number of additional rows/columns
      int m;

      //! Size of augmented Jacobian
      int N;

      //! flag indicating whether A block is zero
      bool isZeroA;

      //! flag indicating whether B block is zero
      bool isZeroB;

      //! flag indicating whether C block is zero
      bool isZeroC;

      //! flag indicating whether F block is zero
      bool isZeroF;

      //! flag indicating whether G block is zero
      bool isZeroG;

      //! flag indicating whether F and A are contiguous
      bool isContiguous;

    };
  } // namespace BorderedSystem
} // namespace LOCA

#endif
