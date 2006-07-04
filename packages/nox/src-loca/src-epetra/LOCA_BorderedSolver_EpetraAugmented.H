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

#ifndef LOCA_BORDEREDSOLVER_EPETRAAUGMENTED_H
#define LOCA_BORDEREDSOLVER_EPETRAAUGMENTED_H

#include "LOCA_BorderedSolver_AbstractStrategy.H"  // base class

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace MultiContinuation {
    class ConstraintInterfaceMVDX;
  }
  namespace Epetra {
    class Group;
  }
}

namespace LOCA {

  namespace BorderedSolver {
    
    /*! 
     * \brief Bordered system solver strategy based on augmenting the Jacobian
     * operator.
     */
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
     * by forming an augmented Epetra_Operator representing
     * \f[
     *     \begin{bmatrix}
     *          J & A    \\
     *        B^T & C
     *     \end{bmatrix}
     * \f]
     * by creating a new Epetra_Map for the additional equations.
     */
    class EpetraAugmented : public LOCA::BorderedSolver::AbstractStrategy {

    public:

      //! Constructor. 
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list 
       * \param solverParams [in] Bordered solver parameters.  Currently none
       * are referenced.
       */
      EpetraAugmented(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams);

      //! Destructor
      virtual ~EpetraAugmented();

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

      //! Intialize solver for a solve
      /*!
       * This should be called after setMatrixBlocks(), but before 
       * applyInverse().
       */
      virtual NOX::Abstract::Group::ReturnType 
      initForSolve();

      //! Intialize solver for a transpose solve
      /*!
       * This should be called after setMatrixBlocks(), but before 
       * applyInverseTranspose().
       */
      virtual NOX::Abstract::Group::ReturnType 
      initForTransposeSolve();
  
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
       * \brief Solves the extended system using the technique described
       * above.
       */
      /*!
       * The \em params argument is the linear solver parameters. If
       * \em isZeroF or \em isZeroG is true, than the corresponding
       * \em F or \em G pointers may be NULL.
       *
       * Note that if either the A or B blocks are zero, the system is
       * solved using a simple block elimination scheme instead of the
       * Householder scheme.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyInverse(Teuchos::ParameterList& params,
		   const NOX::Abstract::MultiVector* F,
		   const NOX::Abstract::MultiVector::DenseMatrix* G,
		   NOX::Abstract::MultiVector& X,
		   NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    private:

      //! Private to prohibit copying
      EpetraAugmented(const EpetraAugmented&);

      //! Private to prohibit copying
      EpetraAugmented& operator = (const EpetraAugmented&);

    protected:

      //! Global data object
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Solver parameters
      Teuchos::RefCountPtr<Teuchos::ParameterList> solverParams;

      //! Pointer to group storing J
      Teuchos::RefCountPtr<LOCA::Epetra::Group> grp;

      //! Pointer to A block
      Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> A;

      //! Pointer to B block
      Teuchos::RefCountPtr<const NOX::Abstract::MultiVector> B;

      //! Pointer to C block
      Teuchos::RefCountPtr<const NOX::Abstract::MultiVector::DenseMatrix> C;

      //! Pointer to constraint interface
      Teuchos::RefCountPtr<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> constraints;

      //! Number of constraint equations
      int numConstraints;

      //! flag indicating whether A block is zero
      bool isZeroA;

      //! flag indicating whether B block is zero
      bool isZeroB;

      //! flag indicating whether C block is zero
      bool isZeroC;

    };
  } // namespace BorderedSolver
} // namespace LOCA

#endif
