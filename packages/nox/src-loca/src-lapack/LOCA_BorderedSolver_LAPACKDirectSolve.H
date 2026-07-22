// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_BORDEREDSOLVER_LAPACKDIRECTSOLVE_H
#define LOCA_BORDEREDSOLVER_LAPACKDIRECTSOLVE_H

#include "LOCA_BorderedSolver_AbstractStrategy.H"  // base class
#include "NOX_LAPACK_LinearSolver.H"                     // class data element

#include "Teuchos_ConfigDefs.hpp"              // for complex<double>

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

  namespace BorderedSolver {

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
    class LAPACKDirectSolve : public LOCA::BorderedSolver::AbstractStrategy {

    public:

      //! Constructor.
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list
       * \param solverParams [in] Bordered solver parameters.  Currently none
       * are referenced.
       */
      LAPACKDirectSolve(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

      //! Destructor
      virtual ~LAPACKDirectSolve();

      //! Set blocks
      /*!
       * The \c blockA or \c blockC pointer may be null if either is zero.
       * Whether block B is zero will be determined by querying \c blockB
       * via ConstraintInterface::isConstraintDerivativesXZero.
       */
      virtual void setMatrixBlocks(
          const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& op,
      const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
      const Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterface>& blockB,
      const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC);

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
       * \brief Solves the extended system as defined above using
       * Gaussian Elimination
       */
      /*!
       * The \em params argument is the linear solver parameters. If
       * \em isZeroF or \em isZeroG is true, than the corresponding
       * \em F or \em G pointers may be NULL.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyInverse(Teuchos::ParameterList& params,
           const NOX::Abstract::MultiVector* F,
           const NOX::Abstract::MultiVector::DenseMatrix* G,
           NOX::Abstract::MultiVector& X,
           NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*!
       * \brief Solves the transpose of the extended system as defined above
       */
      /*!
       * The \em params argument is the linear solver parameters.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyInverseTranspose(Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector* F,
                const NOX::Abstract::MultiVector::DenseMatrix* G,
                NOX::Abstract::MultiVector& X,
                NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    protected:

      //! Solve linear system
      NOX::Abstract::Group::ReturnType
      solve(bool trans,
        Teuchos::ParameterList& params,
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
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Solver parameters
      Teuchos::RCP<Teuchos::ParameterList> solverParams;

      //! Pointer to group
      Teuchos::RCP<const LOCA::LAPACK::Group> grp;

      //! Pointer to operator
      Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator> op;

      //! Pointer to A block
      Teuchos::RCP<const NOX::Abstract::MultiVector> A;

      //! Pointer to B block
      Teuchos::RCP<const LOCA::MultiContinuation::ConstraintInterfaceMVDX> B;

      //! Pointer to C block
      Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix> C;

      //! The augmented Jacobian matrix solver
      Teuchos::RCP< NOX::LAPACK::LinearSolver<double> > augJacSolver;

#ifdef HAVE_TEUCHOS_COMPLEX
      //! The augmented complex matrix solver
      Teuchos::RCP< NOX::LAPACK::LinearSolver< std::complex<double> > > augComplexSolver;
#endif

      //! Matrix dimension
      int n;

      //! Number of additional rows/columns
      int m;

      //! Size of augmented matrix
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

      //! flag indicating whether we are solving the complex matrix or not
      bool isComplex;

    };
  } // namespace BorderedSolver
} // namespace LOCA

#endif
