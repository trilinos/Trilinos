// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_BORDEREDSOLVER_ABSTRACTSTRATEGY_H
#define LOCA_BORDEREDSOLVER_ABSTRACTSTRATEGY_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"           // for ReturnType, CopyType
#include "NOX_Abstract_MultiVector.H"     // for DenseMatrix

// forward declarations
namespace LOCA {
  namespace MultiContinuation {
    class ConstraintInterface;
    class ConstraintInterfaceMVDX;
  }
  namespace BorderedSolver {
    class AbstractOperator;
  }
}

namespace LOCA {

  //! %Strategies for solving bordered systems of equations
  namespace BorderedSolver {

    //! Abstract interface class for solving bordered sets of linear equations
    /*!
     * Abstract interface for solving systems of equations of the form
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
     * where \f$J\f$ is an \f$n\times n\f$ matrix, \f$A\f$ and \f$B\f$ are
     * \f$n\times m\f$, \f$C\f$ is \f$m\times m\f$, \f$X\f$ and \f$F\f$ are
     * \f$n\times p\f$ and \f$Y\f$ and \f$G\f$ are \f$m\times p\f$.  The action
     * of \f$J\f$ and its inverse are represnted by a
     * LOCA::BorderedSolver::AbstractOperator
     * while \f$A\f$ is a NOX::Abstract::MultiVector and \f$B\f$, \f$C\f$ are
     * represtend by the solution and parameter components of the derivative
     * of a constraint contained in
     * LOCA::MultiContinuation::ConstraintInterface.  All classes that
     * implement a method for
     * computing solutions to this system of equations should be derived from
     * this class.  Constructors for derived classes should be of the form:
     *
     * \code
     * class Derived : public AbstractStrategy {
     * public:
     *   Derived(
     *   const Teuchos::RCP<LOCA::GlobalData>& global_data,
     *   const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams);
     *   ...
     * };
     * \endcode
     *
     * where \c global_data is the LOCA global data object, \c topParams is
     * the parsed top-level parameter list, and \c solverParams
     * is a parameter list of bordered-solver parameters.
     *
     * This class and its children follow the Strategy pattern as defined
     * in Erich Gamma, et al. "Design Patterns:  Elements of Reusable
     * Object-Oriented Software." Addison Wesley, Boston, MA, 1995.
     */
    class AbstractStrategy {

    public:

      //! Constructor
      AbstractStrategy() {}

      //! Destructor
      virtual ~AbstractStrategy() {}

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
      const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC) = 0;

      //! Set blocks with multivector constraint
      /*!
       * This is a version of setMatrixBlocks that takes a multivector for
       * \c blockB.  This method has a default implementation to generate
       * a LOCA::MultiContinuation::MultiVecConstraint from \c blockB which
       * is then passed to the setMatrixBlocks() method.
       */
      virtual void setMatrixBlocksMultiVecConstraint(
      const Teuchos::RCP<const LOCA::BorderedSolver::AbstractOperator>& op,
      const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockA,
      const Teuchos::RCP<const NOX::Abstract::MultiVector>& blockB,
      const Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>& blockC);

      //! Intialize solver for a solve
      /*!
       * This should be called after setMatrixBlocks(), but before
       * applyInverse().
       */
      virtual NOX::Abstract::Group::ReturnType
      initForSolve() = 0;

      //! Intialize solver for a transpose solve
      /*!
       * This should be called after setMatrixBlocks(), but before
       * applyInverseTranspose().
       */
      virtual NOX::Abstract::Group::ReturnType
      initForTransposeSolve() = 0;

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
       *     \end{bmatrix}
       * \f]
       * where \f$U\f$ is \f$n\times p\f$, \f$V\f$ is \f$m\times p\f$ and the
       * other blocks are as defined above.
       */
      virtual NOX::Abstract::Group::ReturnType
      apply(const NOX::Abstract::MultiVector& X,
        const NOX::Abstract::MultiVector::DenseMatrix& Y,
        NOX::Abstract::MultiVector& U,
        NOX::Abstract::MultiVector::DenseMatrix& V) const = 0;

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
       *        A^T & C^T
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y
       *     \end{bmatrix}
       * \f]
       * where \f$U\f$ is \f$n\times p\f$, \f$V\f$ is \f$m\times p\f$ and the
       * other blocks are as defined above.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyTranspose(const NOX::Abstract::MultiVector& X,
             const NOX::Abstract::MultiVector::DenseMatrix& Y,
             NOX::Abstract::MultiVector& U,
             NOX::Abstract::MultiVector::DenseMatrix& V) const = 0;

      /*!
       * \brief Solves the extended system as defined above
       */
      /*!
       * The \em params argument is the linear solver parameters.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyInverse(Teuchos::ParameterList& params,
           const NOX::Abstract::MultiVector* F,
           const NOX::Abstract::MultiVector::DenseMatrix* G,
           NOX::Abstract::MultiVector& X,
           NOX::Abstract::MultiVector::DenseMatrix& Y) const = 0;

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
                NOX::Abstract::MultiVector::DenseMatrix& Y) const = 0;

    private:

      //! Private to prohibit copying
      AbstractStrategy(const AbstractStrategy&);

      //! Private to prohibit copying
      AbstractStrategy& operator = (const AbstractStrategy&);

    }; // Class AbstractStrategy

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_ABSTRACTSTRATEGY_H
