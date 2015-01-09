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

#ifndef LOCA_BORDEREDSOLVER_NESTED_H
#define LOCA_BORDEREDSOLVER_NESTED_H

#include "LOCA_BorderedSolver_AbstractStrategy.H"  // base class

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace BorderedSystem {
    class AbstractGroup;
  }
}

namespace LOCA {

  namespace BorderedSolver {

    //! Bordered system solver strategy for nested bordered systems
    /*!
     * This class implements a bordered solver strategy for the bordered
     * system
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
     * when \f$J\f$ itself has this block form.  It combines the blocks
     * for \f$A\f$, \f$B\f$, and \f$C\f$ and then instantiates a solver
     * as specified by the "Nested Bordered Solver" sublist of the
     * \c solverParams pass through the constructor.  This sublist should
     * specify the "Bordered Solver Method" for the solver as well as
     * any other parameters for that method, and any method that can
     * be instantiated through the LOCA::Factory is available.
     *
     * Note that the operator representing \f$J\f$ must implement the
     * LOCA::BorderedSolver::BorderedOperator interface, and the constraint
     * object representing \f$B\f$ must be of type
     * LOCA::MultiContinuation::ConstraintInterfaceMVDX.
     */
    class Nested :
      public LOCA::BorderedSolver::AbstractStrategy {

    public:

      //! Constructor.
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list
       * \param solverParams [in] Bordered solver parameters as described above
       */
      Nested(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

      //! Destructor
      virtual ~Nested();

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
       * \brief Solves the extended system as defined above using bordering
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
       * using bordering
       */
      /*!
       * The \em params argument is the linear solver parameters. If
       * \em isZeroF or \em isZeroG is true, than the corresponding
       * \em F or \em G pointers may be NULL.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyInverseTranspose(Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector* F,
                const NOX::Abstract::MultiVector::DenseMatrix* G,
                NOX::Abstract::MultiVector& X,
                NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    private:

      //! Private to prohibit copying
      Nested(const Nested&);

      //! Private to prohibit copying
      Nested& operator = (const Nested&);

    protected:

      //! Global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Solver parameters
      Teuchos::RCP<Teuchos::ParameterList> solverParams;

      //! Underlying solver
      Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> solver;

      //! Pointer to group storing J
      Teuchos::RCP<const LOCA::BorderedSystem::AbstractGroup> grp;

      //! Pointer to unbordered group
      Teuchos::RCP<const NOX::Abstract::Group> unbordered_grp;

      //! Width for bordered rows/columns
      int myWidth;

      //! Underling width
      int underlyingWidth;

      //! Number of my constraints
      int numConstraints;

    };
  } // namespace BorderedSolver
} // namespace LOCA

#endif
