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

#ifndef LOCA_BORDEREDSOLVER_UPPERTRIANGULARBLOCKELIMINATION_H
#define LOCA_BORDEREDSOLVER_UPPERTRIANGULARBLOCKELIMINATION_H

#include "LOCA_BorderedSolver_AbstractStrategy.H"  // base class

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace BorderedSolver {

    //! Block elimination strategy for solving a block upper-triangular system
    /*!
     * This class solves the extended system of equations
     * \f[
     *     \begin{bmatrix}
     *        op(J) &   A    \\
     *        0     & op(C)
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
     * via block elimination:
     * \f[
     *     \begin{aligned}
     *        Y &= op(C)^{-1}G \\
     *        X &= op(J)^{-1}(F - A Y)
     *     \end{aligned}
     * \f]
     * where \f$op\f$ represents either the identity operation or the
     * transpose.  \f$C\f$ must be nonzero, while \f$A\f$, \f$F\f$ or \f$G\f$
     * may be zero.  The solve
     * for the non-transposed system is implemented by the solve() method,
     * while the solve for the transposed system is implemented by the
     * solveTranspose() method.
    */
    class UpperTriangularBlockElimination {

    public:

      //! Constructor.
      /*!
       * \param global_data [in] Global data object
       */
      UpperTriangularBlockElimination(
     const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~UpperTriangularBlockElimination();

      /*!
       * \brief Solves the extended system as described above.
       */
      /*!
       * Either A, F, or G may be zero by passing NULL.
       */
      NOX::Abstract::Group::ReturnType
      solve(Teuchos::ParameterList& params,
        const LOCA::BorderedSolver::AbstractOperator& op,
        const NOX::Abstract::MultiVector* A,
        const NOX::Abstract::MultiVector::DenseMatrix& C,
        const NOX::Abstract::MultiVector* F,
        const NOX::Abstract::MultiVector::DenseMatrix* G,
        NOX::Abstract::MultiVector& X,
        NOX::Abstract::MultiVector::DenseMatrix& Y) const;

      /*!
       * \brief Solves the extended system using the tranpose of J and C
       * as described above.
       */
      /*!
       * Either A, F, or G may be zero by passing NULL.
       */
      NOX::Abstract::Group::ReturnType
      solveTranspose(Teuchos::ParameterList& params,
             const LOCA::BorderedSolver::AbstractOperator& op,
             const NOX::Abstract::MultiVector* A,
             const NOX::Abstract::MultiVector::DenseMatrix& C,
             const NOX::Abstract::MultiVector* F,
             const NOX::Abstract::MultiVector::DenseMatrix* G,
             NOX::Abstract::MultiVector& X,
             NOX::Abstract::MultiVector::DenseMatrix& Y) const;

    private:

      //! Private to prohibit copying
      UpperTriangularBlockElimination(const UpperTriangularBlockElimination&);

      //! Private to prohibit copying
      UpperTriangularBlockElimination&
      operator = (const UpperTriangularBlockElimination&);

    protected:

      //! Global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // class UpperTriangularBlockElimination
  } // namespace BorderedSolver
} // namespace LOCA

#endif
