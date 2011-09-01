// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//
// @HEADER

/**
  \file   Amesos2_Status.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu Jan 14 08:52:04 2010

  \brief  Container class for status variables.
*/

#ifndef AMESOS2_STATUS_HPP
#define AMESOS2_STATUS_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

#include "Amesos2_TypeDecl.hpp"

namespace Amesos2 {

  template < template <class,class> class ConcreteSolver, class Matrix, class Vector > class SolverCore;

  /**
   * \brief Holds internal status data about the owning Amesos2 solver
   *
   * Provides convenient access to status data for both solvers and
   * outside users.
   *
   * \sa Solver::getStatus()
   */
  class Status {
  public:
    template < template <class,class> class ConcreteSolver, class Matrix, class Vector >
    friend class SolverCore;

    Status()
      : numPreOrder_(0)
      , numSymbolicFact_(0)
      , numNumericFact_(0)
      , numSolve_(0)

      , last_phase_(CLEAN)

      , lu_nnz_(0)
    { }


    /// Default destructor.
    ~Status() { };

    /// Returns the number of pre-orderings performed by the owning solver.
    inline int getNumPreOrder() const
    { return( numPreOrder_ ); }

    /// Returns the number of symbolic factorizations performed by the owning solver.
    inline int getNumSymbolicFact() const
    { return( numSymbolicFact_ ); }

    /// Returns the number of numeric factorizations performed by the owning solver.
    inline int getNumNumericFact() const
    { return( numNumericFact_ ); }

    /// Returns the number of solves performed by the owning solver.
    inline int getNumSolve() const
    { return( numSolve_ ); }

    /// If \c true , then pre-ordering has been performed
    inline bool preOrderingDone() const
    { return( last_phase_ >= PREORDERING ); }

    /// If \c true , then symbolic factorization has been performed
    inline bool symbolicFactorizationDone() const
    { return( last_phase_ >= SYMBFACT ); }

    /// If \c true , then numeric factorization has been performed
    inline bool numericFactorizationDone() const
    { return( last_phase_ >= NUMFACT ); }

    /** \brief Get the number of non-zero entries in the \f$L\f$ and \f$U\f$ factors.
     * 
     * Returns \c 0 if numeric factorization has not yet been performed,
     * or if the solver does not support getting statistics about the
     * factors.
     */
    inline size_t getNnzLU() const
    { return( lu_nnz_ ); }


  private:

    /// Number of pre-ordering phases
    int numPreOrder_;

    /// Number of symbolic factorization phases.
    int numSymbolicFact_;

    /// Number of numeric factorization phases.
    int numNumericFact_;

    /// Number of solves.
    int numSolve_;

    /// The last phase of computation that was performed by the owning solver object
    EPhase last_phase_;

    /// The number of non-zeros in the factors
    size_t lu_nnz_;

  };                              // end class Amesos2::Status


} // end namespace Amesos2

#endif  // AMESOS2_STATUS_HPP
