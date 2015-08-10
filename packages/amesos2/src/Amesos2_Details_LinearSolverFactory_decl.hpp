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

/// \file   Amesos2_Details_LinearSolverFactory_decl.hpp
/// \author Mark Hoemmen
/// \brief  Declaration of Amesos2::Details::LinearSolverFactory.

#ifndef AMESOS2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP
#define AMESOS2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP

#include "Amesos2_config.h"
#include "Trilinos_Details_LinearSolverFactory.hpp"

namespace Amesos2 {
namespace Details {

  /// \class LinearSolverFactory
  /// \brief Interface for a "factory" that creates Amesos2 solvers.
  ///
  /// \tparam MV Type of a (multi)vector, representing either the
  ///   solution(s) X or the right-hand side(s) B of a linear system
  ///   AX=B.  For example, with Tpetra, use a Tpetra::MultiVector
  ///   specialization.  A <i>multivector</i> is a single data structure
  ///   containing zero or more vectors with the same dimensions and
  ///   layout.
  ///
  /// \tparam OP Type of a matrix or linear operator that this Solver
  ///   understands.  For example, for Tpetra, use a Tpetra::Operator
  ///   specialization.  Always use the most abstract interface
  ///   possible; solvers should dynamic_cast to the subclass they
  ///   need.  Also, be consistent: using different classes here
  ///   (e.g., Tpetra::RowMatrix instead of Tpetra::Operator) means
  ///   more expensive explicit template instantiation.
  ///
  /// \tparam NormType Type of the norm of a residual.
  template<class MV, class OP, class NormType>
  class LinearSolverFactory :
    public Trilinos::Details::LinearSolverFactory<MV, OP, NormType> {
  public:
    /// \brief Get an instance of a Amesos2 solver.
    ///
    /// The solver is wrapped in a Trilinos::Details::LinearSolver
    /// interface.
    ///
    /// \param solverName [in] The solver's name.  Names are case
    ///   sensitive.
    /// \return A pointer to the solver, if the name was valid; else,
    ///   a null pointer (Teuchos::null).
    virtual Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, NormType> >
    getLinearSolver (const std::string& solverName);
  };

} // namespace Details
} // namespace Amesos2

#endif // AMESOS2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP
