/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

/// \file   Ifpack2_Details_LinearSolverFactory_decl.hpp
/// \author Mark Hoemmen
/// \brief  Declaration of Ifpack2::Details::LinearSolverFactory.

#ifndef IFPACK2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP
#define IFPACK2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
#include "Tpetra_Operator.hpp"

namespace Ifpack2 {
namespace Details {

/// \class LinearSolverFactory
/// \brief Interface for a "factory" that creates Ifpack2 solvers.
///
/// We use Tpetra's template parameters here, instead of MV, OP, and
/// NormType, because this is not a public-facing class.  We also want
/// to avoid mix-ups between MV and OP.
template<class SC, class LO, class GO, class NT>
class LinearSolverFactory :
  public Trilinos::Details::LinearSolverFactory<Tpetra::MultiVector<SC, LO, GO, NT>,
                                                Tpetra::Operator<SC, LO, GO, NT>,
                                                typename Tpetra::MultiVector<SC, LO, GO, NT>::mag_type>
{
public:
  typedef Trilinos::Details::LinearSolver<Tpetra::MultiVector<SC, LO, GO, NT>,
                                          Tpetra::Operator<SC, LO, GO, NT>,
                                          typename Tpetra::MultiVector<SC, LO, GO, NT>::mag_type> solver_type;

  /// \brief Get an instance of a Ifpack2 solver.
  ///
  /// The solver is wrapped in a Trilinos::Details::LinearSolver
  /// interface.
  ///
  /// \param solverName [in] The solver's name.  Names are case
  ///   sensitive.
  /// \return A pointer to the solver, if the name was valid; else,
  ///   a null pointer (Teuchos::null).
  virtual Teuchos::RCP<solver_type>
  getLinearSolver (const std::string& solverName);

  /// \brief Register this LinearSolverFactory with the central registry.
  ///
  /// Register this LinearSolverFactory with the central registry, for
  /// the given SC, LO, GO, NT template parameters.  This will let any
  /// clients of Trilinos::Details::getLinearSolver create Ifpack2
  /// solvers with those template parameters.
  ///
  /// You may call this function multiple times; it will only have an
  /// effect the first time (it is idempotent).
  ///
  /// Users do not normally have to call this function.  Ifpack2
  /// automatically registers its LinearSolverFactory with the central
  /// repository, for all enabled template parameter combinations.
  static void registerLinearSolverFactory ();
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_LINEARSOLVERFACTORY_DECL_HPP
