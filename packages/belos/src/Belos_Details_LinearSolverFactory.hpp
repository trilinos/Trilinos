//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef BELOS_DETAILS_LINEARSOLVERFACTORY_HPP
#define BELOS_DETAILS_LINEARSOLVERFACTORY_HPP

/// \file Belos_Details_LinearSolverFactory.hpp
/// \brief Implementation of Trilinos::Details::LinearSolverFactory.

#include "BelosSolverFactory.hpp"
#include "Belos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"

namespace Belos {
namespace Details {

/// \class LinearSolver
/// \brief Belos' implementation of Trilinos::Details::LinearSolverFactory.
template<class MV, class OP, class ScalarType, class NormType>
class LinearSolverFactory :
    public Trilinos::Details::LinearSolverFactory<MV, OP, NormType>
{
public:
  /// \brief Get an instance of a Belos solver.
  ///
  /// The solver is wrapped in a Trilinos::Details::LinearSolver
  /// interface.
  ///
  /// \param solverName [in] The solver's name.  Names are case
  ///   sensitive.
  /// \return A pointer to the solver, if the name was valid; else,
  ///   a null pointer (Teuchos::null).
  virtual Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, NormType> >
  getLinearSolver (const std::string& solverName)
  {
    using Teuchos::rcp;
    return rcp (new Belos::Details::LinearSolver<MV, OP, ScalarType, NormType> (solverName));
  }

  /// \brief Register this LinearSolverFactory with the central registry.
  ///
  /// Register this LinearSolverFactory with the central registry, for
  /// the given SC, LO, GO, NT template parameters.  This will let any
  /// clients of Trilinos::Details::getLinearSolver create Belos
  /// solvers with those template parameters.
  ///
  /// You may call this function multiple times; it will only have an
  /// effect the first time (it is idempotent).
  ///
  /// Users do not normally have to call this function.  Belos
  /// automatically registers its LinearSolverFactory with the central
  /// repository, for all enabled template parameter combinations.
  static void registerLinearSolverFactory ()
  {
    typedef Belos::Details::LinearSolverFactory<MV, OP, ScalarType, NormType> this_type;

#ifdef HAVE_TEUCHOSCORE_CXX11
    typedef std::shared_ptr<this_type> ptr_type;
    //typedef std::shared_ptr<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#else
    typedef Teuchos::RCP<this_type> ptr_type;
    //typedef Teuchos::RCP<Trilinos::Details::LinearSolverFactory<MV, OP> > base_ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

    ptr_type factory (new this_type ());
    Trilinos::Details::registerLinearSolverFactory<MV, OP, NormType> ("Belos", factory);
  }
};

} // namespace Details
} // namespace Belos

#endif /* BELOS_DETAILS_LINEARSOLVERFACTORY_HPP */
