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

/// \file Amesos2_Details_registerLinearSolverFactory.hpp
/// \author Mark Hoemmen
/// \brief Declaration of Amesos2::Details::registerLinearSolverFactory.

#ifndef AMESOS2_DETAILS_REGISTERLINEARSOLVERFACTORY_HPP
#define AMESOS2_DETAILS_REGISTERLINEARSOLVERFACTORY_HPP

#include "Amesos2_config.h"

namespace Amesos2 {
namespace Details {

/// \brief Register Amesos2's LinearSolverFactory with the central
///   repository, for all enabled combinations of template parameters.
///
/// For all combinations of template parameters that Amesos2 enables,
/// register Amesos2::Details::LinearSolverFactory with the central
/// repository.  This will let any clients of
/// Trilinos::Details::getLinearSolver create Amesos2 solvers with
/// those template parameters.
///
/// You may call this function multiple times; it will only have an
/// effect the first time (it is idempotent).
///
/// Users do not normally have to call this function, but see Bug
/// 6392.  Amesos2 tries its best to register its LinearSolverFactory
/// automatically with the central repository, for all enabled
/// template parameter combinations.  You may have to call this
/// function if your C++ compiler does not support the necessary
/// features for automatic registration to work, or if Trilinos was
/// configured with automatic registration turned off
/// (<tt>Trilinos_ENABLE_LINEAR_SOLVER_FACTORY_REGISTRATION</tt> was
/// set to <tt>OFF</tt>).  It never hurts to invoke this function
/// manually, though.
///
/// If you need to register Amesos2's LinearSolverFactory for a set of
/// template parameters that is <i>not</i> enabled, see
/// Amesos2_Details_LinearSolverFactory.hpp (in this directory).
void registerLinearSolverFactory ();

} // namespace Details
} // namespace Amesos2

namespace { // (anonymous)

// \class RegisterLinearSolverFactory
// \brief Register Amesos2's solver factory/ies with the central registry.
//
// \warning NOT FOR USERS.  ONLY FOR USE IN THIS FILE.
//
// Invoke this class' constructor to register Amesos2's solver
// factory/ies with the central registry, for all template parameter
// combinations that Amesos2 enabled.  You need not keep the instance
// of the class around; the constructor has a side effect if it
// returns.  (This is the C++ way of doing
// <tt>__attribute__((constructor))</tt>, without actually requiring
// the syntax extension.)
class RegisterLinearSolverFactory {
public:
  RegisterLinearSolverFactory () {
    Amesos2::Details::registerLinearSolverFactory ();
  }
};

// Creating an instance of RegisterLinearSolverFactory invokes its
// constructor, which has the side effect of calling
// Amesos2::Details::registerLinearSolverFactory().
RegisterLinearSolverFactory registerIt;

} // namespace (anonymous)

#endif // AMESOS2_DETAILS_REGISTERLINEARSOLVERFACTORY_HPP
