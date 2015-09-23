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

#ifndef BELOS_DETAILS_REGISTERLINEARSOLVERFACTORY_HPP
#define BELOS_DETAILS_REGISTERLINEARSOLVERFACTORY_HPP

/// \file Belos_Details_registerLinearSolverFactory.hpp
/// \brief Declaration of Belos::Details::registerLinearSolverFactory.

#include "BelosConfigDefs.hpp"

namespace Belos {
namespace Details {

/// \brief Register Belos' LinearSolverFactory with the central
///   repository, for all enabled combinations of template parameters.
///
/// For all combinations of template parameters that Belos enables,
/// register Belos::Details::LinearSolverFactory with the central
/// repository.  This will let any clients of
/// Trilinos::Details::getLinearSolver create Belos solvers with those
/// template parameters.
///
/// You may call this function multiple times; it will only have an
/// effect the first time (it is idempotent).
///
/// Users do not normally have to call this function, but see Bug
/// 6392.  Belos tries its best to register its LinearSolverFactory
/// automatically with the central repository, for all enabled
/// template parameter combinations.  You may have to call this
/// function if your C++ compiler does not support the necessary
/// features for automatic registration to work, or if Trilinos was
/// configured with automatic registration turned off
/// (<tt>Trilinos_ENABLE_LINEAR_SOLVER_FACTORY_REGISTRATION</tt> was
/// set to <tt>OFF</tt>).  It never hurts to invoke this function
/// manually, though.
///
/// If you need to register Belos's LinearSolverFactory for a set of
/// template parameters that is <i>not</i> enabled, see
/// Belos_Details_LinearSolverFactory.hpp (in this directory).
///
/// \warning FIXME (mfh 23 Aug 2015) This currently only works if the
///   compiler understands GCC's __attribute__((weak)) syntax.  See
///   the comments in Belos_Details_registerLinearSolverFactory.cpp in
///   this directory.  As a work-around, you may invoke
///   Belos::Details::Tpetra::registerLinearSolverFactory() for the
///   Tpetra specialization, or
///   Belos::Details::Epetra::registerLinearSolverFactory() for the
///   Epetra specialization.  Either of these requires an extern
///   declaration in your code.
void registerLinearSolverFactory ();

} // namespace Details
} // namespace Belos

namespace { // (anonymous)

// \class RegisterLinearSolverFactory
// \brief Register Belos' solver factory/ies with the central registry.
//
// \warning NOT FOR USERS.  ONLY FOR USE IN THIS FILE.
//
// Invoke this class' constructor to register Belos's solver
// factory/ies with the central registry, for all template parameter
// combinations that Belos enabled.  You need not keep the instance of
// the class around; the constructor has a side effect if it returns.
// (This is the C++ way of doing
// <tt>__attribute__((constructor))</tt>, without actually requiring
// the syntax extension.)
class RegisterLinearSolverFactory {
public:
  RegisterLinearSolverFactory () {
    Belos::Details::registerLinearSolverFactory ();
  }
};

// Creating an instance of RegisterLinearSolverFactory invokes its
// constructor, which has the side effect of calling
// Belos::Details::registerLinearSolverFactory().
RegisterLinearSolverFactory registerIt;

} // namespace (anonymous)

#endif /* BELOS_DETAILS_REGISTERLINEARSOLVERFACTORY_HPP */
