// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
/// manually, though, since it is idempotent.
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
