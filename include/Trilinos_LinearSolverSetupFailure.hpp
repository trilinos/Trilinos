// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TRILINOS_LINEAR_SOLVER_SETUP_FAILURE_HPP
#define TRILINOS_LINEAR_SOLVER_SETUP_FAILURE_HPP

/// \file Trilinos_LinearSolverSetupFailure.hpp
/// \brief Declaration and definition of
///   Trilinos::LinearSolverSetupFailure exception class.

#include "TeuchosRemainder_config.h"
#include <stdexcept>
#include <string>

/// \namespace Trilinos
/// \brief Namespace of things generally useful to many Trilinos packages
namespace Trilinos {

/// \brief Exception reporting failure in setting up a linear solver.
///
/// This class exists so that NOX and Stratimikos can share a common
/// exception type, without introducing a dependency between the
/// packages.  See Trilinos GitHub Issue #1608.
///
/// "Failure to set up a linear solver" includes preconditioners and
/// sparse direct solvers.  It refers to the setup phase, not the
/// actual solve, for solvers that separate setup and solve into two
/// separate interfaces or function calls.  Setup does <i>not</i>
/// refer to a "failure to converge when performing an iterative
/// linear solve."
///
/// Here are some examples of legitimate "failures to set up a linear
/// solver":
///
/// <ul>
/// <li> Sparse LU factorization (either complete or incomplete) found
///      that the matrix is singular. </li>
/// <li> When setting up a block (a.k.a. physics-based, a.k.a. Teko)
///      preconditioner, setup found that one of the blocks is
///      singular. </li>
/// </ul>
///
/// Any code that throws this exception is responsible for ensuring
/// consistent behavior across multiple MPI processes.  Code that
/// catches this exception may assume that if any process in the
/// relevant communicator throws this exception, then all processes in
/// the relevant communicator will throw this exception, in such a way
/// as to ensure correct, deadlock-free MPI behavior.  For example,
/// suppose that a preconditioner implements additive Schwarz domain
/// decomposition, with one subdomain per MPI process.  Suppose that
/// the preconditioner uses a sparse direct solver on each subdomain.
/// If the sparse direct solver discovers that one subdomain has a
/// singular matrix, then the preconditioner is responsible for
/// ensuring that all processes in the matrix's communicator make a
/// consistent choice to throw the exception.
///
/// One would expect to catch and recover from this exception in a
/// nonlinear solver or time integrator.  A nonlinear solver might
/// recover by reducing the nonlinear step size.  A time integrator
/// might recover by reducing the time step size.  In general,
/// recovery from this exception would likely involve slowing down and
/// making choices in favor of stability.
///
/// In general, it's suboptimal coding practice to use exceptions for
/// control flow with nonexceptional conditions.  Better practice is
/// to return with a useful error code.  If a solver can detect
/// conditions like "the matrix is singular," I consider that a
/// nonexceptional condition.  However, it looks like NOX and
/// Stratimikos don't have a better way to communicate, so we'll do
/// our best for now by including this exception.
class LinearSolverSetupFailure : public std::runtime_error {
public:
  LinearSolverSetupFailure (const std::string& msg) :
    std::runtime_error (msg) {}
  virtual ~LinearSolverSetupFailure () = default;
};

} // namespace Trilinos

#endif // TRILINOS_LINEAR_SOLVER_SETUP_FAILURE_HPP
