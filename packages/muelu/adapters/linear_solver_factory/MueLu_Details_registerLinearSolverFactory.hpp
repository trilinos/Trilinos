// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file MueLu_Details_registerLinearSolverFactory.hpp
/// \authors Alicia Klinvex and Mark Hoemmen
/// \brief Declaration of MueLu::Details::registerLinearSolverFactory.

#ifndef MUELU_DETAILS_REGISTERLINEARSOLVERFACTORY_HPP
#define MUELU_DETAILS_REGISTERLINEARSOLVERFACTORY_HPP

#include "MueLu_config.hpp"

namespace MueLu {
namespace Details {

/// \brief Register MueLu's LinearSolverFactory with the central
///   repository, for all enabled combinations of template parameters.
///
/// For all combinations of template parameters that MueLu enables,
/// register MueLu::Details::LinearSolverFactory with the central
/// repository.  This will let any clients of
/// Trilinos::Details::getLinearSolver create MueLu solvers with
/// those template parameters.
///
/// You may call this function multiple times; it will only have an
/// effect the first time (it is idempotent).
///
/// Users do not normally have to call this function, but see Bug
/// 6392.  MueLu tries its best to register its LinearSolverFactory
/// automatically with the central repository, for all enabled
/// template parameter combinations.  You may have to call this
/// function if your C++ compiler does not support the necessary
/// features for automatic registration to work, or if Trilinos was
/// configured with automatic registration turned off
/// (<tt>Trilinos ENABLE_LINEAR_SOLVER_FACTORY_REGISTRATION</tt> was
/// set to <tt>OFF</tt>).  It never hurts to invoke this function
/// manually, though.
///
/// If you need to register MueLu's LinearSolverFactory for a set of
/// template parameters that is <i>not</i> enabled, see
/// MueLu_Details_LinearSolverFactory.hpp (in this directory).
void registerLinearSolverFactory();

}  // namespace Details
}  // namespace MueLu

namespace {  // (anonymous)

// \class RegisterLinearSolverFactory
// \brief Register MueLu's solver factory/ies with the central registry.
//
// \warning NOT FOR USERS.  ONLY FOR USE IN THIS FILE.
//
// Invoke this class' constructor to register MueLu's solver
// factory/ies with the central registry, for all template paramete
// combinations that MueLu enabled.  You need not keep the instance
// of the class around; the constructor has a side effect if it
// returns.  (This is the C++ way of doing
// <tt>__attribute__((constructor))</tt>, without actually requiring
// the syntax extension.)
class RegisterLinearSolverFactory {
 public:
  RegisterLinearSolverFactory() {
    MueLu::Details::registerLinearSolverFactory();
  }
};

// Creating an instance of RegisterLinearSolverFactory invokes its
// constructor, which has the side effect of calling
// MueLu::Details::registerLinearSolverFactory().
RegisterLinearSolverFactory registerIt;

}  // namespace

#endif  // MUELU_DETAILS_REGISTERLINEARSOLVERFACTORY_HPP