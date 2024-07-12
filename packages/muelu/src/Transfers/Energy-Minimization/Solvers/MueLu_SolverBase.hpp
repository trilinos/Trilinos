// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SOLVERBASE_HPP
#define MUELU_SOLVERBASE_HPP

#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_Constraint_fwd.hpp"

namespace MueLu {

/*!
  @class SolverBase class.
  @brief Base class for energy-minimization iterative solvers
  */

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class SolverBase : public BaseClass {
#undef MUELU_SOLVERBASE_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //@{ Constructors/Destructors.
  SolverBase() {}

  virtual ~SolverBase() {}
  //@}

  //! @name Iterate methods.
  //@{

  /*! Iterate
   \param A -- Energy matrix
   \param C -- Constraints
   \param P0 -- Initial guess
   \param P -- Resulting prolongator
   \note We assume that P0 has correct sparsity pattern
   */
  virtual void Iterate(const Matrix& A, const Constraint& C, const Matrix& P0, RCP<Matrix>& P) const = 0;

  //@}

};  // class SolverBase

}  // namespace MueLu

#define MUELU_SOLVERBASE_SHORT

#endif  // ifndef MUELU_SOLVERBASE_HPP
