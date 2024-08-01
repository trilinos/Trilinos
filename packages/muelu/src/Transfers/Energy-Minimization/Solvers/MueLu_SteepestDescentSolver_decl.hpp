// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_STEEPESTDESCENTSOLVER_DECL_HPP
#define MUELU_STEEPESTDESCENTSOLVER_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_CrsMatrixFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_SolverBase.hpp"
#include "MueLu_Constraint_fwd.hpp"

namespace MueLu {

/*!
  @class SteepestDescentSolver class.
  @brief Implements steepest descent algorithm for energy-minimization
  */

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class SteepestDescentSolver : public SolverBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_STEEPESTDESCENTSOLVER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  /*! Constructor
    \param Its -- Number of performed iterations
    \param StepLength -- The modifier of the step length
    */
  SteepestDescentSolver(size_t Its, SC StepLength = Teuchos::ScalarTraits<Scalar>::one());

  //@}

  //! @name Iterate methods.
  //@{

  //! Iterate
  void Iterate(const Matrix& A, const Constraint& C, const Matrix& P0, RCP<Matrix>& P) const;

  //@}

 private:
  size_t nIts_;    //!< Number of performed iterations
  SC stepLength_;  //!< Modifier of the step length
};

}  // namespace MueLu

#define MUELU_STEEPESTDESCENTSOLVER_SHORT
#endif  // MUELU_STEEPESTDESCENTSOLVER_DECL_HPP
