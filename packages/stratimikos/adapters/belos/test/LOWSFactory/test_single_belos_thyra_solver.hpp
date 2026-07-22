// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEST_SINGLE_BELOS_THYRA_SOLVER_HPP
#define TEST_SINGLE_BELOS_THYRA_SOLVER_HPP

#include "BelosConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"

namespace Teuchos { class ParameterList; }

namespace Thyra {

/** \brief Testing function for a single belos solver with a single matrix.
 *
 */
bool test_single_belos_thyra_solver(
  const std::string                       matrixFile
  ,const bool                             testTranspose
  ,const bool                             usePreconditioner
  ,const int                              numRhs
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxResid
  ,const double                           maxSolutionError
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::ParameterList                 *belosLOWSFPL
  ,Teuchos::ParameterList                 *precPL
  ,Teuchos::FancyOStream                  *out
  );

} // namespace Thyra

#endif // TEST_SINGLE_BELOS_THYRA_SOLVER_HPP
