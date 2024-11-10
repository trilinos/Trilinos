// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEST_SINGLE_AZTECOO_THYRA_SOLVER_HPP
#define TEST_SINGLE_AZTECOO_THYRA_SOLVER_HPP

#include "AztecOO_config.h"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"

namespace Teuchos { class ParameterList; }

namespace Thyra {

/** \brief Testing function for a single aztecoo solver with a single matrix.
 *
 */
bool test_single_aztecoo_thyra_solver(
  const std::string                       matrixFile
  ,const bool                             testTranspose
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxResid
  ,const double                           maxSolutionError
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::ParameterList                 *paramList
  ,Teuchos::FancyOStream                  *out
  );

} // namespace Thyra

#endif // TEST_SINGLE_AZTECOO_THYRA_SOLVER_HPP
