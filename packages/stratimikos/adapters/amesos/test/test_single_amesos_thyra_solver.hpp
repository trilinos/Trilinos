// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEST_SINGLE_AMESOS_THYRA_SOLVER_HPP
#define TEST_SINGLE_AMESOS_THYRA_SOLVER_HPP

#include "Thyra_AmesosTypes.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Thyra {

/** \brief Testing function for a single amesos solver with a single matrix.
 *
 */
bool test_single_amesos_thyra_solver(
  const std::string                       matrixFile
  ,Teuchos::ParameterList                 *amesosLOWSFPL
  ,const bool                             testTranspose
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxError
  ,const double                           maxResid
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::FancyOStream                  *out
  );

} // namespace Thyra

#endif // TEST_SINGLE_AMESOS_THYRA_SOLVER_HPP
