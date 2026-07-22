// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEST_SINGLE_AMESOS2_TPETRA_SOLVER_HPP
#define TEST_SINGLE_AMESOS2_TPETRA_SOLVER_HPP

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;

  template<class OrdinalType>
  class Comm;
}

namespace Thyra {

/** \brief Testing function for a single belos solver with a single matrix.
 *
 */
bool test_single_amesos2_tpetra_solver(
  const std::string                       matrixFile
  ,const int                              numRhs
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxResid
  ,const double                           maxSolutionError
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::ParameterList                 *amesos2LOWSFPL
  ,Teuchos::FancyOStream                  *out
  ,const Teuchos::RCP<const Teuchos::Comm<int> >& comm
  );

} // namespace Thyra

#endif // TEST_SINGLE_AMESOS2_TPETRA_SOLVER_HPP
