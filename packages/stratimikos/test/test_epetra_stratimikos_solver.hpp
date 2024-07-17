// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEST_SINGLE_EPETRA_STRATIMIKOS_SOLVER_HPP
#define TEST_SINGLE_EPETRA_STRATIMIKOS_SOLVER_HPP

#include "Stratimikos_Config.h"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"

namespace Teuchos { class ParameterList; }

namespace Thyra {

/** \brief Testing function for a single epetra stratimikos solver for a single
 * matrix.
 *
 * \ingroup stratimikos_testing_grp
 */
bool test_epetra_stratimikos_solver(
  Teuchos::ParameterList                  *paramList
  ,const bool                             dumpAll
  ,Teuchos::FancyOStream                  *out
  );

} // namespace Thyra

#endif // TEST_SINGLE_EPETRA_STRATIMIKOS_SOLVER_HPP
