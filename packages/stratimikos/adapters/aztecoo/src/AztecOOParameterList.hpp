// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SET_AZTECOO_PARAMETER_LIST_HPP
#define THYRA_SET_AZTECOO_PARAMETER_LIST_HPP

#include "AztecOO.h"
#include "Teuchos_RCP.hpp"

/** \brief Setup an AztecOO solver object with a set of parameters.
 *
 * ToDo: Finish documentation!
 */
void setAztecOOParameters(
  Teuchos::ParameterList  *pl
  ,AztecOO                *solver
  );

/** \brief Return the list of all valid AztecOO parameters (to validate
 * against).
 *
 * ToDo: Finish documentation!
 */
Teuchos::RCP<const Teuchos::ParameterList> getValidAztecOOParameters();

#endif // THYRA_SET_AZTECOO_PARAMETER_LIST_HPP
