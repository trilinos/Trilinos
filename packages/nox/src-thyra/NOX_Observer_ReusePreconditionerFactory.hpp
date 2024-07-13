// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_OBSERVER_REUSE_PRECONDITIONER_FACTORY_HPP
#define NOX_OBSERVER_REUSE_PRECONDITIONER_FACTORY_HPP

#include "NOX_Config.h"
#include "NOX_Observer_ReusePreconditioner.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace NOX {

  /** Non-member function to create an observer that will reuse the preconditioner across multiple linear solves and time steps.

      When to update the preconditioner is controlled by the ParameterList arguments below.

      @param "Update prec at start of nonlinear solve" (bool) Update preconditioner at the start of each nonlinear solve. Defaults to true.
      @param "Update prec after this many nonlinear iterations" (int) Update preconditioner after this many nonlinear iterations. Setting to 0 disables. Defaults to 0 (disabled).
      @param "Update prec after this many stalled linear solves" (int) Update prec after this many stalled linear solves. Setting to 0 disables. Defaults to 0 (disabled). 
      @param "Max linear iterations for stall" (int) Maximum number of linear iterations that triggers a nonlinear iteration to be declared stalled. Defaults to 50.
   */
  Teuchos::RCP<NOX::Observer> createReusePreconditionerObserver(Teuchos::ParameterList& pl);

}

#endif
