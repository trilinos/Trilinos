// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_SOLVER_SOLVERUTILS_H
#define NOX_SOLVER_SOLVERUTILS_H

#include "Teuchos_RCP.hpp"
#include "NOX_StatusTest_Generic.H"

namespace Teuchos {
  class ParameterList;
}

namespace NOX {

  class Observer;

  namespace Solver{

    /*! \brief Nonmember method that returns the status test check type.

        This object parses the "Solver Options" parameter list for a
        key "Status Test Check Type" of type <std::string> with
        possible values:

    - "Complete"
    - "Minimal" (default)
    - "None"

    These options correspond to the NOX::StatusTest::CheckType.
    Please follow the link for this object for more information.
    */
    NOX::StatusTest::CheckType
    parseStatusTestCheckType(Teuchos::ParameterList& solver_options_list);

    /** \brief Parses the "Solver Options" sublist for a user defined observer. */
    Teuchos::RCP<NOX::Observer> parseObserver(Teuchos::ParameterList& solver_options_list);

    /*! \brief Nonmember method that validates and sets the defaults for the "Solver Options" sublist. */
    void
    validateSolverOptionsSublist(Teuchos::ParameterList& solver_options_list);

  }  // namespace Solver
}  // namespace NOX
#endif
