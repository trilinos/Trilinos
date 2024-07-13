// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_OBSERVER_PRINT_HPP
#define NOX_OBSERVER_PRINT_HPP

#include "NOX_Common.H"
#include "Teuchos_RCP.hpp"
#include "NOX_Observer.hpp"

namespace NOX {

  class Utils;

  /** \brief A NOX::Observer that provides summary solver output.

      This object demonstrates how to tailor the output from NOX for a
      particular user application. Users can disable the default
      output from NOX and use the NOX::Observer to output information
      in their own format.
  */
  class ObserverPrint : public NOX::Observer {
  public:
    ObserverPrint(const Teuchos::RCP<NOX::Utils>& printUtils);
    // Derived methods
    void runPreIterate(const NOX::Solver::Generic& solver);
    void runPostIterate(const NOX::Solver::Generic& solver);
  private:
    void printStep(const NOX::Solver::Generic& solver);
    Teuchos::RCP<NOX::Utils> os_;
  };
}

#endif
