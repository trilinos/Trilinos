// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_OBSERVER_VECTOR_HPP
#define NOX_OBSERVER_VECTOR_HPP

#include "NOX_Common.H"
#include "Teuchos_RCP.hpp"
#include "NOX_Observer.hpp"
#include <vector>

namespace NOX {

  /** \brief Concrete implementation of NOX::Observer that stores a vector of Observers.

      The intent of this object to to aggregate a set of Observer objects.
  */
  class ObserverVector : public NOX::Observer {
  public:
    // Derived methods
    void runPreIterate(const NOX::Solver::Generic& solver);
    void runPostIterate(const NOX::Solver::Generic& solver);
    void runPreSolve(const NOX::Solver::Generic& solver);
    void runPostSolve(const NOX::Solver::Generic& solver);
    void runPreSolutionUpdate(const NOX::Abstract::Vector& update, const NOX::Solver::Generic& solver);
    void runPostSolutionUpdate(const NOX::Solver::Generic& solver);
    void runPreLineSearch(const NOX::Solver::Generic& solver);
    void runPostLineSearch(const NOX::Solver::Generic& solver);

    //! Add observer to end of vector.
    void pushBack(const Teuchos::RCP<NOX::Observer>& observer);
    //! Remove observer from end of vector.
    void popBack();
    //! Clear the vector of observers.
    void clear();

  private:
    //! std::vector of observer objects
    std::vector<Teuchos::RCP<NOX::Observer>> vec_;
  };
}

#endif
