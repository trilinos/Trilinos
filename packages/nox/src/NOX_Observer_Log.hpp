// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_OBSERVER_LOG_HPP
#define NOX_OBSERVER_LOG_HPP

#include "NOX_Common.H"
#include "NOX_Observer.hpp"
#include <vector>
#include <string>

namespace NOX {
  
  //! Logs observer calls. Useful for unit testing and debugging.
  class ObserverLog : public NOX::Observer {

  private:
    int pre_it_count_;
    int post_it_count_;
    int pre_solve_count_;
    int post_solve_count_;
    int pre_solution_update_count_;
    int post_solution_update_count_;
    int pre_linesearch_count_;
    int post_linesearch_count_;
    bool log_call_order_;
    std::vector<std::string> call_order_;

  public:
    ObserverLog(const bool log_call_order = true);
    ~ObserverLog();
    void runPreIterate(const NOX::Solver::Generic& solver);
    void runPostIterate(const NOX::Solver::Generic& solver);
    void runPreSolve(const NOX::Solver::Generic& solver);
    void runPostSolve(const NOX::Solver::Generic& solver);
    void runPreSolutionUpdate(const NOX::Abstract::Vector& update,
                              const NOX::Solver::Generic& solver);
    void runPostSolutionUpdate(const NOX::Solver::Generic& solver);
    void runPreLineSearch(const NOX::Solver::Generic& solver);
    void runPostLineSearch(const NOX::Solver::Generic& solver);
    int preIterateCount() const;
    int postIterateCount() const;
    int preSolveCount() const;
    int postSolveCount() const;
    int preSolutionUpdateCount() const;
    int postSolutionUpdateCount() const;
    int preLineSearchCount() const;
    int postLineSearchCount() const;
    const std::vector<std::string>& getCallOrder() const;
  };

}

#endif
