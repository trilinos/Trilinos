// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Observer_Log.hpp"

namespace NOX {

  ObserverLog::ObserverLog(const bool log_call_order) :
    pre_it_count_(0),
    post_it_count_(0),
    pre_solve_count_(0),
    post_solve_count_(0),
    pre_solution_update_count_(0),
    post_solution_update_count_(0),
    pre_linesearch_count_(0),
    post_linesearch_count_(0),
    log_call_order_(log_call_order)
  {}

  ObserverLog::~ObserverLog() {}

  void ObserverLog::
  runPreIterate(const NOX::Solver::Generic& /* solver */)
  {
    pre_it_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPreIterate");
  }

  void ObserverLog::runPostIterate(const NOX::Solver::Generic& /* solver */)
  {
    post_it_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPostIterate");
  }

  void ObserverLog::runPreSolve(const NOX::Solver::Generic& /* solver */)
  {
    pre_solve_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPreSolve");
  }

  void ObserverLog::runPostSolve(const NOX::Solver::Generic& /* solver */)
  {
    post_solve_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPostSolve");
  }

  void ObserverLog::runPreSolutionUpdate(const NOX::Abstract::Vector& /* update */,
                                         const NOX::Solver::Generic& /* solver */)
  {
    pre_solution_update_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPreSolutionUpdate");
  }

  void ObserverLog::runPostSolutionUpdate(const NOX::Solver::Generic& /* solver */)
  {
    post_solution_update_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPostSolutionUpdate");
  }

  void ObserverLog::runPreLineSearch(const NOX::Solver::Generic& /* solver */)
  {
    pre_linesearch_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPreLineSearch");
  }

  void ObserverLog::runPostLineSearch(const NOX::Solver::Generic& /* solver */)
  {
    post_linesearch_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPostLineSearch");
  }

  int ObserverLog::preIterateCount() const
  {return pre_it_count_;}

  int ObserverLog::postIterateCount() const
  {return post_it_count_;}

  int ObserverLog::preSolveCount() const
  {return pre_solve_count_;}

  int ObserverLog::postSolveCount() const
  {return post_solve_count_;}

  int ObserverLog::preSolutionUpdateCount() const
  {return pre_solution_update_count_;}

  int ObserverLog::postSolutionUpdateCount() const
  {return post_solution_update_count_;}

  int ObserverLog::preLineSearchCount() const
  {return pre_linesearch_count_;}

  int ObserverLog::postLineSearchCount() const
  {return post_linesearch_count_;}

  const std::vector<std::string>& ObserverLog::getCallOrder() const
  {return call_order_;}

}
