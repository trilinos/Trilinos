//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
//@HEADER

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
  runPreIterate(const NOX::Solver::Generic& solver)
  {
    pre_it_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPreIterate");
  }

  void ObserverLog::runPostIterate(const NOX::Solver::Generic& solver)
  {
    post_it_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPostIterate");
  }

  void ObserverLog::runPreSolve(const NOX::Solver::Generic& solver)
  {
    pre_solve_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPreSolve");
  }

  void ObserverLog::runPostSolve(const NOX::Solver::Generic& solver)
  {
    post_solve_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPostSolve");
  }

  void ObserverLog::runPreSolutionUpdate(const NOX::Abstract::Vector& update,
                                         const NOX::Solver::Generic& solver)
  {
    pre_solution_update_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPreSolutionUpdate");
  }

  void ObserverLog::runPostSolutionUpdate(const NOX::Solver::Generic& solver)
  {
    post_solution_update_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPostSolutionUpdate");
  }

  void ObserverLog::runPreLineSearch(const NOX::Solver::Generic& solver)
  {
    pre_linesearch_count_ += 1;
    if (log_call_order_)
      call_order_.push_back("runPreLineSearch");
  }

  void ObserverLog::runPostLineSearch(const NOX::Solver::Generic& solver)
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
