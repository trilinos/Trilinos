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
