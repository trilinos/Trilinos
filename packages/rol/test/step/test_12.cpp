// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test line search.
*/

#define USE_HESSVEC 0

#include "ROL_Step.hpp"
#include "ROL_GetTestProblems.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;
  ROL::Ptr<ROL::Vector<RealT> > x0, x;
  ROL::Ptr<ROL::OptimizationProblem<RealT> > problem;
  ROL::GetTestProblem<RealT>(problem,x0,x,ROL::TESTOPTPROBLEM_HS1);
  ROL::Ptr<ROL::Objective<RealT> > obj = problem->getObjective();
  ROL::Ptr<ROL::BoundConstraint<RealT> > bnd = problem->getBoundConstraint();
  ROL::AlgorithmState<RealT> algo_state;

  // *** Test body.
  ROL::Step<RealT> step;
  int thrown = 0;
  try {
    try {
      step.compute(*x0,*x,*obj,*bnd,algo_state);
    }
    catch (ROL::Exception::NotImplemented exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    try {
      step.update(*x0,*x,*obj,*bnd,algo_state);
    }
    catch (ROL::Exception::NotImplemented exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    try {
      step.printName();
    }
    catch (ROL::Exception::NotImplemented exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    try {
      step.printHeader();
    }
    catch (ROL::Exception::NotImplemented exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    try {
      step.print(algo_state,true);
    }
    catch (ROL::Exception::NotImplemented exc) {
      *outStream << exc.what() << std::endl;
      thrown++;
    };
    errorFlag = (thrown==5) ? 0 : 1;
  }
  catch (std::logic_error err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;

}

