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
    \brief Test StatusTest input mechanism in OptimizationSolver.
*/

#include "ROL_TestObjectives.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>

typedef double RealT;

template<class Real>
class myStatusTest : public ROL::StatusTest<Real> {
private:
  const Real ftol_;
  Real f0_;

public:
  myStatusTest(const Real ftol) : ftol_(ftol), f0_(ROL::ROL_INF<Real>()) {}

  bool check( ROL::AlgorithmState<Real> &state ) {
    const Real zero(0);
    Real snorm = state.snorm;
    Real fval  = state.value;
    Real diff  = std::abs(f0_-fval);
    f0_ = fval;
    return ((diff < ftol_ && snorm > zero) ? false : true);
  }
};

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Test body.

  try {
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Setup optimization problem
    Teuchos::RCP<ROL::Vector<RealT> > x0, z;
    Teuchos::RCP<ROL::Objective<RealT> > obj;
    Teuchos::RCP<ROL::BoundConstraint<RealT> > bnd;
    ROL::getTestObjectives<RealT>(obj,bnd,x0,z,ROL::TESTOPTPROBLEM_HS1);
    Teuchos::RCP<ROL::Vector<RealT> > x = x0->clone(); x->set(*x0);
    ROL::OptimizationProblem<RealT> optProblem(obj,x,bnd);

    // Get Dimension of Problem
    int dim = x0->dimension(); 
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit", 2*dim);

    // Check Derivatives
    optProblem.check(*outStream);

    // Error Vector
    Teuchos::RCP<ROL::Vector<RealT> > e = x0->clone();
    e->zero();

    // Setup optimization solver
    parlist->sublist("Status Test").set("Gradient Tolerance",static_cast<RealT>(1e-6));
    parlist->sublist("Step").set("Type", "Interior Point");
    ROL::OptimizationSolver<RealT> optSolver(optProblem,*parlist);
    Teuchos::RCP<ROL::StatusTest<RealT> > myStatus
      = Teuchos::rcp( new myStatusTest<RealT>(1.e-8) );
    optSolver.solve(*outStream, myStatus, false);

    // Compute Error
    e->set(*x);
    e->axpy(static_cast<RealT>(-1),*z);
    RealT enorm = e->norm();
    *outStream << std::endl << "Norm of Error: " << enorm << std::endl;

    RealT tol = static_cast<RealT>(1e-3)*std::max(z->norm(),static_cast<RealT>(1));
    errorFlag += ((enorm < tol) ? 0 : 1);
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

