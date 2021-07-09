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

#define USE_HESSVEC 1

#include "ROL_GetTestProblems.hpp"
#include "ROL_TypeU_TrustRegionAlgorithm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
//#include <fenv.h>

typedef double RealT;

int main(int argc, char *argv[]) {
  //feenableexcept(FE_INVALID | FE_OVERFLOW);

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Test body.

  try {

    std::string filename = "input.xml";
    
    auto parlist = ROL::getParametersFromXmlFile( filename );
    parlist->sublist("General").set("Output Level", iprint);
    ROL::Ptr<ROL::TypeU::Algorithm<RealT>> algo;
    ROL::Ptr<ROL::Vector<RealT>> e, x, x0;

    for ( ROL::ETestOptProblem objFunc = ROL::TESTOPTPROBLEM_ROSENBROCK;
          objFunc < ROL::TESTOPTPROBLEM_LAST; objFunc++ ) {
      std::vector<ROL::Ptr<ROL::Vector<RealT>>> z;
      ROL::Ptr<ROL::OptimizationProblem<RealT>> problem;
      ROL::GetTestProblem<RealT>(problem,x0,z,objFunc);
      e = x0->clone();
      x = x0->clone();
      if (problem->getProblemType() == ROL::TYPE_U
          && objFunc != ROL::TESTOPTPROBLEM_MINIMAX1
          && objFunc != ROL::TESTOPTPROBLEM_MINIMAX2
          && objFunc != ROL::TESTOPTPROBLEM_MINIMAX3) {
        // Get Dimension of Problem
        int dim = x0->dimension(); 
        parlist->sublist("General").sublist("Krylov").set("Iteration Limit", 2*dim);
        parlist->sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Iteration Limit", 50);
        parlist->sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Maximum Spectral Step Size", 1e12);
        parlist->sublist("Step").sublist("Trust Region").sublist("SPG").sublist("Solver").set("Minimum Spectral Step Size", 1e-12);

        parlist->sublist("Step").set("Type","Trust Region");
        for ( ROL::ETrustRegionU desc = ROL::TRUSTREGION_U_CAUCHYPOINT;
              desc < ROL::TRUSTREGION_U_LAST; desc++ ) {

          parlist->sublist("Step").sublist("Trust Region").set("Subproblem Solver", ROL::ETrustRegionUToString(desc));
          *outStream << std::endl << std::endl
                     << ROL::ETestOptProblemToString(objFunc)
                     << std::endl << std::endl
                     << std::endl << std::endl
                     << ROL::ETrustRegionUToString(desc)
                     << std::endl << std::endl;

          // Define Line Search Solver
          algo = ROL::makePtr<ROL::TypeU::TrustRegionAlgorithm<RealT>>(*parlist);

          // Run Solver
          x->set(*x0);
          algo->run(*x,
                    *problem->getObjective(),
                    *outStream);

          // Compute Error
          e->zero();
          RealT err(0);
          for (int i = 0; i < static_cast<int>(z.size()); ++i) {
            e->set(*x);
            e->axpy(-1.0,*z[i]);
            if (i == 0) {
              err = e->norm();
            }
            else {
              err = std::min(err,e->norm());
            }
          }
          *outStream << std::endl << "Norm of Error: " << err << std::endl;
          //errorFlag += (int)(e.norm() < std::sqrt(ROL::ROL_EPSILON<RealT>())); 
        }
      }
    }
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;

}
