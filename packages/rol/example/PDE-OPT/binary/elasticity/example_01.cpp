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

/*! \file  example_01.cpp
    \brief Shows how to solve the stuctural topology optimization problem.

*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_PEBBL_BranchAndBound.hpp"

#include "opfactory.hpp"
#include "transform.hpp"
#include "branching.hpp"
#include "branchHelper.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int>> comm
    = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();
  if ((iprint> 0) && (myRank == 0)) {
    outStream = ROL::makePtrFromRef(std::cout);
  }
  else {
    outStream = ROL::makePtrFromRef(bhs);
  }
  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input_ex01.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const std::string example = parlist->sublist("Problem").get("Example", "Default");
    int probDim               = parlist->sublist("Problem").get("Problem Dimension", 2);
    RealT maxWeight           = parlist->sublist("Problem").get("Maximum Weight Fraction", 1.0);
    if (example == "2D Wheel"                   ||
        example == "2D Truss"                   ||
        example == "2D Cantilever with 1 Load"  ||
        example == "2D Cantilever with 3 Loads" ||
        example == "2D Beams"                   ||
        example == "2D Carrier Plate") {
      probDim = 2;
    }
    else if (example == "3D Cantilever") {
      probDim = 3;
    }
    *outStream << std::endl;
    *outStream << "Problem Data"            << std::endl;
    *outStream << "  Example:             " << example   << std::endl;
    *outStream << "  Dimension:           " << probDim   << std::endl;
    *outStream << "  Max Weight Fraction: " << maxWeight << std::endl;
    *outStream << std::endl;

    // Set up and solve.
    ROL::Ptr<ElasticityFactory<RealT>> factory
      = ROL::makePtr<ElasticityFactory<RealT>>(*parlist,comm,outStream);
    bool derivCheck = parlist->sublist("Problem").get("Check derivatives",false);
    if (derivCheck) {
      factory->check();
      ROL::Ptr<ROL::PEBBL::IntegerProblem<RealT>> prob = factory->build();
      ROL::Ptr<TpetraMultiMatIntegerTransformation<RealT>> trans
        = ROL::makePtr<TpetraMultiMatIntegerTransformation<RealT>>();
      trans->add({  1, RealT(1)});
      trans->add({ 17, RealT(0)});
      trans->add({107, RealT(1)});
      prob->setTransformation(trans);
      prob->finalize(false,true,*outStream);
      prob->check(true,*outStream);
    }
    Teuchos::Time algoTimer("Algorithm Time", true);
    bool binary = parlist->sublist("Problem").get("Binary",true);
    ROL::Ptr<ROL::Vector<RealT>> z;
    if (!binary) {
      ROL::Ptr<ROL::PEBBL::IntegerProblem<RealT>> problem = factory->build();
      problem->finalize(false,true,*outStream);
      ROL::Solver<RealT> solver(problem,*parlist);
      solver.solve(*outStream);
      z = problem->getPrimalOptimizationVector();
    }
    else {
      RealT intTol  = parlist->sublist("Problem").get("Integrality Tolerance",1e-6);
      int method    = parlist->sublist("Problem").get("Branching Method",0);
      int verbosity = parlist->sublist("Problem").get("BB Output Level",0);
      ROL::Ptr<TpetraMultiMatBranchHelper<RealT>> bHelper
        = ROL::makePtr<TpetraMultiMatBranchHelper<RealT>>(intTol,method);
      ROL::Ptr<MultiMatBranching<RealT>> branching
        = ROL::makePtr<MultiMatBranching<RealT>>(factory,parlist,bHelper,verbosity,outStream);
      ROL::PEBBL::BranchAndBound<RealT> pebbl(branching);
      pebbl.solve(argc,argv,*outStream);
      z = pebbl.getSolution()->clone();
      z->set(*pebbl.getSolution());
    }
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    factory->print(*z);
    
    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
