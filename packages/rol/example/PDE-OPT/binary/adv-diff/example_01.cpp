// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the binary advection-diffusion control problem.
*/

#include "Teuchos_GlobalMPISession.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
//#include <fenv.h>

#include "ROL_Stream.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_PEBBL_BranchAndBound.hpp"
#include "ROL_PEBBL_TeuchosBranchHelper.hpp"
#include "branchHelper.hpp"

#include "opfactory.hpp"

int main(int argc, char *argv[]) {
  //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int myRank = comm->getRank();
  ROL::Ptr<std::ostream> outStream = ROL::makeStreamPtr( std::cout, (argc > 1) && (myRank==0) );

  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input_ex01.xml";
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile(filename);
    Teuchos::RCP<Teuchos::ParameterList> tparlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, tparlist.ptr() );

    /*************************************************************************/
    /***************** BUILD OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    bool useQP = parlist->sublist("Problem").get("Use Quadratic Program",false);
    ROL::Ptr<ROL::PEBBL::IntegerProblemFactory<RealT>> factory;
    ROL::Ptr<ROL::PEBBL::BranchHelper<RealT>> bHelper;
    if (!useQP) {
      RealT intTol = parlist->sublist("Problem").get("Integrality Tolerance",1e-6);
      int method = parlist->sublist("Problem").get("Branching Method",0);
      factory = ROL::makePtr<BinaryAdvDiffFactory<RealT>>(*parlist,comm,outStream);
      bHelper = ROL::makePtr<AdvDiffBranchHelper<RealT>>(intTol,method);
    }
    else {
      factory = ROL::makePtr<BinaryAdvDiffQPFactory<RealT>>(*parlist,comm,outStream);
      bHelper = ROL::makePtr<ROL::PEBBL::TeuchosBranchHelper<int,RealT>>();
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::PEBBL::BranchAndBound<RealT> pebbl(factory,parlist,bHelper,3,outStream);
    pebbl.solve(argc,argv,*outStream);

    //*outStream << "OPTIMAL CONTROLS" << std::endl;
    //for (int i=0; i<controlDim; ++i) {
    //  *outStream << (*z_ptr)[i] << "  ";
    //}
    //*outStream << std::endl;
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
