// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_05.hpp"
#include "ROL_PEBBL_BranchAndBound.hpp"
#include "ROL_PEBBL_StdBranchHelper.hpp"

typedef double RealT;

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  try {
    /**********************************************************************************************/
    /************************* GET PROBLEM PARAMETERS *********************************************/
    /**********************************************************************************************/
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile("input.xml");
    /**********************************************************************************************/
    /************************* CONSTRUCT PROBLEM FACTORY ******************************************/
    /**********************************************************************************************/
    ROL::Ptr<Test05Factory<RealT>> factory = ROL::makePtr<Test05Factory<RealT>>(*parlist);
    /**********************************************************************************************/
    /************************* SOLVE **************************************************************/
    /**********************************************************************************************/
    ROL::Ptr<ROL::PEBBL::StdBranchHelper<RealT>> bHelper
      = ROL::makePtr<ROL::PEBBL::StdBranchHelper<RealT>>();
    ROL::PEBBL::BranchAndBound<RealT> pebbl(factory,parlist,bHelper,3,outStream);
    pebbl.solve(argc,argv,*outStream);
    *outStream << "Solution Vector:" << std::endl << "  ";
    pebbl.getSolution()->print(*outStream);
    *outStream << std::endl;
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
