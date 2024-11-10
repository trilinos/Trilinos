// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_DistributionFactory.hpp"

typedef double RealT;

int main(int argc, char* argv[]) {

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

  try {
    ROL::Ptr<ROL::Distribution<RealT> > dist;

    // Get ROL parameterlist
    std::string filename = "input_02.xml";

    auto parlist = ROL::getParametersFromXmlFile( filename ); 
 
    for (ROL::EDistribution ed = ROL::DISTRIBUTION_ARCSINE; ed != ROL::DISTRIBUTION_LAST; ed++) {
      *outStream << ROL::EDistributionToString(ed) << std::endl << std::endl;
      parlist->sublist("SOL").sublist("Distribution").set("Name",ROL::EDistributionToString(ed));
      dist = ROL::DistributionFactory<RealT>(*parlist);
      dist->test(*outStream);
    }
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
