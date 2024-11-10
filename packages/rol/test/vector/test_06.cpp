// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_06.cpp
    \brief Test ProfiledVector interface.
*/

#include "ROL_ProfiledVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Zakharov.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Stream.hpp"

#include "Teuchos_GlobalMPISession.hpp"

typedef int    OrdinalT;
typedef double RealT;

template<>
ROL::VectorFunctionCalls<int>
ROL::ProfiledVector<int,RealT>::functionCalls_ = ROL::VectorFunctionCalls<int>();

int main(int argc, char *argv[]) {


  using ROL::ParameterList;

  typedef std::vector<RealT>                  vector;

  typedef ROL::Vector<RealT>                  V;
  typedef ROL::StdVector<RealT>               SV;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);

  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {
    // Dimension of the optimization vector
    int dim = 10;

    std::string paramfile = "parameters.xml";
    auto parlist = ROL::getParametersFromXmlFile(paramfile);

    // Define algorithm.
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);

    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(dim,1.0);
    ROL::Ptr<vector> k_ptr = ROL::makePtr<vector>(dim);

    for(int i=0;i<dim;++i) {
      (*k_ptr)[i] = 1.0 + i;
    }

    ROL::Ptr<V> xs = ROL::makePtr<SV>(x_ptr);
    ROL::Ptr<V> ks = ROL::makePtr<SV>(k_ptr);

    // Create ProfiledVector objects
    ROL::ProfiledVector<int,RealT> xpf(xs);
    ROL::Ptr<V> kpf = ROL::makePtr<ROL::ProfiledVector<int,RealT>>(ks);

    ROL::ZOO::Objective_Zakharov<RealT> obj(kpf);

    // Run algorithm.
    algo.run(xpf, obj, true, *outStream);

    // Display number of function calls
    ROL::printVectorFunctionCalls(xpf, *outStream);
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
