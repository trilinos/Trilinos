// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_01.hpp"

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
    /************************* CONSTRUCT ROL ALGORITHM ********************************************/
    /**********************************************************************************************/
    // Get ROL parameterlist
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    /**********************************************************************************************/
    /************************* CONSTRUCT VECTORS **************************************************/
    /**********************************************************************************************/
    // Build control vectors
    int N = 10;
    ROL::Ptr<std::vector<RealT>> x_ptr    = ROL::makePtr<std::vector<RealT>>(N,0.0);
    ROL::Ptr<ROL::Vector<RealT>> x        = ROL::makePtr<ROL::StdVector<RealT>>(x_ptr);
    std::vector<RealT> solution(N,0.0);
    /**********************************************************************************************/
    /************************* CONSTRUCT OBJECTIVE FUNCTION ***************************************/
    /**********************************************************************************************/
    std::vector<RealT> alpha(N,0.0);
    *outStream << std::endl;
    *outStream << "alpha =";
    for (int i = 0; i < N; ++i) {
      alpha[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      *outStream << "  " << alpha[i];
    }
    *outStream << std::endl << std::endl;
    ROL::Ptr<ROL::Objective<RealT>> obj
      = ROL::makePtr<Objective_SimpleBinary<RealT>>(alpha);
    /**********************************************************************************************/
    /************************* ENUMERATE SOLUTION *************************************************/
    /**********************************************************************************************/
    int budget = 3;
    for (int i = 0; i < budget; ++i) {
      (*x_ptr)[N-1-i] = 1.0;
    }
    // There are only 120 possible solutions
    RealT minval(ROL::ROL_INF<RealT>()), val(0), tol(1e-8);
    do {
      val = obj->value(*x,tol);
      if (val < minval) {
        minval = val;
        solution.assign(x_ptr->begin(),x_ptr->end());
      }
      //for (int i = 0; i < N; ++i) {
      //  *outStream << "  " << (*x_ptr)[i];
      //}
      //*outStream << std::endl << "  val = " << val << std::endl;
    } while (std::next_permutation(x_ptr->begin(),x_ptr->end()));

    *outStream << "Minimum Value:  " << minval << std::endl;
    *outStream << "Optimal Solution:";
    for (int i = 0; i < N; ++i) {
      *outStream << "  " << solution[i];
    }
    *outStream << std::endl << std::endl;

    RealT sum(0);
    std::sort(alpha.begin(),alpha.end());
    for (int i = 0; i < budget; ++i) {
      sum += alpha[i];
    }
    sum *= static_cast<RealT>(0.5);
    errorFlag += (std::abs(minval-sum)<ROL::ROL_EPSILON<RealT>() ? 0 : 1);
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
