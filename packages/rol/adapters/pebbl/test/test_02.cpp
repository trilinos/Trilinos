// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "ROL_PEBBL_IntegerConstraint.hpp"
#include "ROL_StdVector.hpp"

typedef double RealT;

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Teuchos::DefaultComm<int>::getComm();

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0 && Teuchos::rank<int>(*comm)==0)
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
    ROL::Ptr<std::vector<RealT>> d_ptr    = ROL::makePtr<std::vector<RealT>>(N,0.0);
    ROL::Ptr<ROL::Vector<RealT>> d        = ROL::makePtr<ROL::StdVector<RealT>>(d_ptr);
    for (int i = 0; i < N; ++i) {
      (*x_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      (*d_ptr)[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
    }
    /**********************************************************************************************/
    /************************* CONSTRUCT CONSTRAINT ***********************************************/
    /**********************************************************************************************/
    RealT tol(1e-8);
    ROL::PEBBL::IntegerConstraint<RealT> econ1;
    std::pair<int,RealT> Ip(2,0.0);
    econ1.add(Ip);
    ROL::Ptr<ROL::Vector<RealT>> j1 = econ1.makeConstraintVector(); 
    econ1.value(*j1,*x,tol);
    RealT j1norm = j1->norm();
    *outStream << "Constraint Violation: " << j1norm << std::endl;
    econ1.checkApplyJacobian(*x,*d,*j1,true,*outStream);

    *outStream << std::endl << std::endl;
    std::map<int,RealT> Im;
    Im.insert(std::pair<int,RealT>(3,1.0));
    Im.insert(std::pair<int,RealT>(9,1.0));
    econ1.add(Im);
    ROL::Ptr<ROL::Vector<RealT>> j2 = econ1.makeConstraintVector(); 
    econ1.value(*j2,*x,tol);
    RealT j2norm = j2->norm();
    *outStream << "Constraint Violation: " << j2norm << std::endl;
    econ1.checkApplyJacobian(*x,*d,*j2,true,*outStream);

    *outStream << std::endl << std::endl;
    ROL::PEBBL::IntegerConstraint<RealT> econ2(econ1);
    ROL::Ptr<ROL::Vector<RealT>> j3 = econ2.makeConstraintVector(); 
    econ2.value(*j3,*x,tol);
    RealT j3norm = j3->norm();
    *outStream << "Constraint Violation: " << j3norm << std::endl;
    econ2.checkApplyJacobian(*x,*d,*j3,true,*outStream);

    *outStream << std::endl << std::endl;
    errorFlag += (std::abs(j2norm-j3norm)<tol ? 0 : 1);
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
