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
    ROL::Ptr<std::vector<RealT>> emul_ptr = ROL::makePtr<std::vector<RealT>>(1,0.0);
    ROL::Ptr<ROL::Vector<RealT>> emul     = ROL::makePtr<ROL::StdVector<RealT>>(emul_ptr);
    /**********************************************************************************************/
    /************************* CONSTRUCT OBJECTIVE FUNCTION ***************************************/
    /**********************************************************************************************/
    std::vector<RealT> alpha(N,0.0);
    *outStream << std::endl << std::endl;
    *outStream << "alpha =";
    for (int i = 0; i < N; ++i) {
      alpha[i] = static_cast<RealT>(rand())/static_cast<RealT>(RAND_MAX);
      *outStream << "  " << alpha[i];
    }
    *outStream << std::endl << std::endl;
    ROL::Ptr<ROL::Objective<RealT>> obj
      = ROL::makePtr<Objective_SimpleBinary<RealT>>(alpha);
    /**********************************************************************************************/
    /************************* CONSTRUCT CONSTRAINT ***********************************************/
    /**********************************************************************************************/
    int budget = 3;
    ROL::Ptr<ROL::Constraint<RealT>> econ
      = ROL::makePtr<Constraint_SimpleBinary<RealT>>(budget);
    /**********************************************************************************************/
    /************************* CONSTRUCT BOUND CONSTRAINT *****************************************/
    /**********************************************************************************************/
    ROL::Ptr<std::vector<RealT>> xl_ptr = ROL::makePtr<std::vector<RealT>>(N,0.0);
    ROL::Ptr<ROL::Vector<RealT>> xl     = ROL::makePtr<ROL::StdVector<RealT>>(xl_ptr);
    ROL::Ptr<std::vector<RealT>> xu_ptr = ROL::makePtr<std::vector<RealT>>(N,1.0);
    ROL::Ptr<ROL::Vector<RealT>> xu     = ROL::makePtr<ROL::StdVector<RealT>>(xu_ptr);
    ROL::Ptr<ROL::BoundConstraint<RealT>> bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(xl,xu);
    /**********************************************************************************************/
    /************************* SOLVE **************************************************************/
    /**********************************************************************************************/
    ROL::Ptr<ROL::Problem<RealT>>
      problem = ROL::makePtr<ROL::Problem<RealT>>(obj,x);
    problem->addBoundConstraint(bnd);
    problem->addLinearConstraint("Linear",econ,emul);
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);
    ROL::Solver<RealT> solver(problem,*parlist);
    *outStream << "Solve problem with no fixed binary variables"
               << std::endl << std::endl;
    clock_t start = clock();
    solver.solve(*outStream);
    *outStream << "Optimization time: "
               << (RealT)(clock()-start)/(RealT)CLOCKS_PER_SEC
               << " seconds." << std::endl;
    *outStream << std::endl;
    RealT sum(0);
    *outStream << "x =";
    for (int i = 0; i < N; ++i) {
      *outStream << "  " << (*x_ptr)[i];
      sum += (*x_ptr)[i];
    }
    *outStream << std::endl << std::endl;
    *outStream << "Sum(x) = " << sum << "  Budget = " << budget;
    *outStream << std::endl << std::endl;
    /**********************************************************************************************/
    /************************* ADD BINARY CONSTRAINTS AND SOLVE ***********************************/
    /**********************************************************************************************/
    std::map<int,RealT> fixed;
    fixed.insert(std::pair<int,RealT>(2,0.0));
    fixed.insert(std::pair<int,RealT>(5,0.0));
    fixed.insert(std::pair<int,RealT>(3,1.0));
    fixed.insert(std::pair<int,RealT>(9,1.0));
    ROL::Ptr<ROL::PEBBL::IntegerConstraint<RealT>> econ_bin
      = ROL::makePtr<ROL::PEBBL::IntegerConstraint<RealT>>();
    econ_bin->add(fixed);
    ROL::Ptr<ROL::Vector<RealT>> emul_bin
      = econ_bin->makeConstraintVector();
    problem->edit();
    problem->addLinearConstraint("Integer",econ_bin,emul_bin);
    problem->finalize(false,true,*outStream);
    problem->check(true,*outStream);
    ROL::Solver<RealT> solver_bin(problem,*parlist);
    *outStream << "Solve problem with {2,5} set to 0 and {3,9} set to 1"
               << std::endl << std::endl;
    clock_t start_bin = clock();
    solver_bin.solve(*outStream);
    *outStream << "Optimization time: "
               << (RealT)(clock()-start_bin)/(RealT)CLOCKS_PER_SEC
               << " seconds." << std::endl;
    *outStream << std::endl;
    RealT sum_bin(0);
    *outStream << "x =";
    for (int i = 0; i < N; ++i) {
      *outStream << "  " << (*x_ptr)[i];
      sum_bin += (*x_ptr)[i];
    }
    *outStream << std::endl << std::endl;
    *outStream << "Sum(x) = " << sum_bin << "  Budget = " << budget;
    *outStream << std::endl << std::endl;

    errorFlag += (std::abs((*x_ptr)[2]-0.0)<std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 0 : 1);
    errorFlag += (std::abs((*x_ptr)[5]-0.0)<std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 0 : 1);
    errorFlag += (std::abs((*x_ptr)[3]-1.0)<std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 0 : 1);
    errorFlag += (std::abs((*x_ptr)[9]-1.0)<std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 0 : 1);
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
