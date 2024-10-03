// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_01.hpp"
#include "ROL_PEBBL_StdIntegerTransformation.hpp"
#include "ROL_PEBBL_BuildTransformation.hpp"

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
    ROL::Ptr<std::vector<RealT>> xbin_ptr = ROL::makePtr<std::vector<RealT>>(N,0.0);
    ROL::Ptr<ROL::Vector<RealT>> xbin     = ROL::makePtr<ROL::StdVector<RealT>>(xbin_ptr);
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
    ROL::Ptr<ROL::Objective<RealT> > obj
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
    ROL::Ptr<ROL::PEBBL::StdIntegerTransformation<RealT>> trans_bin
      = ROL::makePtr<ROL::PEBBL::StdIntegerTransformation<RealT>>();
    std::map<int,RealT> fixed;
    fixed.insert(std::pair<int,RealT>(2,0.0));
    fixed.insert(std::pair<int,RealT>(5,0.0));
    fixed.insert(std::pair<int,RealT>(3,1.0));
    fixed.insert(std::pair<int,RealT>(9,1.0));
    trans_bin->add(fixed);
    RealT tol(1e-8);
    trans_bin->value(*xbin,*x,tol);
    ROL::Ptr<ROL::PEBBL::BuildTransformation<RealT>> build
      = ROL::makePtr<ROL::PEBBL::BuildTransformation<RealT>>(trans_bin,x);
    ROL::Ptr<ROL::Objective<RealT>>   obj_trans = build->transform(obj);
    ROL::Ptr<ROL::Constraint<RealT>> econ_trans = build->transform(econ);
    ROL::Ptr<ROL::Problem<RealT>>
      problem_bin = ROL::makePtr<ROL::Problem<RealT>>(obj_trans,xbin);
    problem_bin->addBoundConstraint(bnd);
    problem_bin->addLinearConstraint("Linear",econ_trans,emul);
    problem_bin->finalize(false,true,*outStream);
    problem_bin->check(true,*outStream);
    ROL::Solver<RealT> solver_bin(problem_bin,*parlist);
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
      *outStream << "  " << (*xbin_ptr)[i];
      sum_bin += (*xbin_ptr)[i];
    }
    *outStream << std::endl << std::endl;
    *outStream << "Sum(x) = " << sum_bin << "  Budget = " << budget;
    *outStream << std::endl << std::endl;

    errorFlag += (std::abs((*xbin_ptr)[2]-0.0)<std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 0 : 1);
    errorFlag += (std::abs((*xbin_ptr)[5]-0.0)<std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 0 : 1);
    errorFlag += (std::abs((*xbin_ptr)[3]-1.0)<std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 0 : 1);
    errorFlag += (std::abs((*xbin_ptr)[9]-1.0)<std::sqrt(ROL::ROL_EPSILON<RealT>()) ? 0 : 1);
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
