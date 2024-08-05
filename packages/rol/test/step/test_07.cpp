// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_07.cpp
    \brief Interior Point test using Hock & Schittkowski problem 32.
*/

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_HS32.hpp"
#include "ROL_Algorithm.hpp"


typedef double RealT;

int main(int argc, char *argv[]) {

  
   

  typedef std::vector<RealT>            vec;
  typedef ROL::StdVector<RealT>         SV;
  typedef ROL::Ptr<ROL::Vector<RealT> >      ROL::PtrV;

//  typedef ROL::PartitionedVector<RealT> PV;


  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag = 0;

  try {

    int xopt_dim  = 3; // Dimension of optimization vectors
    int ce_dim = 1;    // Dimension of equality constraint
    int ci_dim = 4;    // Dimension of inequality constraint

    // Exact solution
    ROL::Ptr<vec> x_exact_ptr = ROL::makePtr<vec>(xopt_dim,0.0);
    (*x_exact_ptr)[xopt_dim-1] = 1.0;

    ROL::Ptr<vec> xopt_ptr = ROL::makePtr<vec>(xopt_dim,0.0); // Optimization variables

    ROL::Ptr<vec> le_ptr  = ROL::makePtr<vec>(ce_dim,0.0);    // Equality multiplier
    ROL::Ptr<vec> li_ptr  = ROL::makePtr<vec>(ci_dim,0.0);    // Inequality multiplier
     
    // Feasible initial guess
    (*xopt_ptr)[0] = 0.1;
    (*xopt_ptr)[1] = 0.7;
    (*xopt_ptr)[2] = 0.2;

    ROL::PtrV xopt = ROL::makePtr<SV>(xopt_ptr);
    ROL::PtrV le  = ROL::makePtr<SV>(le_ptr);
    ROL::PtrV li  = ROL::makePtr<SV>(li_ptr);

    using ROL::ZOO::Objective_HS32;
    using ROL::ZOO::EqualityConstraint_HS32;
    using ROL::ZOO::InequalityConstraint_HS32;    

    ROL::Ptr<ROL::Objective<RealT> > obj_hs32 = ROL::makePtr<Objective_HS32<RealT>>(); 
    ROL::Ptr<ROL::EqualityConstraint<RealT> > eqcon_hs32 = ROL::makePtr<EqualityConstraint_HS32<RealT>>();
    ROL::Ptr<ROL::InequalityConstraint<RealT> > incon_hs32 = ROL::makePtr<InequalityConstraint_HS32<RealT>>();
    
    
    std::string stepname = "Interior Point";

    RealT mu = 0.1;            // Initial penalty parameter
    RealT factor = 0.1;        // Penalty reduction factor

    // Set solver parameters
    parlist->sublist("Step").sublist("Interior Point").set("Initial Barrier Penalty",mu);
    parlist->sublist("Step").sublist("Interior Point").set("Minimium Barrier Penalty",1e-8);
    parlist->sublist("Step").sublist("Interior Point").set("Barrier Penalty Reduction Factor",factor);
    parlist->sublist("Step").sublist("Interior Point").set("Subproblem Iteration Limit",30);

    parlist->sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Nominal Relative Tolerance",1.e-4);
    parlist->sublist("Step").sublist("Composite Step").sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist->sublist("Step").sublist("Composite Step").sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist->sublist("Step").sublist("Composite Step").sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist->sublist("Step").sublist("Composite Step").set("Output Level",0);

    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist->sublist("Status Test").set("Constraint Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Iteration Limit",100);

    ROL::OptimizationProblem<RealT> problem( obj_hs32, xopt, eqcon_hs32, le, incon_hs32, li, parlist);  

    // Define algorithm.
    ROL::Ptr<ROL::Algorithm<RealT> > algo;    
    algo = ROL::makePtr<ROL::Algorithm<RealT>>(stepname,*parlist);

    algo->run(problem,true,*outStream);   
  
    *outStream << std::endl << std::setw(20) << "Computed Minimizer" << std::setw(20) << "Exact Minimizer" << std::endl;
    for( int i=0;i<xopt_dim;++i ) {   
      *outStream << std::setw(20) << (*xopt_ptr)[i] << std::setw(20) << (*x_exact_ptr)[i] << std::endl;
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

