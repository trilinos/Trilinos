// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_11.cpp
    \brief Interior Point test using Hock & Schittkowski problem 24.
*/

#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_HS24.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_ObjectiveFromBoundConstraint.hpp"

int main(int argc, char *argv[]) {

  
   

  typedef double RealT;

  typedef ROL::Vector<RealT>               V;
  typedef ROL::BoundConstraint<RealT>      BC;
  typedef ROL::Objective<RealT>            OBJ;
  typedef ROL::InequalityConstraint<RealT> INEQ; 

   

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

    ROL::ZOO::getHS24<RealT> HS24;
    ROL::Ptr<V>    x     = HS24.getInitialGuess();
    ROL::Ptr<V>    xs    = HS24.getSolution();
    ROL::Ptr<V>    inmul = HS24.getInequalityMultiplier();    

    ROL::Ptr<BC>   bnd   = HS24.getBoundConstraint();
    ROL::Ptr<OBJ>  obj   = HS24.getObjective();
    ROL::Ptr<INEQ> incon = HS24.getInequalityConstraint();
    ROL::Ptr<BC>   inbnd = HS24.getSlackBoundConstraint();
   
    

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

    // Define Optimization Problem 
    ROL::OptimizationProblem<RealT> problem( obj, x, bnd, incon, inmul, inbnd );

    ROL::Ptr<V> d = x->clone();
    RandomizeVector(*d);

//    problem.checkObjectiveGradient(*d); 
//    problem.checkObjectiveHessVec(*d); 

    // Define algorithm.
    ROL::Ptr<ROL::Algorithm<RealT> > algo;    
    algo = ROL::makePtr<ROL::Algorithm<RealT>>(stepname,*parlist);

    algo->run(problem,true,*outStream);   

    x->axpy(-1.0,*xs);

    if( x->norm()>= 1e-4 )
    {
      ++errorFlag;
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
