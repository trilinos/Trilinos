// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test ConstraintFromObjective and ObjectiveFromConstraint
           which will help test ConstraintManager and OptimizationProblem
*/

#include "ROL_SimpleEqConstrained.hpp"

#include "ROL_ConstraintFromObjective.hpp"
#include "ROL_ObjectiveFromConstraint.hpp"

#include "ROL_OptimizationProblem.hpp"
#include "ROL_OptimizationSolver.hpp"


#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  using Teuchos::RCP; using Teuchos::rcp;

  using Obj     = ROL::Objective<RealT>;
  using Con     = ROL::Constraint<RealT>;
  using V       = ROL::Vector<RealT>;
  using StdV    = ROL::StdVector<RealT>;
  using ScalarV = ROL::SingletonVector<RealT>;

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  int errorFlag  = 0;

  RealT tol = std::sqrt(ROL::ROL_EPSILON<RealT>());

  // *** Test body.

  try {

    int xdim = 5;
    int cdim = 3;
    
    // Optimization vector
    RCP<V> x  = rcp( new StdV( rcp( new std::vector<RealT>(xdim) ) ) );
    RCP<V> c  = rcp( new StdV( rcp( new std::vector<RealT>(cdim) ) ) );
    RCP<V> e0 = c->basis(0);
    RCP<V> e1 = c->basis(1);
    RCP<V> e2 = c->basis(2);
 
    // Exact solution
    RCP<V> sol   = x->clone();   
    RCP<V> error = x->clone();

    RCP<Obj> obj = Teuchos::null; 
    RCP<Con> con = Teuchos::null;
    
    ROL::ZOO::getSimpleEqConstrained<RealT,StdV,StdV,StdV,StdV>( obj, con, *x, *sol );

    error->set(*sol);

    // Extract constraint components to make objectives
    RCP<Obj> obj_0 = rcp( new ROL::ObjectiveFromConstraint<RealT>( con, *e0 ) );
    RCP<Obj> obj_1 = rcp( new ROL::ObjectiveFromConstraint<RealT>( con, *e1 ) );
    RCP<Obj> obj_2 = rcp( new ROL::ObjectiveFromConstraint<RealT>( con, *e2 ) );

    // Create separate constraints from the objectives
    RCP<Con> con_0 = rcp( new ROL::ConstraintFromObjective<RealT>( obj_0 ) );
    RCP<Con> con_1 = rcp( new ROL::ConstraintFromObjective<RealT>( obj_1 ) );
    RCP<Con> con_2 = rcp( new ROL::ConstraintFromObjective<RealT>( obj_2 ) );
    
    std::vector<RCP<Con> > con_array;
    con_array.push_back(con_0);
    con_array.push_back(con_1);
    con_array.push_back(con_2);

    // Lagrange multipliers
    RCP<V> l0 = rcp( new ScalarV(0) );
    RCP<V> l1 = rcp( new ScalarV(0) );
    RCP<V> l2 = rcp( new ScalarV(0) );
  
    std::vector<RCP<V> > l_array;
    l_array.push_back(l0);
    l_array.push_back(l1);
    l_array.push_back(l2);
   
    ROL::OptimizationProblem<RealT> opt( obj,             // Objective
                                         x,               // Optimization vector
                                         Teuchos::null,   // No bound constraint
                                         con_array,       // Array of scalar equality constraints
                                         l_array);        // Array of scalar lagrange multipliers
 
    opt.check(*outStream);

    // Define algorithm.
    Teuchos::ParameterList parlist;
    std::string stepname = "Composite Step";
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Nominal Relative Tolerance",1.e-4);
    parlist.sublist("Step").sublist(stepname).sublist("Optimality System Solver").set("Fix Tolerance",true);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Iteration Limit",20);
    parlist.sublist("Step").sublist(stepname).sublist("Tangential Subproblem Solver").set("Relative Tolerance",1e-2);
    parlist.sublist("Step").sublist(stepname).set("Output Level",0);
    parlist.sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Constraint Tolerance",1.e-12);
    parlist.sublist("Status Test").set("Step Tolerance",1.e-18);
    parlist.sublist("Status Test").set("Iteration Limit",100);

    ROL::OptimizationSolver<RealT> solver( opt, parlist );

    solver.solve( *outStream );

    error->axpy(-1.0,*x);
    RealT error_norm = error->norm();
    
    *outStream << "\n\n Relative norm of final optimization vector error: " << error_norm << std::endl;
    
    if(error_norm>tol) 
      ++errorFlag;

  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return 0;

}

