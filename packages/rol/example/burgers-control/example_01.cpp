// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve an optimal control problem constrained by
           steady Burgers' equation with bound constraints.
*/

#include "ROL_TypeB_PrimalDualActiveSetAlgorithm.hpp"
#include "ROL_TypeB_LinMoreAlgorithm.hpp"
#include "ROL_Bounds.hpp"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>

#include "ROL_Stream.hpp"

#include "example_01.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>    vector;
  typedef ROL::Vector<RealT>    V;
  typedef ROL::StdVector<RealT> SV;

  typedef typename vector::size_type uint;

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

  // *** Example body.

  try {
    // Initialize objective function.
    uint nx     = 1028;  // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    Objective_BurgersControl<RealT> obj(alpha,nx);
    // Initialize iteration vectors.
    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(nx+2, 1.0);
    ROL::Ptr<vector> y_ptr = ROL::makePtr<vector>(nx+2, 0.0);
    for (uint i=0; i<nx+2; i++) {
      (*x_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
      (*y_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }

    SV x(x_ptr);
    SV y(y_ptr);

    // Check derivatives.
    obj.checkGradient(x,x,y,true,*outStream);
    obj.checkHessVec(x,x,y,true,*outStream);

    // Initialize Constraints
    ROL::Ptr<vector> l_ptr = ROL::makePtr<vector>(nx+2,0.0);
    ROL::Ptr<vector> u_ptr = ROL::makePtr<vector>(nx+2,1.0);
    ROL::Ptr<V> lo = ROL::makePtr<SV>(l_ptr);
    ROL::Ptr<V> up = ROL::makePtr<SV>(u_ptr);

    ROL::Bounds<RealT> bcon(lo,up);

    // Primal dual active set.
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

    // Krylov parameters.
    parlist->sublist("General").sublist("Krylov").set("Absolute Tolerance",1.e-8);
    parlist->sublist("General").sublist("Krylov").set("Relative Tolerance",1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit",50);
    // PDAS parameters.
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-10);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-8);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit", 10);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",(alpha>0.0)?alpha:1.e-4);
    // Status test parameters.
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-16);
    parlist->sublist("Status Test").set("Iteration Limit",100);
    // Set initial guess.
    x.zero();
    {
      // Define algorithm.
      ROL::TypeB::PrimalDualActiveSetAlgorithm<RealT> algo(*parlist);
      // Run algorithm.
      algo.run(x, obj, bcon, *outStream);
    }
    // Output control to file.
    std::ofstream file_pdas;
    file_pdas.open("control_PDAS.txt");
    for ( unsigned i = 0; i < (unsigned)nx+2; i++ ) {
      file_pdas << (*x_ptr)[i] << "\n";
    }
    file_pdas.close();

    // Projected Newton.
    parlist->sublist("General").sublist("Krylov").set("Absolute Tolerance",1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Relative Tolerance",1.e-2);
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit",50);
    // Set initial guess.
    y.zero();
    {
      // Define algorithm.
      ROL::TypeB::LinMoreAlgorithm<RealT> algo(*parlist);
      // Run Algorithm
      algo.run(y,obj,bcon,*outStream);
    }
    // Output control to file.
    std::ofstream file_tr;
    file_tr.open("control_TR.txt");
    for ( unsigned i = 0; i < (unsigned)nx+2; i++ ) {
      file_tr << (*y_ptr)[i] << "\n";
    }
    file_tr.close();
    // Output state to file.
    std::vector<RealT> u(nx,0.0);
    std::vector<RealT> param(4,0.0);
    obj.solve_state(u,*x_ptr,param);
    std::ofstream file;
    file.open("state.txt");
    for (unsigned i=0; i<(unsigned)nx; i++) {
      file << i/((RealT)(nx+1)) << "  " << u[i] << "\n";
    }
    file.close();

    // Compute error between PDAS and Lin-More solutions.
    ROL::Ptr<ROL::Vector<RealT> > diff = x.clone();
    diff->set(x);
    diff->axpy(-1.0,y);
    RealT error = diff->norm();
    *outStream << "\nError between PDAS solution and TR solution is " << error << "\n";
    errorFlag = ((error > 1e2*std::sqrt(ROL::ROL_EPSILON<RealT>())) ? 1 : 0);
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

