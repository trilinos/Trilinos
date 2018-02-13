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

/*! \file  example_01.cpp
    \brief Shows how to solve an optimal control problem constrained by 
           steady Burgers' equation with bound constraints.
*/

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
  Teuchos::oblackholestream bhs; // outputs nothing
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

    // Check deriatives.
    obj.checkGradient(x,x,y,true,*outStream);
    obj.checkHessVec(x,x,y,true,*outStream);

    // Initialize Constraints
    ROL::Ptr<vector> l_ptr = ROL::makePtr<vector>(nx+2,0.0);
    ROL::Ptr<vector> u_ptr = ROL::makePtr<vector>(nx+2,1.0);
    ROL::Ptr<V> lo = ROL::makePtr<SV>(l_ptr);
    ROL::Ptr<V> up = ROL::makePtr<SV>(u_ptr); 
      
    ROL::Bounds<RealT> icon(lo,up);

    // ROL components.
    ROL::Ptr<ROL::Algorithm<RealT> > algo;

    // Primal dual active set.
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
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
    // Define algorithm.
    algo = ROL::makePtr<ROL::Algorithm<RealT>>("Primal Dual Active Set",*parlist,false);
    // Run algorithm.
    x.zero();
    algo->run(x, obj, icon, true, *outStream);
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
    // Define algorithm.
    algo = ROL::makePtr<ROL::Algorithm<RealT>>("Trust Region",*parlist,false);
    // Run Algorithm
    y.zero();
    algo->run(y,obj,icon,true,*outStream);
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
    // Compute error 
    ROL::Ptr<ROL::Vector<RealT> > diff = x.clone();
    diff->set(x);
    diff->axpy(-1.0,y);
    RealT error = diff->norm();
    *outStream << "\nError between PDAS solution and TR solution is " << error << "\n";
    errorFlag = ((error > 1e2*std::sqrt(ROL::ROL_EPSILON<RealT>())) ? 1 : 0);
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

