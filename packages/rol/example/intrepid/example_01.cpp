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

#include<iostream>

#include "example_01.hpp"

#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_ConstraintStatusTest.hpp"

// Minimize (1/2)*||u-utarget||^2+(gamma/2)||z||^2
// 
// subject to the BVP constraint 
//
// -u"+(zu)^2 = f
// u'(0)=u'(L)=0
//
// The target function is x^2*(L-x)^2

using namespace ROL;

typedef double             RealT;
typedef Vector<RealT>      V;
typedef StdVector<RealT>   SV;
typedef std::vector<RealT> vec;


int main(int argc, char *argv[]) {

  
  
     
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    ROL::makePtrFromRef(std::cout);
  else
    ROL::makePtrFromRef(bhs);

  int errorFlag = 0;
 
  try {

    int numCells = 20;
    int numFields = 2;
    RealT domainLength = 1.0;
    RealT h = domainLength/RealT(numCells);
    RealT gamma = 1e-4;

    int nDoF = numCells*(numFields-1)+1;

    // Create discretization
    ROL::Ptr<Discretization<RealT>> disc = ROL::makePtr<Discretization<RealT>>(numCells,numFields,domainLength);
  
    ROL::Ptr<vec> u_ptr   = ROL::makePtr<vec>(nDoF,1.0);      // Simulation vector 
    ROL::Ptr<vec> z_ptr   = ROL::makePtr<vec>(nDoF,1.0);      // Optimization vector 
    ROL::Ptr<vec> yu_ptr  = ROL::makePtr<vec>(nDoF,0.0);      // Test vector in U
    ROL::Ptr<vec> yz_ptr  = ROL::makePtr<vec>(nDoF,0.0);      // Test vector in Z

    ROL::Ptr<vec> gu_ptr  = ROL::makePtr<vec>(nDoF,0.0);      // Gradient w.r.t. Sim vector
    ROL::Ptr<vec> gz_ptr  = ROL::makePtr<vec>(nDoF,0.0);      // Gradient w.r.t. Opt vector

    ROL::Ptr<vec> utarget_ptr = ROL::makePtr<vec>(nDoF,1.0);  // Target vector

    ROL::Ptr<vec> v_ptr   = ROL::makePtr<vec>(nDoF,1.0);      // Simulation vector 
    ROL::Ptr<vec> w_ptr   = ROL::makePtr<vec>(nDoF,1.0);      // Optimization vector 

    ROL::Ptr<vec> c_ptr  = ROL::makePtr<vec>(nDoF,0.0);       // Constraint vector
    ROL::Ptr<vec> l_ptr  = ROL::makePtr<vec>(nDoF,0.0);       // Lagrange multiplier

    // -----------------------
    // Begin derivative checks
    // -----------------------

    RealT left = -1e0, right = 1e0; 
    for(int i=0;i<nDoF;++i) {

      (*v_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*w_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*yu_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;
      (*yz_ptr)[i] = ( (RealT)rand() / (RealT)RAND_MAX ) * (right - left) + left;

      RealT x = i*h; // Grid points

      (*utarget_ptr)[i] = x*x*(domainLength-x)*(domainLength-x);
    }

    // Make ROL::StdVector 
    SV u(u_ptr);
    SV z(z_ptr);
    SV gu(u_ptr);
    SV gz(z_ptr);
    SV yu(yu_ptr);
    SV yz(yz_ptr);
    SV v(v_ptr);
    SV w(w_ptr);
    SV c(c_ptr); 
    SV l(l_ptr);

    ROL::Ptr<V> utarget = ROL::makePtr<SV>(utarget_ptr); 

    ROL::Ptr<V> up   = &u,false;
    ROL::Ptr<V> zp   = &z,false;
    ROL::Ptr<V> gup  = &gu,false;
    ROL::Ptr<V> gzp  = &gz,false;
    ROL::Ptr<V> yup  = &yu,false;
    ROL::Ptr<V> yzp  = &yz,false;

    Vector_SimOpt<RealT> uz(up,zp);
    Vector_SimOpt<RealT> g(gup,gzp);
    Vector_SimOpt<RealT> y(yup,yzp);

    // Tracking Objective
    ROL::Ptr<Objective_SimOpt<RealT>> obj = ROL::makePtr<TrackingObjective<RealT>>(disc,utarget,gamma);

    // Constraint
    ROL::Ptr<Constraint_SimOpt<RealT>> con = ROL::makePtr<BVPConstraint<RealT>>(disc);
 
    obj->checkGradient(uz,y,true,*outStream);
    obj->checkHessVec(uz,y,true,*outStream);

    con->checkApplyJacobian(uz,y,c,true,*outStream);
    con->checkApplyAdjointHessian(uz,yz,y,uz,true,*outStream);

    con->checkInverseJacobian_1(c,yu,u,z,true,*outStream);
    con->checkInverseAdjointJacobian_1(c,yu,u,z,true,*outStream);

    con->checkAdjointConsistencyJacobian_1(w,v,u,z,true,*outStream);
    con->checkAdjointConsistencyJacobian_2(w,v,u,z,true,*outStream);
 
    // --------------------
    // End derivative checks 
    // --------------------

  
    // ----------------
    // Run optimization 
    // ----------------  

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
    parlist.sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist.sublist("Status Test").set("Iteration Limit",100);
    Algorithm<RealT> algo(stepname, parlist);

    // Run algorithm.
    algo.run(uz,g,l,c,*obj,*con,true,*outStream);

  } 
  catch ( std::logic_error err ) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

