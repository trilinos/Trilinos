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

/*! \file  test_09.cpp
    \brief Shows how to use the Constraint_Partitioned interface
           to solve Hock & Schittkowski's problem 39
*/

#include "ROL_HS39.hpp" 

#include "ROL_RandomVector.hpp"
#include "ROL_Constraint_Partitioned.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Algorithm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>

typedef double RealT;


int main(int argc, char *argv[]) {

  typedef std::vector<RealT>               vector;
  typedef ROL::Vector<RealT>               V;
  typedef ROL::StdVector<RealT>            SV;
  typedef ROL::Objective<RealT>            OBJ;
  typedef ROL::Constraint<RealT>           EC;  

  typedef typename vector::size_type       uint;

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

    uint xdim = 4;
    uint cdim = 1;

    ROL::Ptr<vector> x_exact_ptr = ROL::makePtr<vector>(xdim);
    (*x_exact_ptr)[0] = 1.0;
    (*x_exact_ptr)[1] = 1.0;

    ROL::Ptr<V> x     = ROL::makePtr<SV>( ROL::makePtr<vector>(xdim, 0.0) );
    ROL::Ptr<V> d     = ROL::makePtr<SV>( ROL::makePtr<vector>(xdim, 0.0) );
    ROL::Ptr<V> xtest = ROL::makePtr<SV>( ROL::makePtr<vector>(xdim, 0.0) );

    ROL::Ptr<V> c1    = ROL::makePtr<SV>( ROL::makePtr<vector>(cdim, 1.0) );
    ROL::Ptr<V> c2    = ROL::makePtr<SV>( ROL::makePtr<vector>(cdim, 1.0) );
    ROL::Ptr<V> l1    = ROL::makePtr<SV>( ROL::makePtr<vector>(cdim, 1.0) );
    ROL::Ptr<V> l2    = ROL::makePtr<SV>( ROL::makePtr<vector>(cdim, 1.0) );

    ROL::Ptr<V> c    =  ROL::CreatePartitionedVector( c1, c2 );
    ROL::Ptr<V> l    =  ROL::CreatePartitionedVector( l1, l2 );
  


    SV x_exact( x_exact_ptr );
 
    // Initial guess from H&S 39
    x->applyUnary(ROL::Elementwise::Fill<RealT>(2.0));
 
    ROL::RandomizeVector(*d, -1.0, 1.0 ); 
    ROL::RandomizeVector(*xtest, -1.0, 1.0 ); 
    
    ROL::Ptr<OBJ> obj  = ROL::makePtr<ROL::ZOO::Objective_HS39<RealT>>();
    ROL::Ptr<EC>  con1 = ROL::makePtr<ROL::ZOO::Constraint_HS39a<RealT>>();
    ROL::Ptr<EC>  con2 = ROL::makePtr<ROL::ZOO::Constraint_HS39b<RealT>>();
    std::vector<ROL::Ptr<EC> > cvec(2); cvec[0] = con1; cvec[1] = con2;
   
    ROL::Ptr<EC>  con = ROL::makePtr<ROL::Constraint_Partitioned<RealT>>(cvec);

    *outStream << "Checking objective" << std::endl;
    obj->checkGradient(*x,*d,true,*outStream); 

    *outStream << "\nChecking first equality constraint" << std::endl;
    con1->checkApplyJacobian( *xtest, *d, *c1 , true, *outStream );
    con1->checkApplyAdjointJacobian( *xtest, *l1, *c1, *d, true, *outStream );
    con1->checkApplyAdjointHessian( *xtest, *l1, *d, *xtest, true, *outStream );

    *outStream << "\nChecking second equality constraint" << std::endl;
    con2->checkApplyJacobian( *xtest, *d, *c2, true, *outStream );
    con2->checkApplyAdjointJacobian( *xtest, *l2, *c2, *d, true, *outStream );
    con2->checkApplyAdjointHessian( *xtest, *l2, *d, *xtest, true, *outStream );

    *outStream << "\nChecking partitioned equality constraint" << std::endl;
    con->checkApplyJacobian( *xtest, *d, *c, true, *outStream );
    con->checkApplyAdjointJacobian( *xtest, *l, *c, *d, true, *outStream );
    con->checkApplyAdjointHessian( *xtest, *l, *d, *xtest, true, *outStream );

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
    ROL::Algorithm<RealT> algo(stepname, parlist);

    algo.run(*x,x->dual(),*l,*c,*obj,*con,true,*outStream);

    x->axpy(-1.0,x_exact);

    if( x->norm() > 1e-6 ) {
      ++errorFlag;  
    }

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

