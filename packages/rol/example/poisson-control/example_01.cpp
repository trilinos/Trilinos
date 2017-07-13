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
    \brief Shows how to solve the inverse Poisson problem using trust-region
           methods with dense Hessian diagnostics.
*/

#define USE_HESSVEC 1

#include "ROL_PoissonControl.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_BoundConstraint.hpp"


template <class Real>
class StatusTest_PDAS : public ROL::StatusTest<Real> {
private:

  Real gtol_;
  Real stol_;
  int  max_iter_;

public:

  virtual ~StatusTest_PDAS() {}

  StatusTest_PDAS( Real gtol = 1.e-6, Real stol = 1.e-12, int max_iter = 100 ) :
    gtol_(gtol), stol_(stol), max_iter_(max_iter) {}

  /** \brief Check algorithm status.
  */
  virtual bool check( ROL::AlgorithmState<Real> &state ) {
     if ( (state.gnorm > this->gtol_) &&
          (state.snorm > this->stol_) &&
          (state.iter  < this->max_iter_) ) {
       return true;
     }
     else {
       if ( state.iter < 2 ) {
         return true;
       }
       return false;
     }
  }

};

typedef double RealT;

int main(int argc, char *argv[]) {

  typedef std::vector<RealT>     vector;
  typedef ROL::Vector<RealT>     V;
  typedef ROL::StdVector<RealT>  SV;
 
  typedef typename vector::size_type uint;
  
  using Teuchos::RCP;  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.

  try {
    uint dim = 256; // Set problem dimension.
    RealT alpha = 1.e-4;
    ROL::ZOO::Objective_PoissonControl<RealT> obj(alpha);

    RCP<vector> l_rcp = rcp( new vector(dim) );
    RCP<vector> u_rcp = rcp( new vector(dim) );

    RCP<V> lo = rcp( new SV(l_rcp) );
    RCP<V> up = rcp( new SV(u_rcp) );

    for ( uint i = 0; i < dim; i++ ) {
      if ( i < dim/3.0  ||  i > 2*dim/3.0 ) {
        (*l_rcp)[i] = 0.0; 
        (*u_rcp)[i] = 0.25;
      }
      else {
        (*l_rcp)[i] = 0.75;
        (*u_rcp)[i] = 1.0;
      }
    }
    ROL::BoundConstraint<RealT> icon(lo,up);

    // Primal dual active set.
    std::string filename = "input.xml";
    RCP<Teuchos::ParameterList> parlist = rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Krylov parameters.
    parlist->sublist("General").sublist("Krylov").set("Absolute Tolerance",1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Relative Tolerance",1.e-2);
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit",50);

    // PDAS parameters.
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-8);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-6);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit", 1);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",(alpha>0.0)?alpha:1.e-4);

    // Status test parameters.
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Iteration Limit",100);

    // Define algorithm.
    RCP<ROL::Algorithm<RealT> > algo = rcp(new ROL::Algorithm<RealT>("Primal Dual Active Set",*parlist,false));

    // Iteration vector.
    RCP<vector> x_rcp = rcp( new vector(dim, 0.0) );
    SV x(x_rcp);

    // Run algorithm.
    x.zero();
    algo->run(x, obj, icon, true, *outStream);
    std::ofstream file;
    file.open("control_PDAS.txt");

    for ( uint i = 0; i < dim; i++ ) {
      file << (*x_rcp)[i] << "\n";
    }
    file.close();

    // Projected Newton.
    // re-load parameters
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    // Set algorithm.
    algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
    // Iteration vector.
    RCP<vector> y_rcp = rcp( new vector(dim, 0.0) );
    SV y(y_rcp);

    // Run Algorithm
    y.zero();
    algo->run(y, obj, icon, true, *outStream);

    std::ofstream file_tr;
    file_tr.open("control_TR.txt");
    for ( uint i = 0; i < dim; i++ ) {
      file_tr << (*y_rcp)[i] << "\n";
    }
    file_tr.close();
   
    RCP<V> error = x.clone();
    error->set(x);
    error->axpy(-1.0,y);
    *outStream << "\nError between PDAS solution and TR solution is " << error->norm() << "\n";
    errorFlag = ((error->norm() > 1e2*std::sqrt(ROL::ROL_EPSILON<RealT>())) ? 1 : 0);
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

