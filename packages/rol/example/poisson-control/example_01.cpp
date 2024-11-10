// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the inverse Poisson problem using trust-region
           methods with dense Hessian diagnostics.
*/

#define USE_HESSVEC 1

#include "ROL_Bounds.hpp"
#include "ROL_PoissonControl.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Types.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>
#include <algorithm>



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
    uint dim = 256; // Set problem dimension.
    RealT alpha = 1.e-4;
    ROL::ZOO::Objective_PoissonControl<RealT> obj(alpha);

    ROL::Ptr<vector> l_ptr = ROL::makePtr<vector>(dim);
    ROL::Ptr<vector> u_ptr = ROL::makePtr<vector>(dim);

    ROL::Ptr<V> lo = ROL::makePtr<SV>(l_ptr);
    ROL::Ptr<V> up = ROL::makePtr<SV>(u_ptr);

    for ( uint i = 0; i < dim; i++ ) {
      if ( i < dim/3.0  ||  i > 2*dim/3.0 ) {
        (*l_ptr)[i] = 0.0; 
        (*u_ptr)[i] = 0.25;
      }
      else {
        (*l_ptr)[i] = 0.75;
        (*u_ptr)[i] = 1.0;
      }
    }
    ROL::Bounds<RealT> icon(lo,up);

    // Primal dual active set.
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

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
    ROL::Ptr<ROL::Step<RealT>>       step   = ROL::makePtr<ROL::PrimalDualActiveSetStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>> status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Ptr<ROL::Algorithm<RealT>>  algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);

    // Iteration vector.
    ROL::Ptr<vector> x_ptr = ROL::makePtr<vector>(dim, 0.0);
    SV x(x_ptr);

    // Run algorithm.
    x.zero();
    algo->run(x, obj, icon, true, *outStream);
    std::ofstream file;
    file.open("control_PDAS.txt");

    for ( uint i = 0; i < dim; i++ ) {
      file << (*x_ptr)[i] << "\n";
    }
    file.close();

    // Projected Newton.
    // re-load parameters
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    // Set algorithm.
    step   = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    algo   = ROL::makePtr<ROL::Algorithm<RealT>>(step,status,false);
    // Iteration vector.
    ROL::Ptr<vector> y_ptr = ROL::makePtr<vector>(dim, 0.0);
    SV y(y_ptr);

    // Run Algorithm
    y.zero();
    algo->run(y, obj, icon, true, *outStream);

    std::ofstream file_tr;
    file_tr.open("control_TR.txt");
    for ( uint i = 0; i < dim; i++ ) {
      file_tr << (*y_ptr)[i] << "\n";
    }
    file_tr.close();
   
    ROL::Ptr<V> error = x.clone();
    error->set(x);
    error->axpy(-1.0,y);
    *outStream << "\nError between PDAS solution and TR solution is " << error->norm() << "\n";
    errorFlag = ((error->norm() > 1e2*std::sqrt(ROL::ROL_EPSILON<RealT>())) ? 1 : 0);
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

