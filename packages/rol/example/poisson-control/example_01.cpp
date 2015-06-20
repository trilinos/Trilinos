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
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Types.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_BoundConstraint.hpp"

template<class Real>
class BoundConstraint_PoissonControl : public ROL::BoundConstraint<Real> {
private:
  int dim_;
  std::vector<Real> x_lo_;
  std::vector<Real> x_up_;
  Real min_diff_;
public:
  BoundConstraint_PoissonControl(std::vector<Real> &lo, std::vector<Real> &up) {
    dim_ = lo.size();
    x_lo_.clear();
    x_lo_.assign(lo.begin(),lo.end());
    x_up_.clear();
    x_up_.assign(up.begin(),up.end());
    for ( unsigned i = 0; i < (unsigned)dim_; i++ ) {
      if ( i == 0 ) {
        min_diff_ = x_up_[i]-x_lo_[i];
      }
      else {
        min_diff_ = std::min(min_diff_,x_up_[i]-x_lo_[i]);
      }
    }
    min_diff_ *= 0.5;
  }
  bool isFeasible( const ROL::Vector<Real> &x ) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    bool val = true;
    int  cnt = 1;
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( (*ex)[i] >= this->x_lo_[i] && (*ex)[i] <= this->x_up_[i] ) { cnt *= 1; }
      else                                                            { cnt *= 0; }
    }
    if ( cnt == 0 ) { val = false; }
    return val;
  }
  void project( ROL::Vector<Real> &x ) {
    Teuchos::RCP<std::vector<Real> > ex =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(x)).getVector());
    for ( int i = 0; i < this->dim_; i++ ) {
      (*ex)[i] = std::max(this->x_lo_[i],std::min(this->x_up_[i],(*ex)[i]));
    }
  }
  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ) {
        (*ev)[i] = 0.0;
      }
    }
  }
  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn) ) {
        (*ev)[i] = 0.0;
      }
    }
  }
  void pruneLowerActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) ){
        (*ev)[i] = 0.0;
      }
    }
  }
  void pruneUpperActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
    Teuchos::RCP<const std::vector<Real> > ex =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<const std::vector<Real> > eg =
      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(g))).getVector();
    Teuchos::RCP<std::vector<Real> > ev =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
    Real epsn = std::min(eps,this->min_diff_);
    for ( int i = 0; i < this->dim_; i++ ) {
      if ( ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
        (*ev)[i] = 0.0;
      }
    }
  }
  void setVectorToUpperBound( ROL::Vector<Real> &u ) {
    Teuchos::RCP<std::vector<Real> > us = Teuchos::rcp( new std::vector<Real>(this->dim_,0.0) );
    us->assign(this->x_up_.begin(),this->x_up_.end()); 
    Teuchos::RCP<ROL::Vector<Real> > up = Teuchos::rcp( new ROL::StdVector<Real>(us) );
    u.set(*up);
  }
  void setVectorToLowerBound( ROL::Vector<Real> &l ) {
    Teuchos::RCP<std::vector<Real> > ls = Teuchos::rcp( new std::vector<Real>(this->dim_,0.0) );
    ls->assign(this->x_lo_.begin(),this->x_lo_.end()); 
    Teuchos::RCP<ROL::Vector<Real> > lp = Teuchos::rcp( new ROL::StdVector<Real>(ls) );
    l.set(*lp);
  }
//  void pruneActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real eps) {
//    Teuchos::RCP<const std::vector<Real> > ex =
//      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
//    Teuchos::RCP<std::vector<Real> > ev =
//      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
//    Real epsn = std::min(eps,this->min_diff_);
//    for ( int i = 0; i < this->dim_; i++ ) {
//      if ( ((*ex)[i] <= this->x_lo_[i]+epsn) ||
//           ((*ex)[i] >= this->x_up_[i]-epsn) ) {
//        (*ev)[i] = 0.0;
//      }
//    }
//  }
//  void pruneActive(ROL::Vector<Real> &v, const ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real eps) {
//    Teuchos::RCP<const std::vector<Real> > ex =
//      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
//    Teuchos::RCP<const std::vector<Real> > eg =
//      (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(g))).getVector();
//    Teuchos::RCP<std::vector<Real> > ev =
//      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(v)).getVector());
//    Real epsn = std::min(eps,this->min_diff_);
//    for ( int i = 0; i < this->dim_; i++ ) {
//      if ( ((*ex)[i] <= this->x_lo_[i]+epsn && (*eg)[i] > 0.0) ||
//           ((*ex)[i] >= this->x_up_[i]-epsn && (*eg)[i] < 0.0) ) {
//        (*ev)[i] = 0.0;
//      }
//    }
//  }
};

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
    int dim = 256; // Set problem dimension.
    RealT alpha = 1.e-4;
    ROL::ZOO::Objective_PoissonControl<RealT> obj(alpha);
    std::vector<RealT> lo(dim);
    std::vector<RealT> up(dim);
    for ( unsigned i = 0; i < (unsigned)dim; i++ ) {
      if ( i < (unsigned)dim/3  ||  i > 2*(unsigned)dim/3 ) {
        lo[i] = 0.0; 
        up[i] = 0.25;
      }
      else {
        lo[i] = 0.75;
        up[i] = 1.0;
      }
    }
    BoundConstraint_PoissonControl<RealT> icon(lo,up);

    Teuchos::ParameterList parlist;
    // Krylov parameters.
    parlist.set("Absolute Krylov Tolerance",              1.e-4);
    parlist.set("Relative Krylov Tolerance",              1.e-2);
    parlist.set("Maximum Number of Krylov Iterations",    50);
    // PDAS parameters.
    parlist.set("PDAS Relative Step Tolerance",           1.e-8);
    parlist.set("PDAS Relative Gradient Tolerance",       1.e-6);
    parlist.set("PDAS Maximum Number of Iterations",      1);
    parlist.set("PDAS Dual Scaling",                      alpha);      
    // Define step.
    ROL::PrimalDualActiveSetStep<RealT> step_pdas(parlist);

    // Define status test.
    RealT gtol  = 1e-12;  // norm of gradient tolerance
    RealT stol  = 1e-14;  // norm of step tolerance
    int   maxit = 100;    // maximum number of iterations
    StatusTest_PDAS<RealT> status_pdas(gtol, stol, maxit);    

    // Define algorithm.
    ROL::DefaultAlgorithm<RealT> algo_pdas(step_pdas,status_pdas,false);

    // Iteration vector.
    Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    // Set initial guess.
    for (int i=0; i<dim; i++) {
      (*x_rcp)[i] = 0.0;
    }
    ROL::StdVector<RealT> x(x_rcp);

    // Run algorithm.
    algo_pdas.run(x, obj, icon, true);

    std::ofstream file;
    file.open("control_PDAS.txt");
    for ( unsigned i = 0; i < (unsigned)dim; i++ ) {
      file << (*x_rcp)[i] << "\n";
    }
    file.close();

    // Use Projected Methods
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist_tr = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist_tr) );
    ROL::TrustRegionStep<RealT> step_tr(*parlist_tr);
    ROL::StatusTest<RealT> status_tr(gtol,stol,maxit);
    ROL::DefaultAlgorithm<RealT> algo_tr(step_tr,status_tr,false);
    // Iteration vector.
    Teuchos::RCP<std::vector<RealT> > y_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
    // Set initial guess.
    for (int i=0; i<dim; i++) {
      (*y_rcp)[i] = 0.0;
    }
    ROL::StdVector<RealT> y(y_rcp);
    // Run Algorithm
    algo_tr.run(y,obj,icon,true);

    std::ofstream file_tr;
    file_tr.open("control_TR.txt");
    for ( unsigned i = 0; i < (unsigned)dim; i++ ) {
      file_tr << (*y_rcp)[i] << "\n";
    }
    file_tr.close();
   
    Teuchos::RCP<ROL::Vector<RealT> > error = x.clone();
    error->set(x);
    error->axpy(-1.0,y);
    std::cout << "\nError between PDAS solution and TR solution is " << error->norm() << "\n";
        
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

