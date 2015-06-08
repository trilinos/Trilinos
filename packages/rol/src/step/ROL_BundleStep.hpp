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

#ifndef ROL_BUNDLE_STEP_H
#define ROL_BUNDLE_STEP_H

#include "ROL_Bundle.hpp"
#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_LineSearch.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

/** @ingroup step_group
    \class ROL::BundleStep
    \brief Provides the interface to compute bundle trust-region steps.
*/

namespace ROL {

template <class Real>
class BundleStep : public Step<Real> {
private:
  // Bundle
  Teuchos::RCP<Bundle<Real> >     bundle_;     // Bundle of subgradients and linearization errors
  Teuchos::RCP<LineSearch<Real> > lineSearch_; // Line-search object for nonconvex problems

  // Dual cutting plane solution
  unsigned QPiter_;                        // Number of QP solver iterations
  unsigned QPmaxit_;                       // Maximum number of QP iterations
  Real QPtol_;                             // QP subproblem tolerance

  // Step flag
  int step_flag_;                          // Whether serious or null step

  // Additional storage
  Teuchos::RCP<Vector<Real> > y_;

  // Updated iterate storage
  Real linErrNew_;
  Real valueNew_;

  // Aggregate subgradients, linearizations, and distance measures
  Teuchos::RCP<Vector<Real> > aggSubGradNew_;  // New aggregate subgradient
  Real aggSubGradOldNorm_;                          // Old aggregate subgradient norm
  Real aggLinErrNew_;                               // New aggregate linearization error
  Real aggLinErrOld_;                               // Old aggregate linearization error
  Real aggDistMeasNew_;                             // New aggregate distance measure

  // Algorithmic parameters
  Real T_;
  Real tol_;
  Real m1_;
  Real m2_;
  Real m3_;
  Real nu_;

  // Line-search parameters
  int ls_maxit_;

  bool first_print_;
  bool isConvex_;

  Real ftol_;

  void LineSearchFactory(Teuchos::ParameterList &parlist) {
    ls_maxit_ = parlist.get("Maximum Number of Function Evaluations",20);
    ELineSearch els = StringToELineSearch(parlist.get("Linesearch Type","Cubic Interpolation"));
    switch(els) {
      case LINESEARCH_ITERATIONSCALING:
        lineSearch_ = Teuchos::rcp( new IterationScaling<Real>(parlist) );     break;
      case LINESEARCH_PATHBASEDTARGETLEVEL:
        lineSearch_ = Teuchos::rcp( new PathBasedTargetLevel<Real>(parlist) ); break;
      case LINESEARCH_BACKTRACKING:
        lineSearch_ = Teuchos::rcp( new BackTracking<Real>(parlist) );         break;
      case LINESEARCH_BISECTION:
        lineSearch_ = Teuchos::rcp( new Bisection<Real>(parlist) );            break;
      case LINESEARCH_BRENTS:
        lineSearch_ = Teuchos::rcp( new Brents<Real>(parlist) );               break;
      case LINESEARCH_GOLDENSECTION:
        lineSearch_ = Teuchos::rcp( new GoldenSection<Real>(parlist) );        break;
      case LINESEARCH_CUBICINTERP:
      default:
        lineSearch_ = Teuchos::rcp( new CubicInterp<Real>(parlist) );          break;
    }
  }

public:

  BundleStep(Teuchos::ParameterList &parlist)
    : QPiter_(0), step_flag_(0), aggLinErrNew_(0.0), aggLinErrOld_(0.0), aggDistMeasNew_(0.0),
      first_print_(true), ftol_(ROL_EPSILON) {
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    state->searchSize = parlist.get("Bundle Step: Initial Trust-Region Parameter",1.e3);
    T_                = parlist.get("Bundle Step: Maximum Trust-Region Parameter",1.e8); 
    tol_              = parlist.get("Bundle Step: Epsilon Solution Tolerance",1.e-6);
    m1_               = parlist.get("Bundle Step: Upper Threshold for Serious Step",0.1);
    m2_               = parlist.get("Bundle Step: Lower Threshold for Serious Step",0.2);
    m3_               = parlist.get("Bundle Step: Upper Threshold for Null Step",0.9);
    nu_               = parlist.get("Bundle Step: Tolerance for Trust-Region Parameter",1.e-3);

    // Initialize bundle
    Real coeff        = parlist.get("Bundle Step: Distance Measure Coefficient",0.0);
    unsigned maxSize  = parlist.get("Bundle Step: Maximum Bundle Size",200);
    unsigned remSize  = parlist.get("Bundle Step: Removal Size for Bundle Update",2);
    bundle_   = Teuchos::rcp(new Bundle<Real>(maxSize,coeff,remSize));
    isConvex_ = ((coeff == 0.0) ? true : false);

    // Initialize QP solver 
    QPtol_   = parlist.get("Bundle Step: Cutting Plane Solver Tolerance",1.e-8);
    QPmaxit_ = parlist.get("Bundle Step: Cutting Plane Solver Maximum Number of Iterations",1000);

    // Initialize Line Search
    lineSearch_ = Teuchos::null;
    if ( !isConvex_ ) {
      LineSearchFactory(parlist);
    }
  }

  void initialize( Vector<Real> &x, const Vector<Real> &g, 
                   Objective<Real> &obj, BoundConstraint<Real> &con, 
                   AlgorithmState<Real> &algo_state ) { 
    // Call default initializer, but maintain current searchSize
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    Real searchSize = state->searchSize;
    Step<Real>::initialize(x,x,g,obj,con,algo_state);
    state->searchSize = searchSize;
    // Initialize bundle
    bundle_->initialize(*(state->gradientVec));
    // Initialize storage for updated iterate
    y_ = x.clone();
    // Initialize storage for aggregate subgradients
    aggSubGradNew_ = x.clone();
    aggSubGradOldNorm_ = algo_state.gnorm;
  }

  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, 
                BoundConstraint<Real> &con, AlgorithmState<Real> &algo_state ) { 
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    first_print_ = false;                     // Print header only on first serious step
    QPiter_ = (step_flag_ ? 0 : QPiter_);     // Reset QPiter only on serious steps
    Real v = 0.0, l = 0.0, u = T_, gd = 0.0;  // Scalar storage
    bool flag = true;
    while (flag) {
      /*************************************************************/
      /******** Solve Dual Cutting Plane QP Problem ****************/
      /*************************************************************/
      QPiter_ += bundle_->solveDual(state->searchSize,QPmaxit_,QPtol_);  // Solve QP subproblem
      bundle_->aggregate(*aggSubGradNew_,aggLinErrNew_,aggDistMeasNew_); // Compute aggregate info
      algo_state.aggregateGradientNorm = aggSubGradNew_->norm();         // Aggregate subgradient norm
      /*************************************************************/
      /******** Construct Cutting Plane Solution *******************/
      /*************************************************************/
      v = -state->searchSize*std::pow(algo_state.aggregateGradientNorm,2.0)-aggLinErrNew_; // CP objective value
      s.set(*aggSubGradNew_); s.scale(-state->searchSize);                        // CP solution
      algo_state.snorm = state->searchSize*algo_state.aggregateGradientNorm;      // Step norm
      /*************************************************************/
      /******** Decide Whether Step is Serious or Null *************/
      /*************************************************************/
      if (std::max(algo_state.aggregateGradientNorm,aggLinErrNew_) <= tol_) {
        // Current iterate is already epsilon optimal!
        s.zero(); algo_state.snorm = 0.0;
        flag = false;
        step_flag_ = 1;
        algo_state.flag = true;
        break;
      }
      else {
        // Current iterate is not epsilon optimal.
        y_->set(x); y_->plus(s);                       // y is the candidate iterate
        obj.update(*y_,true,algo_state.iter);          // Update objective at y
        valueNew_ = obj.value(*y_,ftol_);              // Compute objective value at y
        algo_state.nfval++;
        obj.gradient(*(state->gradientVec),*y_,ftol_); // Compute objective (sub)gradient at y
        algo_state.ngrad++;
        // Compute new linearization error and distance measure
        gd = s.dot(*(state->gradientVec));
        linErrNew_ = algo_state.value - (valueNew_ - gd); // Linearization error
        // Determine whether to take a serious or null step
        bool SS1  = (valueNew_-algo_state.value <  m1_*v);
        //bool NS1  = (valueNew_-algo_state.value >= m1_*v);
        bool NS2a = (bundle_->computeAlpha(algo_state.snorm,linErrNew_) <= m3_*aggLinErrOld_);
        bool NS2b = (std::abs(algo_state.value-valueNew_) <= aggSubGradOldNorm_ + aggLinErrOld_);
        if ( isConvex_ ) {
          /************* Convex objective ****************/
          if ( SS1 ) {
            bool SS2 = (gd >= m2_*v || state->searchSize >= T_-nu_);
            if ( SS2 ) { // Serious Step
              step_flag_ = 1;
              flag       = false;
              break;
            }
            else { // Increase trust-region radius
              l = state->searchSize;
              state->searchSize = 0.5*(u+l);
            }
          }
          else {
            if ( NS2a || NS2b ) { // Null step
              s.zero(); algo_state.snorm = 0.0;
              step_flag_ = 0;
              flag       = false;
              break;
            }
            else { // Decrease trust-region radius
              u = state->searchSize;
              state->searchSize = 0.5*(u+l);
            }
          }
        }
        else {
          /************* Nonconvex objective *************/
          bool NS3 = (gd - bundle_->computeAlpha(algo_state.snorm,linErrNew_) >= m2_*v);
          if ( SS1 ) { // Serious step
            step_flag_ = 1;
            flag       = false;
            break;
          }
          else {
            if ( NS2a || NS2b ) {
              if ( NS3 ) { // Null step
                s.zero();
                step_flag_ = 0;
                flag       = false;
                break;
              }
              else {
                if ( NS2b ) { // Line search
                  Real alpha = 0.0;
                  int ls_nfval = 0, ls_ngrad = 0;
                  lineSearch_->run(alpha,valueNew_,ls_nfval,ls_ngrad,gd,s,x,obj,con);
                  if ( ls_nfval == ls_maxit_ ) { // Null step
                    s.zero();
                    step_flag_ = 0;
                    flag       = false;
                    break;
                  }
                  else { // Serious step
                    s.scale(alpha);
                    step_flag_ = 1;
                    flag       = false;
                    break;
                  }
                }
                else { // Decrease trust-region radius
                  u = state->searchSize;
                  state->searchSize = 0.5*(u+l);
                }
              }
            }
            else { // Decrease trust-region radius
              u = state->searchSize;
              state->searchSize = 0.5*(u+l);
            }
          }
        }
      }
    } // End While
    /*************************************************************/
    /******** Update Algorithm State *****************************/
    /*************************************************************/
    algo_state.aggregateModelError = aggLinErrNew_;
    aggSubGradOldNorm_ = algo_state.aggregateGradientNorm;
    aggLinErrOld_      = aggLinErrNew_;
  } // End Compute

  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, 
               BoundConstraint<Real> &con, AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    if ( !algo_state.flag ) {
      /*************************************************************/
      /******** Reset Bundle If Maximum Size Reached ***************/
      /*************************************************************/
      bundle_->reset(*aggSubGradNew_,aggLinErrNew_,algo_state.snorm);
      /*************************************************************/
      /******** Update Bundle with Step Information ****************/
      /*************************************************************/
      if ( step_flag_ ) {
        // Serious step was taken
        x.plus(s);                        // Update iterate
        Real valueOld = algo_state.value; // Store previous objective value
        algo_state.value = valueNew_;     // Add new objective value to state
        bundle_->update(step_flag_,valueNew_-valueOld,algo_state.snorm,*(state->gradientVec),s);
      }
      else {
        // Null step was taken
        bundle_->update(step_flag_,linErrNew_,algo_state.snorm,*(state->gradientVec),s);
      }
    }
    /*************************************************************/
    /******** Update Algorithm State *****************************/
    /*************************************************************/
    algo_state.iterateVec->set(x);
    algo_state.gnorm = (state->gradientVec)->norm();
    if ( step_flag_ ) {
      algo_state.iter++;
    }
  } // End Update

  std::string printHeader( void ) const {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6) << std::left << "iter";
    hist << std::setw(15) << std::left << "value";
    hist << std::setw(15) << std::left << "gnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    hist << std::setw(15) << std::left << "znorm";
    hist << std::setw(15) << std::left << "alpha";
    hist << std::setw(15) << std::left << "TRparam";
    hist << std::setw(10) << std::left << "QPiter";
    hist << "\n";
    return hist.str();
  }

  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << "Bundle Trust-Region Algorithm \n";
    return hist.str();
  }

  std::string print( AlgorithmState<Real> &algo_state, bool print_header = false ) const {
    Teuchos::RCP<const StepState<Real> > state = Step<Real>::getStepState();
    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 && first_print_ ) {
      hist << printName();
      if ( print_header ) {
        hist << printHeader();
      }
      hist << "  ";
      hist << std::setw(6) << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << "\n";
    }
    if ( step_flag_ && algo_state.iter > 0 ) {
      if ( print_header ) {
        hist << printHeader();
      }
      else {
        hist << "  ";
        hist << std::setw(6) << std::left << algo_state.iter;
        hist << std::setw(15) << std::left << algo_state.value;
        hist << std::setw(15) << std::left << algo_state.gnorm;
        hist << std::setw(15) << std::left << algo_state.snorm;
        hist << std::setw(10) << std::left << algo_state.nfval;
        hist << std::setw(10) << std::left << algo_state.ngrad;
        hist << std::setw(15) << std::left << algo_state.aggregateGradientNorm;
        hist << std::setw(15) << std::left << algo_state.aggregateModelError;
        hist << std::setw(15) << std::left << state->searchSize;
        hist << std::setw(10) << std::left << QPiter_;
        hist << "\n";
      }
    }
    return hist.str();
  };

}; // class BundleStep

} // namespace ROL

#endif
