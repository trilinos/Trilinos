// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BUNDLE_STEP_H
#define ROL_BUNDLE_STEP_H

#include "ROL_Bundle_AS.hpp"
#include "ROL_Bundle_TT.hpp"
#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_LineSearch.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"

/** @ingroup step_group
    \class ROL::BundleStep
    \brief Provides the interface to compute bundle trust-region steps.
*/

namespace ROL {

template <class Real>
class BundleStep : public Step<Real> {
private:
  // Bundle
  ROL::Ptr<Bundle<Real> >     bundle_;     // Bundle of subgradients and linearization errors
  ROL::Ptr<LineSearch<Real> > lineSearch_; // Line-search object for nonconvex problems

  // Dual cutting plane solution
  unsigned QPiter_;                        // Number of QP solver iterations
  unsigned QPmaxit_;                       // Maximum number of QP iterations
  Real QPtol_;                             // QP subproblem tolerance

  // Step flag
  int step_flag_;                          // Whether serious or null step

  // Additional storage
  ROL::Ptr<Vector<Real> > y_;

  // Updated iterate storage
  Real linErrNew_;
  Real valueNew_;

  // Aggregate subgradients, linearizations, and distance measures
  ROL::Ptr<Vector<Real> > aggSubGradNew_;  // New aggregate subgradient
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

  int verbosity_;

public:

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

  BundleStep(ROL::ParameterList &parlist)
    : bundle_(ROL::nullPtr), lineSearch_(ROL::nullPtr),
      QPiter_(0), QPmaxit_(0), QPtol_(0), step_flag_(0),
      y_(ROL::nullPtr), linErrNew_(0), valueNew_(0),
      aggSubGradNew_(ROL::nullPtr), aggSubGradOldNorm_(0),
      aggLinErrNew_(0), aggLinErrOld_(0), aggDistMeasNew_(0),
      T_(0), tol_(0), m1_(0), m2_(0), m3_(0), nu_(0),
      ls_maxit_(0), first_print_(true), isConvex_(false),
      ftol_(ROL_EPSILON<Real>()) {
    Real zero(0), two(2), oem3(1.e-3), oem6(1.e-6), oem8(1.e-8), p1(0.1), p2(0.2), p9(0.9), oe3(1.e3), oe8(1.e8);
    ROL::Ptr<StepState<Real> > state = Step<Real>::getState();
    state->searchSize = parlist.sublist("Step").sublist("Bundle").get("Initial Trust-Region Parameter", oe3);
    T_   = parlist.sublist("Step").sublist("Bundle").get("Maximum Trust-Region Parameter",       oe8); 
    tol_ = parlist.sublist("Step").sublist("Bundle").get("Epsilon Solution Tolerance",           oem6);
    m1_  = parlist.sublist("Step").sublist("Bundle").get("Upper Threshold for Serious Step",     p1);
    m2_  = parlist.sublist("Step").sublist("Bundle").get("Lower Threshold for Serious Step",     p2);
    m3_  = parlist.sublist("Step").sublist("Bundle").get("Upper Threshold for Null Step",        p9);
    nu_  = parlist.sublist("Step").sublist("Bundle").get("Tolerance for Trust-Region Parameter", oem3);

    // Initialize bundle
    Real coeff        = parlist.sublist("Step").sublist("Bundle").get("Distance Measure Coefficient",   zero);
    Real omega        = parlist.sublist("Step").sublist("Bundle").get("Locality Measure Coefficient",   two);
    unsigned maxSize  = parlist.sublist("Step").sublist("Bundle").get("Maximum Bundle Size",            200);
    unsigned remSize  = parlist.sublist("Step").sublist("Bundle").get("Removal Size for Bundle Update", 2);
    if ( parlist.sublist("Step").sublist("Bundle").get("Cutting Plane Solver",0) == 1 ) {
      bundle_ = ROL::makePtr<Bundle_TT<Real>>(maxSize,coeff,omega,remSize);
      //bundle_ = ROL::makePtr<Bundle_AS<Real>>(maxSize,coeff,omega,remSize);
    }
    else {
      bundle_ = ROL::makePtr<Bundle_AS<Real>>(maxSize,coeff,omega,remSize);
    }
    isConvex_ = ((coeff == zero) ? true : false);

    // Initialize QP solver 
    QPtol_   = parlist.sublist("Step").sublist("Bundle").get("Cutting Plane Tolerance", oem8);
    QPmaxit_ = parlist.sublist("Step").sublist("Bundle").get("Cutting Plane Iteration Limit", 1000);

    // Initialize Line Search
    ls_maxit_
      = parlist.sublist("Step").sublist("Line Search").get("Maximum Number of Function Evaluations",20);
    if ( !isConvex_ ) {
      lineSearch_ = LineSearchFactory<Real>(parlist);
    }

    // Get verbosity level
    verbosity_ = parlist.sublist("General").get("Print Verbosity", 0);
  }

  void initialize( Vector<Real> &x, const Vector<Real> &g, 
                   Objective<Real> &obj, BoundConstraint<Real> &con, 
                   AlgorithmState<Real> &algo_state ) { 
    // Call default initializer, but maintain current searchSize
    ROL::Ptr<StepState<Real> > state = Step<Real>::getState();
    Real searchSize = state->searchSize;
    Step<Real>::initialize(x,x,g,obj,con,algo_state);
    state->searchSize = searchSize;
    // Initialize bundle
    bundle_->initialize(*(state->gradientVec));
    // Initialize storage for updated iterate
    y_ = x.clone();
    // Initialize storage for aggregate subgradients
    aggSubGradNew_     = g.clone();
    aggSubGradOldNorm_ = algo_state.gnorm;
    // Initialize line search
    if ( !isConvex_ ) {
      lineSearch_->initialize(x,x,g,obj,con);
    }
  }

  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, 
                BoundConstraint<Real> &con, AlgorithmState<Real> &algo_state ) { 
    ROL::Ptr<StepState<Real> > state = Step<Real>::getState();
    first_print_ = false;                     // Print header only on first serious step
    QPiter_ = (step_flag_==1 ? 0 : QPiter_);  // Reset QPiter only on serious steps
    Real v(0), l(0), u = T_, gd(0);           // Scalar storage
    Real zero(0), two(2), half(0.5);
    bool flag = true;
    while (flag) {
      /*************************************************************/
      /******** Solve Dual Cutting Plane QP Problem ****************/
      /*************************************************************/
      QPiter_ += bundle_->solveDual(state->searchSize,QPmaxit_,QPtol_);  // Solve QP subproblem
      bundle_->aggregate(*aggSubGradNew_,aggLinErrNew_,aggDistMeasNew_); // Compute aggregate info
      algo_state.aggregateGradientNorm = aggSubGradNew_->norm();         // Aggregate subgradient norm
      if (verbosity_ > 0) {
        std::cout << std::endl;
        std::cout << "  Computation of aggregrate quantities" << std::endl;
        std::cout << "    Aggregate subgradient norm:       " << algo_state.aggregateGradientNorm << std::endl;
        std::cout << "    Aggregate linearization error:    " << aggLinErrNew_ << std::endl;
        std::cout << "    Aggregate distance measure:       " << aggDistMeasNew_ << std::endl;
      }
      /*************************************************************/
      /******** Construct Cutting Plane Solution *******************/
      /*************************************************************/
      v = -state->searchSize*std::pow(algo_state.aggregateGradientNorm,two)-aggLinErrNew_; // CP objective value
      s.set(aggSubGradNew_->dual()); s.scale(-state->searchSize);            // CP solution
      algo_state.snorm = state->searchSize*algo_state.aggregateGradientNorm; // Step norm
      if (verbosity_ > 0) {
        std::cout << std::endl;
        std::cout << "  Solve cutting plan subproblem" << std::endl;
        std::cout << "    Cutting plan objective value:     " << v << std::endl;
        std::cout << "    Norm of computed step:            " << algo_state.snorm << std::endl;
        std::cout << "    'Trust-region' radius:            " << state->searchSize << std::endl;
      }
      /*************************************************************/
      /******** Decide Whether Step is Serious or Null *************/
      /*************************************************************/
      if (std::max(algo_state.aggregateGradientNorm,aggLinErrNew_) <= tol_) {
        // Current iterate is already epsilon optimal!
        s.zero(); algo_state.snorm = zero;
        flag = false;
        step_flag_ = 1;
        algo_state.flag = true;
        break;
      }
      else if (std::isnan(algo_state.aggregateGradientNorm)
               || std::isnan(aggLinErrNew_)
               || (std::isnan(aggDistMeasNew_) && !isConvex_)) {
        s.zero(); algo_state.snorm = zero;
        flag = false;
        step_flag_ = 2;
        algo_state.flag = true;
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
        gd = s.dot(state->gradientVec->dual());
        linErrNew_ = algo_state.value - (valueNew_ - gd); // Linearization error
        // Determine whether to take a serious or null step
        Real eps  = static_cast<Real>(10)*ROL_EPSILON<Real>();
        Real del  = eps*std::max(static_cast<Real>(1),std::abs(algo_state.value));
        Real Df   = (valueNew_ - algo_state.value) - del;
        Real Dm   = v - del;
        bool SS1  = false;
        if (std::abs(Df) < eps && std::abs(Dm) < eps) {
          SS1 = true;
        }
        else {
          SS1 = (Df < m1_*Dm);
        }
        //bool SS1  = (valueNew_-algo_state.value <  m1_*v);
        //bool NS1  = (valueNew_-algo_state.value >= m1_*v);
        bool NS2a = (bundle_->computeAlpha(algo_state.snorm,linErrNew_) <= m3_*aggLinErrOld_);
        bool NS2b = (std::abs(algo_state.value-valueNew_) <= aggSubGradOldNorm_ + aggLinErrOld_);
        if (verbosity_ > 0) {
          std::cout << std::endl;
          std::cout << "  Check for serious/null step" << std::endl;
          std::cout << "    Serious step test SS(i):          " << SS1 << std::endl;
          std::cout << "       -> Left hand side:             " << valueNew_-algo_state.value << std::endl;
          std::cout << "       -> Right hand side:            " << m1_*v << std::endl;
          std::cout << "    Null step test NS(iia):           " << NS2a << std::endl;
          std::cout << "       -> Left hand side:             " << bundle_->computeAlpha(algo_state.snorm,linErrNew_) << std::endl;
          std::cout << "       -> Right hand side:            " << m3_*aggLinErrOld_ << std::endl;
          std::cout << "    Null step test NS(iib):           " << NS2b << std::endl;
          std::cout << "       -> Left hand side:             " << std::abs(algo_state.value-valueNew_) << std::endl;
          std::cout << "       -> Right hand side:            " << aggSubGradOldNorm_ + aggLinErrOld_ << std::endl;
        }
        if ( isConvex_ ) {
          /************* Convex objective ****************/
          if ( SS1 ) {
            bool SS2 = (gd >= m2_*v || state->searchSize >= T_-nu_);
            if (verbosity_ > 0) {
              std::cout << "    Serious step test SS(iia):        " << (gd >= m2_*v) << std::endl;
              std::cout << "       -> Left hand side:             " << gd << std::endl;
              std::cout << "       -> Right hand side:            " << m2_*v << std::endl;
              std::cout << "    Serious step test SS(iia):        " << (state->searchSize >= T_-nu_) << std::endl;
              std::cout << "       -> Left hand side:             " << state->searchSize << std::endl;
              std::cout << "       -> Right hand side:            " << T_-nu_ << std::endl;
            }
            if ( SS2 ) { // Serious Step
              step_flag_ = 1;
              flag       = false;
              if (verbosity_ > 0) {
                std::cout << "  Serious step taken" << std::endl;
              }
              break;
            }
            else { // Increase trust-region radius
              l = state->searchSize;
              state->searchSize = half*(u+l);
              if (verbosity_ > 0) {
                std::cout << "    Increase 'trust-region' radius:   " << state->searchSize << std::endl;
              }
            }
          }
          else {
            if ( NS2a || NS2b ) { // Null step
              s.zero(); algo_state.snorm = zero;
              step_flag_ = 0;
              flag       = false;
              if (verbosity_ > 0) {
                std::cout << "  Null step taken" << std::endl;
              }
              break;
            }
            else { // Decrease trust-region radius
              u = state->searchSize;
              state->searchSize = half*(u+l);
              if (verbosity_ > 0) {
                std::cout << "    Decrease 'trust-region' radius:   " << state->searchSize << std::endl;
              }
            }
          }
        }
        else {
          /************* Nonconvex objective *************/
          bool NS3 = (gd - bundle_->computeAlpha(algo_state.snorm,linErrNew_) >= m2_*v);
          if (verbosity_ > 0) {
            std::cout << "    Null step test NS(iii):           " << NS3 << std::endl;
            std::cout << "       -> Left hand side:             " << gd - bundle_->computeAlpha(algo_state.snorm,linErrNew_) << std::endl;
            std::cout << "       -> Right hand side:            " << m2_*v << std::endl;
          }
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
                  Real alpha = zero;
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
                  state->searchSize = half*(u+l);
                }
              }
            }
            else { // Decrease trust-region radius
              u = state->searchSize;
              state->searchSize = half*(u+l);
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
    ROL::Ptr<StepState<Real> > state = Step<Real>::getState();
    state->flag = step_flag_;
    state->SPiter = QPiter_;
    if ( !algo_state.flag ) {
      /*************************************************************/
      /******** Reset Bundle If Maximum Size Reached ***************/
      /*************************************************************/
      bundle_->reset(*aggSubGradNew_,aggLinErrNew_,algo_state.snorm);
      /*************************************************************/
      /******** Update Bundle with Step Information ****************/
      /*************************************************************/
      if ( step_flag_==1 ) {
        // Serious step was taken
        x.plus(s);                        // Update iterate
        Real valueOld = algo_state.value; // Store previous objective value
        algo_state.value = valueNew_;     // Add new objective value to state
        bundle_->update(step_flag_,valueNew_-valueOld,algo_state.snorm,*(state->gradientVec),s);
      }
      else if ( step_flag_==0 ) {
        // Null step was taken
        bundle_->update(step_flag_,linErrNew_,algo_state.snorm,*(state->gradientVec),s);
      }
    }
    /*************************************************************/
    /******** Update Algorithm State *****************************/
    /*************************************************************/
    algo_state.iterateVec->set(x);
    algo_state.gnorm = (state->gradientVec)->norm();
    if ( step_flag_==1 ) {
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
    const ROL::Ptr<const StepState<Real> > state = Step<Real>::getStepState();
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
    if ( step_flag_==1 && algo_state.iter > 0 ) {
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
