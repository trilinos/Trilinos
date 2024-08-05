// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEU_BUNDLEALGORITHM_DEF_H
#define ROL_TYPEU_BUNDLEALGORITHM_DEF_H

#include "ROL_BundleStatusTest.hpp"
#include "ROL_Bundle_U_AS.hpp"
#include "ROL_Bundle_U_TT.hpp"
#include "ROL_LineSearch_U_Factory.hpp"

namespace ROL {
namespace TypeU {


template<typename Real>
BundleAlgorithm<Real>::BundleAlgorithm( ParameterList &parlist,
                                        const Ptr<LineSearch_U<Real>> &lineSearch )
  : Algorithm<Real>(),
    bundle_(ROL::nullPtr), lineSearch_(ROL::nullPtr),
    QPiter_(0), QPmaxit_(0), QPtol_(0), step_flag_(0),
    T_(0), tol_(0), m1_(0), m2_(0), m3_(0), nu_(0),
    ls_maxit_(0), first_print_(true), isConvex_(false) {
  // Set bundle status test
  status_->reset();
  status_->add(makePtr<BundleStatusTest<Real>>(parlist));

  // Parse parameter parlist
  const Real zero(0), two(2), oem3(1.e-3), oem6(1.e-6), oem8(1.e-8);
  const Real p1(0.1), p2(0.2), p9(0.9), oe3(1.e3), oe8(1.e8);
  ParameterList &blist = parlist.sublist("Step").sublist("Bundle");
  state_->searchSize = blist.get("Initial Trust-Region Parameter", oe3);
  T_   = blist.get("Maximum Trust-Region Parameter",       oe8); 
  tol_ = blist.get("Epsilon Solution Tolerance",           oem6);
  m1_  = blist.get("Upper Threshold for Serious Step",     p1);
  m2_  = blist.get("Lower Threshold for Serious Step",     p2);
  m3_  = blist.get("Upper Threshold for Null Step",        p9);
  nu_  = blist.get("Tolerance for Trust-Region Parameter", oem3);

  // Initialize bundle
  Real coeff        = blist.get("Distance Measure Coefficient",   zero);
  Real omega        = blist.get("Locality Measure Coefficient",   two);
  unsigned maxSize  = blist.get("Maximum Bundle Size",            200);
  unsigned remSize  = blist.get("Removal Size for Bundle Update", 2);
  if ( blist.get("Cutting Plane Solver",0) == 1 ) {
    bundle_ = makePtr<Bundle_U_TT<Real>>(maxSize,coeff,omega,remSize);
  }
  else {
    bundle_ = makePtr<Bundle_U_AS<Real>>(maxSize,coeff,omega,remSize);
  }
  isConvex_ = ((coeff == zero) ? true : false);

  // Initialize QP solver 
  QPtol_   = blist.get("Cutting Plane Tolerance", oem8);
  QPmaxit_ = blist.get("Cutting Plane Iteration Limit", 1000);

  // Initialize Line Search
  ParameterList &lslist = parlist.sublist("Step").sublist("Line Search");
  ls_maxit_ = lslist.get("Maximum Number of Function Evaluations",20);
  if ( !isConvex_ && lineSearch_==nullPtr ) {
    lineSearch_ = LineSearchUFactory<Real>(parlist);
  }

  // Get verbosity level
  verbosity_ = parlist.sublist("General").get("Output Level", 0);
  printHeader_ = verbosity_ > 2;
}

template<typename Real>
void BundleAlgorithm<Real>::initialize( const Vector<Real> &x,
                                        const Vector<Real> &g,
                                        Objective<Real>    &obj,
                                        std::ostream &outStream) {
  // Initialize data
  Algorithm<Real>::initialize(x,g);
  if (!isConvex_) {
    lineSearch_->initialize(x,g);
  }
  // Update objective function, get value and gradient
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  obj.update(x,UpdateType::Initial,state_->iter);
  state_->value = obj.value(x,tol);
  state_->nfval++;
  obj.gradient(*state_->gradientVec,x,tol);
  state_->ngrad++;
  state_->gnorm = state_->gradientVec->norm();
  bundle_->initialize(*state_->gradientVec);
}

template<typename Real>
void BundleAlgorithm<Real>::run( Vector<Real>       &x,
                                 const Vector<Real> &g,
                                 Objective<Real>    &obj,
                                 std::ostream       &outStream ) {
  const Real zero(0), two(2), half(0.5);
  // Initialize trust-region data
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  initialize(x,g,obj,outStream);
  Ptr<Vector<Real>> y = x.clone();
  Ptr<Vector<Real>> aggSubGradNew = g.clone();
  Real aggSubGradOldNorm = state_->gnorm;
  Real linErrNew(0), valueNew(0);
  Real aggLinErrNew(0), aggLinErrOld(0), aggDistMeasNew(0);
  Real v(0), l(0), u(T_), gd(0);
  bool flag = true;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    first_print_ = false; // Print header only on first serious step
    QPiter_ = (step_flag_==1 ? 0 : QPiter_); // Reset QPiter only on serious steps
    v = zero; l = zero; u = T_; gd = zero;
    flag = true;
    while (flag) {
      /*************************************************************/
      /******** Solve Dual Cutting Plane QP Problem ****************/
      /*************************************************************/
      QPiter_ += bundle_->solveDual(state_->searchSize,QPmaxit_,QPtol_); // Solve QP subproblem
      bundle_->aggregate(*aggSubGradNew,aggLinErrNew,aggDistMeasNew); // Compute aggregate info
      state_->aggregateGradientNorm = aggSubGradNew->norm();            // Aggregate subgradient norm
      if (verbosity_ > 1) {
        outStream << std::endl;
        outStream << "  Computation of aggregrate quantities" << std::endl;
        outStream << "    Aggregate subgradient norm:       " << state_->aggregateGradientNorm << std::endl;
        outStream << "    Aggregate linearization error:    " << aggLinErrNew << std::endl;
        outStream << "    Aggregate distance measure:       " << aggDistMeasNew << std::endl;
      }
      /*************************************************************/
      /******** Construct Cutting Plane Solution *******************/
      /*************************************************************/
      v = -state_->searchSize*std::pow(state_->aggregateGradientNorm,two)-aggLinErrNew; // CP objective value
      state_->stepVec->set(aggSubGradNew->dual());
      state_->stepVec->scale(-state_->searchSize); // CP solution
      state_->snorm = state_->searchSize*state_->aggregateGradientNorm; // Step norm
      if (verbosity_ > 1) {
        outStream << std::endl;
        outStream << "  Solve cutting plan subproblem" << std::endl;
        outStream << "    Cutting plan objective value:     " << v << std::endl;
        outStream << "    Norm of computed step:            " << state_->snorm << std::endl;
        outStream << "    'Trust-region' radius:            " << state_->searchSize << std::endl;
      }
      /*************************************************************/
      /******** Decide Whether Step is Serious or Null *************/
      /*************************************************************/
      if (std::max(state_->aggregateGradientNorm,aggLinErrNew) <= tol_) {
        // Current iterate is already epsilon optimal!
        state_->stepVec->zero(); state_->snorm = zero;
        flag = false;
        step_flag_ = 1;
        state_->flag = true;
        break;
      }
      else if (std::isnan(state_->aggregateGradientNorm)
               || std::isnan(aggLinErrNew)
               || (std::isnan(aggDistMeasNew) && !isConvex_)) {
        state_->stepVec->zero(); state_->snorm = zero;
        flag = false;
        step_flag_ = 2;
        state_->flag = true;
      }
      else {
        // Current iterate is not epsilon optimal.
        y->set(x); y->plus(*state_->stepVec);      // y is the candidate iterate
        obj.update(*y,UpdateType::Accept,state_->iter); // Update objective at y
        valueNew = obj.value(*y,tol);              // Compute objective value at y
        state_->nfval++;
        obj.gradient(*state_->gradientVec,*y,tol); // Compute objective (sub)gradient at y
        state_->ngrad++;
        // Compute new linearization error and distance measure
        //gd = state_->stepVec->dot(state_->gradientVec->dual());
        gd = state_->stepVec->apply(*state_->gradientVec);
        linErrNew = state_->value - (valueNew - gd); // Linearization error
        // Determine whether to take a serious or null step
        Real eps  = static_cast<Real>(10)*ROL_EPSILON<Real>();
        Real del  = eps*std::max(static_cast<Real>(1),std::abs(state_->value));
        Real Df   = (valueNew - state_->value) - del;
        Real Dm   = v - del;
        bool SS1  = false;
        if (std::abs(Df) < eps && std::abs(Dm) < eps) {
          SS1 = true;
        }
        else {
          SS1 = (Df < m1_*Dm);
        }
        //bool SS1  = (valueNew-state_->value <  m1_*v);
        //bool NS1  = (valueNew-state_->value >= m1_*v);
        bool NS2a = (bundle_->computeAlpha(state_->snorm,linErrNew) <= m3_*aggLinErrOld);
        bool NS2b = (std::abs(state_->value-valueNew) <= aggSubGradOldNorm + aggLinErrOld);
        if (verbosity_ > 1) {
          outStream << std::endl;
          outStream << "  Check for serious/null step" << std::endl;
          outStream << "    Serious step test SS(i):          " << SS1 << std::endl;
          outStream << "       -> Left hand side:             " << valueNew-state_->value << std::endl;
          outStream << "       -> Right hand side:            " << m1_*v << std::endl;
          outStream << "    Null step test NS(iia):           " << NS2a << std::endl;
          outStream << "       -> Left hand side:             " << bundle_->computeAlpha(state_->snorm,linErrNew) << std::endl;
          outStream << "       -> Right hand side:            " << m3_*aggLinErrOld << std::endl;
          outStream << "    Null step test NS(iib):           " << NS2b << std::endl;
          outStream << "       -> Left hand side:             " << std::abs(state_->value-valueNew) << std::endl;
          outStream << "       -> Right hand side:            " << aggSubGradOldNorm + aggLinErrOld << std::endl;
        }
        if ( isConvex_ ) {
          /************* Convex objective ****************/
          if ( SS1 ) {
            bool SS2 = (gd >= m2_*v || state_->searchSize >= T_-nu_);
            if (verbosity_ > 1) {
              outStream << "    Serious step test SS(iia):        " << (gd >= m2_*v) << std::endl;
              outStream << "       -> Left hand side:             " << gd << std::endl;
              outStream << "       -> Right hand side:            " << m2_*v << std::endl;
              outStream << "    Serious step test SS(iia):        " << (state_->searchSize >= T_-nu_) << std::endl;
              outStream << "       -> Left hand side:             " << state_->searchSize << std::endl;
              outStream << "       -> Right hand side:            " << T_-nu_ << std::endl;
            }
            if ( SS2 ) { // Serious Step
              step_flag_ = 1;
              flag       = false;
              if (verbosity_ > 1) {
                outStream << "  Serious step taken" << std::endl;
              }
              break;
            }
            else { // Increase trust-region radius
              l = state_->searchSize;
              state_->searchSize = half*(u+l);
              if (verbosity_ > 1) {
                outStream << "    Increase 'trust-region' radius:   " << state_->searchSize << std::endl;
              }
            }
          }
          else {
            if ( NS2a || NS2b ) { // Null step
              state_->stepVec->zero(); state_->snorm = zero;
              step_flag_ = 0;
              flag       = false;
              if (verbosity_ > 1) {
                outStream << "  Null step taken" << std::endl;
              }
              break;
            }
            else { // Decrease trust-region radius
              u = state_->searchSize;
              state_->searchSize = half*(u+l);
              if (verbosity_ > 1) {
                outStream << "    Decrease 'trust-region' radius:   " << state_->searchSize << std::endl;
              }
            }
          }
        }
        else {
          /************* Nonconvex objective *************/
          bool NS3 = (gd - bundle_->computeAlpha(state_->snorm,linErrNew) >= m2_*v);
          if (verbosity_ > 1) {
            outStream << "    Null step test NS(iii):           " << NS3 << std::endl;
            outStream << "       -> Left hand side:             " << gd - bundle_->computeAlpha(state_->snorm,linErrNew) << std::endl;
            outStream << "       -> Right hand side:            " << m2_*v << std::endl;
          }
          if ( SS1 ) { // Serious step
            step_flag_ = 1;
            flag       = false;
            break;
          }
          else {
            if ( NS2a || NS2b ) {
              if ( NS3 ) { // Null step
                state_->stepVec->zero();
                step_flag_ = 0;
                flag       = false;
                break;
              }
              else {
                if ( NS2b ) { // Line search
                  Real alpha = zero;
                  int ls_nfval = 0, ls_ngrad = 0;
                  lineSearch_->run(alpha,valueNew,ls_nfval,ls_ngrad,gd,*state_->stepVec,x,obj);
                  if ( ls_nfval == ls_maxit_ ) { // Null step
                    state_->stepVec->zero();
                    step_flag_ = 0;
                    flag       = false;
                    break;
                  }
                  else { // Serious step
                    state_->stepVec->scale(alpha);
                    step_flag_ = 1;
                    flag       = false;
                    break;
                  }
                }
                else { // Decrease trust-region radius
                  u = state_->searchSize;
                  state_->searchSize = half*(u+l);
                }
              }
            }
            else { // Decrease trust-region radius
              u = state_->searchSize;
              state_->searchSize = half*(u+l);
            }
          }
        }
      }
    } // End While
    /*************************************************************/
    /******** Update Algorithm State *****************************/
    /*************************************************************/
    state_->aggregateModelError = aggLinErrNew;
    aggSubGradOldNorm = state_->aggregateGradientNorm;
    aggLinErrOld      = aggLinErrNew;

    if ( !state_->flag ) {
      /*************************************************************/
      /******** Reset Bundle If Maximum Size Reached ***************/
      /*************************************************************/
      bundle_->reset(*aggSubGradNew,aggLinErrNew,state_->snorm);
      /*************************************************************/
      /******** Update Bundle with Step Information ****************/
      /*************************************************************/
      if ( step_flag_==1 ) {
        // Serious step was taken
        x.plus(*state_->stepVec);      // Update iterate
        Real valueOld = state_->value; // Store previous objective value
        state_->value = valueNew;      // Add new objective value to state
        bundle_->update(step_flag_,valueNew-valueOld,state_->snorm,*state_->gradientVec,*state_->stepVec);
      }
      else if ( step_flag_==0 ) {
        // Null step was taken
        bundle_->update(step_flag_,linErrNew,state_->snorm,*state_->gradientVec,*state_->stepVec);
      }
    }
    /*************************************************************/
    /******** Update Algorithm State *****************************/
    /*************************************************************/
    state_->iterateVec->set(x);
    state_->gnorm = state_->gradientVec->norm();
    if ( step_flag_==1 ) {
      state_->iter++;
    }

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,printHeader_);
  }
  if (verbosity_ > 0) Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void BundleAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << "  ";
  os << std::setw(6) << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "gnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(10) << std::left << "#fval";
  os << std::setw(10) << std::left << "#grad";
  os << std::setw(15) << std::left << "znorm";
  os << std::setw(15) << std::left << "alpha";
  os << std::setw(15) << std::left << "TRparam";
  os << std::setw(10) << std::left << "QPiter";
  os << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void BundleAlgorithm<Real>::writeName(std::ostream& os) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::endl << "Bundle Trust-Region Algorithm" << std::endl;
  os.flags(osFlags);
}

template<typename Real>
void BundleAlgorithm<Real>::writeOutput( std::ostream& os, bool print_header) const {
  std::ios_base::fmtflags osFlags(os.flags());
  os << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 && first_print_ ) {
    writeName(os);
    if ( print_header ) {
      writeHeader(os);
    }
    os << "  ";
    os << std::setw(6) << std::left << state_->iter;
    os << std::setw(15) << std::left << state_->value;
    os << std::setw(15) << std::left << state_->gnorm;
    os << std::setw(15) << std::left << "---";
    os << std::setw(10) << std::left << state_->nfval;
    os << std::setw(10) << std::left << state_->ngrad;
    os << std::setw(15) << std::left << "---";
    os << std::setw(15) << std::left << "---";
    os << std::setw(15) << std::left << state_->searchSize;
    os << std::setw(10) << std::left << "---";
    os << std::endl;
  }
  if ( step_flag_==1 && state_->iter > 0 ) {
    if ( print_header ) {
      writeHeader(os);
    }
    else {
      os << "  ";
      os << std::setw(6) << std::left << state_->iter;
      os << std::setw(15) << std::left << state_->value;
      os << std::setw(15) << std::left << state_->gnorm;
      os << std::setw(15) << std::left << state_->snorm;
      os << std::setw(10) << std::left << state_->nfval;
      os << std::setw(10) << std::left << state_->ngrad;
      os << std::setw(15) << std::left << state_->aggregateGradientNorm;
      os << std::setw(15) << std::left << state_->aggregateModelError;
      os << std::setw(15) << std::left << state_->searchSize;
      os << std::setw(10) << std::left << QPiter_;
      os << std::endl;
    }
  }
  os.flags(osFlags);
}

} // namespace TypeU
} // namespace ROL

#endif
