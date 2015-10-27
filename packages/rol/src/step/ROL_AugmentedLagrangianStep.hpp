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

#ifndef ROL_AUGMENTEDLAGRANGIANSTEP_H
#define ROL_AUGMENTEDLAGRANGIANSTEP_H

#include "ROL_AugmentedLagrangian.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_Types.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_StatusTest.hpp"
#include "ROL_Step.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "Teuchos_ParameterList.hpp"

/** @ingroup step_group
    \class ROL::AugmentedLagrangianStep
    \brief Provides the interface to compute augmented Lagrangian steps.
*/


namespace ROL {

template <class Real>
class AugmentedLagrangianStep : public Step<Real> {
private:
  Teuchos::RCP<AugmentedLagrangian<Real> > augLag_;
  Teuchos::RCP<Algorithm<Real> > algo_;
  Teuchos::RCP<Vector<Real> > x_; 

  Teuchos::ParameterList parlist_;
  // Lagrange multiplier update
  bool scaleLagrangian_;
  Real minPenaltyReciprocal_;
  Real minPenaltyLowerBound_;
  // Optimality tolerance update
  Real optIncreaseExponent_;
  Real optDecreaseExponent_;
  Real optToleranceInitial_;
  Real optTolerance_;
  // Feasibility tolerance update
  Real feasIncreaseExponent_;
  Real feasDecreaseExponent_;
  Real feasToleranceInitial_;
  Real feasTolerance_;
  // Subproblem information
  bool print_;
  int maxit_;
  int subproblemIter_;
  std::string subStep_;
  Real outerOptTolerance_;
  Real outerFeasTolerance_;

  Real computeGradient(Vector<Real> &g, const Vector<Real> &x,
                       const Real mu, BoundConstraint<Real> &bnd) {
    Real gnorm = 0., tol = std::sqrt(ROL_EPSILON);
    augLag_->gradient(g,x,tol);
    if ( scaleLagrangian_ ) {
      g.scale(mu);
    }
    // Compute norm of projected gradient
    if (bnd.isActivated()) {
      x_->set(x);
      x_->axpy(-1.,g.dual());
      bnd.project(*x_);
      x_->axpy(-1.,x);
      gnorm = x_->norm();
    }
    else {
      gnorm = g.norm();
    }
    return gnorm;
  }

public:
  ~AugmentedLagrangianStep() {}

  AugmentedLagrangianStep(Teuchos::ParameterList &parlist)
    : Step<Real>(), augLag_(Teuchos::null), algo_(Teuchos::null),
      x_(Teuchos::null), parlist_(parlist), subproblemIter_(0) {
    Teuchos::ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    Step<Real>::getState()->searchSize = sublist.get("Initial Penalty Parameter",1.e1);
    // Multiplier update parameters
    scaleLagrangian_      = sublist.get("Use Scaled Augmented Lagrangian",          false);
    minPenaltyLowerBound_ = sublist.get("Penalty Parameter Reciprocal Lower Bound", 0.1);
    minPenaltyReciprocal_ = 0.1;
    // Optimality tolerance update
    optIncreaseExponent_ = sublist.get("Optimality Tolerance Update Exponent",    1.0);
    optDecreaseExponent_ = sublist.get("Optimality Tolerance Decrease Exponent",  1.0);
    optToleranceInitial_ = sublist.get("Initial Optimality Tolerance",            1.0);
    // Feasibility tolerance update    
    feasIncreaseExponent_ = sublist.get("Feasibility Tolerance Update Exponent",   0.1);
    feasDecreaseExponent_ = sublist.get("Feasibility Tolerance Decrease Exponent", 0.9);
    feasToleranceInitial_ = sublist.get("Initial Feasibility Tolerance",           1.0);
    // Subproblem information
    print_   = sublist.get("Print Intermediate Optimization History", false);
    maxit_   = sublist.get("Subproblem Iteration Limit",              1000);
    subStep_ = sublist.get("Subproblem Step Type",                    "Trust Region");
    parlist_.sublist("Status Test").set("Iteration Limit",maxit_);
    // Outer iteration tolerances
    outerFeasTolerance_ = parlist.sublist("Status Test").get("Constraint Tolerance", 1.e-8);
    outerOptTolerance_  = parlist.sublist("Status Test").get("Gradient Tolerance", 1.e-8);
  }

  /** \brief Initialize step with equality constraint.
  */
  void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                   Objective<Real> &obj, EqualityConstraint<Real> &con, BoundConstraint<Real> &bnd,
                   AlgorithmState<Real> &algo_state ) {
    // Initialize step state
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    state->descentVec    = x.clone();
    state->gradientVec   = g.clone();
    state->constraintVec = c.clone();
    // Initialize additional storage
    x_ = x.clone();
    // Initialize the algorithm state
    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;
    // Initialize intermediate stopping tolerances
    minPenaltyReciprocal_ = std::min(1./state->searchSize,minPenaltyLowerBound_);
    optTolerance_  = optToleranceInitial_*std::pow(minPenaltyReciprocal_,optDecreaseExponent_);
    feasTolerance_ = feasToleranceInitial_*std::pow(minPenaltyReciprocal_,feasDecreaseExponent_);
    // Initialize the Augmented Lagrangian
    augLag_ = Teuchos::rcp(new AugmentedLagrangian<Real>(obj,con,x,c,l,state->searchSize,parlist_));
    // Project x onto the feasible set
    if ( bnd.isActivated() ) {
      bnd.project(x);
    }
    bnd.update(x,true,algo_state.iter);
    // Update objective and constraint.
    augLag_->update(x,true,algo_state.iter);
    algo_state.value = augLag_->getObjectiveValue(x);
    algo_state.gnorm = computeGradient(*(state->gradientVec),x,state->searchSize,bnd);
    augLag_->getConstraintVec(*(state->constraintVec),x);
    algo_state.cnorm = (state->constraintVec)->norm();
    // Update evaluation counters
    algo_state.ncval += augLag_->getNumberConstraintEvaluations();
    algo_state.nfval += augLag_->getNumberFunctionEvaluations();
    algo_state.ngrad += augLag_->getNumberGradientEvaluations();
  }

  /** \brief Compute step (equality and bound constraints).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, EqualityConstraint<Real> &con, 
                BoundConstraint<Real> &bnd, AlgorithmState<Real> &algo_state ) {
    parlist_.sublist("Status Test").set("Gradient Tolerance",optTolerance_);
    parlist_.sublist("Status Test").set("Step Tolerance",1.e-6*optTolerance_);
    algo_ = Teuchos::rcp(new Algorithm<Real>(subStep_,parlist_,false));
    x_->set(x);
    if ( bnd.isActivated() ) {
      algo_->run(*x_,*augLag_,bnd,print_);
    }
    else {
      algo_->run(*x_,*augLag_,print_);
    }
    s.set(*x_); s.axpy(-1.,x);
    subproblemIter_ = (algo_->getState())->iter;
  }

  /** \brief Update step, if successful (equality and bound constraints).
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, EqualityConstraint<Real> &con,
               BoundConstraint<Real> &bnd,
               AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    // Update the step and store in state
    x.plus(s);
    algo_state.iterateVec->set(x);
    state->descentVec->set(s);
    algo_state.snorm = s.norm();
    algo_state.iter++;
    // Update objective function and constraints
    augLag_->update(x,true,algo_state.iter);
    bnd.update(x,true,algo_state.iter);
    // Update multipliers
    bool updated = augLag_->updateMultipliers(l,state->searchSize,x,feasTolerance_);
    algo_state.snorm += (updated ? l.norm() + 1. : 0.);
    algo_state.lagmultVec->set(l);
    minPenaltyReciprocal_ = std::min(1./state->searchSize,minPenaltyLowerBound_);
    if ( algo_state.cnorm < feasTolerance_ ) {
      optTolerance_  *= std::pow(minPenaltyReciprocal_,optIncreaseExponent_);
      feasTolerance_ *= std::pow(minPenaltyReciprocal_,feasIncreaseExponent_);
    }
    else {
      optTolerance_  = optToleranceInitial_*std::pow(minPenaltyReciprocal_,optDecreaseExponent_);
      feasTolerance_ = feasToleranceInitial_*std::pow(minPenaltyReciprocal_,feasDecreaseExponent_);
    }
    // Update objective function value
    algo_state.value = augLag_->getObjectiveValue(x);
    // Update constraint value
    augLag_->getConstraintVec(*(state->constraintVec),x);
    algo_state.cnorm = (state->constraintVec)->norm();
    // Compute gradient of the augmented Lagrangian
    algo_state.gnorm = computeGradient(*(state->gradientVec),x,state->searchSize,bnd);
    // Update evaluation counters
    algo_state.nfval += augLag_->getNumberFunctionEvaluations();
    algo_state.ngrad += augLag_->getNumberGradientEvaluations();
    algo_state.ncval += augLag_->getNumberConstraintEvaluations();
    augLag_->reset();
  }

  /** \brief Print iterate header.
  */
  std::string printHeader( void ) const {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "fval";
    hist << std::setw(15) << std::left << "cnorm";
    hist << std::setw(15) << std::left << "gLnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(10) << std::left << "penalty";
    hist << std::setw(10) << std::left << "feasTol";
    hist << std::setw(10) << std::left << "optTol";
    hist << std::setw(8) << std::left << "#fval";
    hist << std::setw(8) << std::left << "#grad";
    hist << std::setw(8) << std::left << "#cval";
    hist << std::setw(8) << std::left << "subIter";
    hist << "\n";
    return hist.str();
  }

  /** \brief Print step name.
  */
  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << " Augmented Lagrangian solver";
    hist << "\n";
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> &algo_state, bool pHeader = false ) const {
    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << printName();
    }
    if ( pHeader ) {
      hist << printHeader();
    }
    if ( algo_state.iter == 0 ) {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << " ";
      hist << std::scientific << std::setprecision(2);
      hist << std::setw(10) << std::left << Step<Real>::getStepState()->searchSize;
      hist << std::setw(10) << std::left << std::max(feasTolerance_,outerFeasTolerance_);
      hist << std::setw(10) << std::left << std::max(optTolerance_,outerOptTolerance_);
      hist << "\n";
    }
    else {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.snorm;
      hist << std::scientific << std::setprecision(2);
      hist << std::setw(10) << std::left << Step<Real>::getStepState()->searchSize;
      hist << std::setw(10) << std::left << feasTolerance_;
      hist << std::setw(10) << std::left << optTolerance_;
      hist << std::scientific << std::setprecision(6);
      hist << std::setw(8) << std::left << algo_state.nfval;
      hist << std::setw(8) << std::left << algo_state.ngrad;
      hist << std::setw(8) << std::left << algo_state.ncval;
      hist << std::setw(8) << std::left << subproblemIter_;
      hist << "\n";
    }
    return hist.str();
  }

  /** \brief Compute step for bound constraints; here only to satisfy the
             interface requirements, does nothing, needs refactoring.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj,
                        BoundConstraint<Real> &con,
                        AlgorithmState<Real> &algo_state ) {}

  /** \brief Update step, for bound constraints; here only to satisfy the
             interface requirements, does nothing, needs refactoring.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj,
                       BoundConstraint<Real> &con,
                       AlgorithmState<Real> &algo_state ) {}

}; // class AugmentedLagrangianStep

} // namespace ROL

#endif
