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
  Teuchos::RCP<Step<Real> > step_;
  Teuchos::RCP<StatusTest<Real> > status_;
  Teuchos::RCP<DefaultAlgorithm<Real> > algo_;
  Teuchos::RCP<Vector<Real> > x_; 

  Teuchos::RCP<Teuchos::ParameterList> parlist_;

  Real tau_;
  Real alpha1_;
  Real alpha2_;
  Real beta1_;
  Real beta2_;
  Real eta0_;
  Real eta1_;
  Real omega0_;
  Real omega1_;
  Real gamma1_;

  bool print_;

  Real eta_;
  Real omega_;
  Real gamma_;

  int maxit_;
  int subproblemIter_;
  int useTR_;

public:
  ~AugmentedLagrangianStep() {}

  AugmentedLagrangianStep(Teuchos::ParameterList &parlist)
    : Step<Real>(), augLag_(Teuchos::null),
      step_(Teuchos::null), status_(Teuchos::null), algo_(Teuchos::null),
      x_(Teuchos::null), parlist_(Teuchos::null),
      tau_(1.e2), alpha1_(1.), alpha2_(0.1), beta1_(1.), beta2_(0.9),
      eta0_(1.), eta1_(1.e-8), omega0_(1.), omega1_(1.e-8), gamma1_(0.1),
      print_(false), eta_(0.), omega_(0.), gamma_(0.),
      maxit_(1000), subproblemIter_(0), useTR_(0) {
    Teuchos::ParameterList& sublist = parlist.sublist("Step").sublist("Augmented Lagrangian");
    Step<Real>::getState()->searchSize = sublist.get("Initial Penalty Parameter",1.e1);
    tau_    = sublist.get("Penalty Parameter Growth Factor",         1.e2);
    alpha1_ = sublist.get("Optimality Tolerance Update Exponent",    1.0);
    alpha2_ = sublist.get("Feasibility Tolerance Update Exponent",   0.1);
    beta1_  = sublist.get("Optimality Tolerance Decrease Exponent",  1.0);
    beta2_  = sublist.get("Feasibility Tolerance Decrease Exponent", 0.9);
    eta0_   = sublist.get("Initial Optimality Tolerance",            1.0);
    omega0_ = sublist.get("Initial Feasibility Tolerance",           1.0);
    gamma1_ = sublist.get("Minimum Penalty Parameter Reciprocal",    0.1);
    print_  = sublist.get("Print Intermediate Optimization History", false);
    maxit_  = sublist.get("Subproblem Iteration Limit",              1000);
    useTR_  = sublist.get("Subproblem Step Type",                    0);

    eta1_   = parlist.sublist("Status Test").get("Gradient Tolerance",   1.e-8);
    omega1_ = sublist.sublist("Status Test").get("Constraint Tolerance", 1.e-8);

    parlist_ = Teuchos::rcp(&parlist,false);
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
    // Initialize intermediate stopping tolerances
    state->searchSize = 10.0;
    gamma_ = std::min(1.0/state->searchSize,gamma1_);
    omega_ = omega0_*std::pow(gamma_,alpha1_);
    eta_   = eta0_*std::pow(gamma_,alpha2_);
    // Project x onto the feasible set
    if ( bnd.isActivated() ) {
      bnd.project(x);
    }
    // Update the Lagrangian
    augLag_ = Teuchos::rcp(new AugmentedLagrangian<Real>(obj,con,x,c));
    augLag_->updateMultipliers(l,state->searchSize);
    // Initialize the algorithm state
    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;
    // Initialize additional storage
    x_ = x.clone();
    // Update objective and constraint.
    Real zerotol = 0.0;
    augLag_->update(x,true,algo_state.iter);
    algo_state.value = augLag_->value(x, zerotol);
    algo_state.value = augLag_->getObjectiveValue();
    algo_state.nfval += augLag_->getNumberFunctionEvaluations();
    augLag_->gradient(*(state->gradientVec), x, zerotol);
    algo_state.ngrad += augLag_->getNumberGradientEvaluations();
    if ( bnd.isActivated() ) {
      x_->set(x);
      x_->axpy(-1.0,(state->gradientVec)->dual());
      bnd.project(*x_);
      x_->axpy(-1.0,x);
      algo_state.gnorm = x_->norm();
    }
    else {
      algo_state.gnorm = (state->gradientVec)->norm();
    }
    state->constraintVec = c.clone();
    augLag_->getConstraintVec(*(state->constraintVec),x);
    algo_state.ncval += augLag_->getNumberConstraintEvaluations();
    algo_state.cnorm = (state->constraintVec)->norm();
  }

  /** \brief Compute step (equality and bound constraints).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, EqualityConstraint<Real> &con, 
                BoundConstraint<Real> &bnd, 
                AlgorithmState<Real> &algo_state ) {
    Real tol = std::max(omega_,omega1_);
    if ( useTR_ == 0 ) {
      step_ = Teuchos::rcp(new TrustRegionStep<Real>(*parlist_));
    }
    else {
      step_ = Teuchos::rcp(new LineSearchStep<Real>(*parlist_));
    }
    status_ = Teuchos::rcp(new StatusTest<Real>(tol,1.e-6*tol,maxit_));
    algo_   = Teuchos::rcp(new DefaultAlgorithm<Real>(*step_,*status_,false));
    x_->set(x);
    algo_->run(*x_,*augLag_,bnd,print_);
    s.set(*x_); s.axpy(-1.0,x);
    subproblemIter_ = (algo_->getState())->iter;
  }

  /** \brief Update step, if successful (equality and bound constraints).
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, EqualityConstraint<Real> &con,
               BoundConstraint<Real> &bnd,
               AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    state->gradientVec->set(*((step_->getStepState())->gradientVec));
    augLag_->getConstraintVec(*(state->constraintVec),x);
    
    x.plus(s);
    algo_state.iter++;
    augLag_->update(x,true,algo_state.iter);
    bnd.update(x,true,algo_state.iter);

    if (bnd.isActivated()) {
      x_->set(x);
      x_->axpy(-1.0,(state->gradientVec)->dual());
      bnd.project(*x_);
      x_->axpy(-1.0,x);
      algo_state.gnorm = x_->norm();
    }
    else {
      algo_state.gnorm = (state->gradientVec)->norm();
    }
    algo_state.nfval += augLag_->getNumberFunctionEvaluations();
    algo_state.ngrad += augLag_->getNumberGradientEvaluations();
    algo_state.ncval += augLag_->getNumberConstraintEvaluations();
    algo_state.value = augLag_->getObjectiveValue();
    algo_state.snorm = s.norm();
    algo_state.cnorm = (state->constraintVec)->norm();
    algo_state.iterateVec->set(x);
    state->descentVec->set(s);
    //if ( algo_state.cnorm < std::max(eta_,eta1_) ) {
    if ( algo_state.cnorm < eta_ ) {
      l.axpy(-state->searchSize,(state->constraintVec)->dual());
      algo_state.snorm += state->searchSize*algo_state.cnorm;

      gamma_  = std::min(1.0/state->searchSize,gamma1_);
      omega_ *= std::pow(gamma_,beta1_);
      eta_   *= std::pow(gamma_,beta2_);
    }
    else {
      state->searchSize *= tau_;
      gamma_  = std::min(1.0/state->searchSize,gamma1_);
      omega_  = omega0_*std::pow(gamma_,alpha1_);
      eta_    = eta0_*std::pow(gamma_,alpha2_);
    }
    algo_state.lagmultVec->set(l);
    augLag_->updateMultipliers(l,state->searchSize);
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
    hist << std::setw(15) << std::left << "penalty";
    hist << std::setw(15) << std::left << "feasTol";
    hist << std::setw(15) << std::left << "optTol";
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
      hist << std::setw(15) << std::left << Step<Real>::getStepState()->searchSize;
      hist << std::setw(15) << std::left << std::max(eta_,eta1_);
      hist << std::setw(15) << std::left << std::max(omega_,omega1_);
      hist << "\n";
    }
    else {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.snorm;
      hist << std::setw(15) << std::left << Step<Real>::getStepState()->searchSize;
      hist << std::setw(15) << std::left << eta_;
      hist << std::setw(15) << std::left << omega_;
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
