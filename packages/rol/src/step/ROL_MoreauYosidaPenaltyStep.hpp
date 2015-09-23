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

#ifndef ROL_MOREAUYOSIDAPENALTYSTEP_H
#define ROL_MOREAUYOSIDAPENALTYSTEP_H

#include "ROL_MoreauYosidaPenalty.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_Types.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_StatusTestSQP.hpp"
#include "ROL_CompositeStepSQP.hpp"
#include "Teuchos_ParameterList.hpp"

/** @ingroup step_group
    \class ROL::MoreauYosidaPenaltyStep
    \brief Provides the interface to compute augmented Lagrangian steps.
*/


namespace ROL {

template <class Real>
class MoreauYosidaPenaltyStep : public Step<Real> {
private:
  Teuchos::RCP<MoreauYosidaPenalty<Real> > myPen_;
  Teuchos::RCP<Step<Real> > step_;
  Teuchos::RCP<StatusTest<Real> > status_;
  Teuchos::RCP<DefaultAlgorithm<Real> > algo_;
  Teuchos::RCP<Vector<Real> > x_; 
  Teuchos::RCP<Vector<Real> > l_; 

  Real tau_;
  Real eta_;
  Real omega_;
  bool print_;
  int maxit_;

  Teuchos::RCP<Teuchos::ParameterList> parlist_;
  int subproblemIter_;

public:
  ~MoreauYosidaPenaltyStep() {}

  MoreauYosidaPenaltyStep(Teuchos::ParameterList &parlist)
    : Step<Real>(), myPen_(Teuchos::null),
      step_(Teuchos::null), status_(Teuchos::null), algo_(Teuchos::null),
      x_(Teuchos::null), l_(Teuchos::null),
      tau_(10.), eta_(1.e-8), omega_(1.e-8), print_(false), maxit_(1000),
      parlist_(Teuchos::null), subproblemIter_(0) {
    Teuchos::ParameterList& steplist = parlist.sublist("Step").sublist("Moreau-Yosida Penalty");
    Step<Real>::getState()->searchSize = steplist.get("Initial Penalty Parameter",10.0);
    tau_   = steplist.get("Penalty Parameter Growth Factor",10.0);
    eta_   = steplist.sublist("Subproblem").get("Optimality Tolerance",1.e-8);
    omega_ = steplist.sublist("Subproblem").get("Feasibility Tolerance",1.e-8);
    maxit_ = steplist.sublist("Subproblem").get("Iteration Limit",1000);
    print_ = steplist.sublist("Subproblem").get("Print History",false);

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
    // Project x onto the feasible set
    if ( bnd.isActivated() ) {
      bnd.project(x);
    }
    // Update the Lagrangian
    myPen_ = Teuchos::rcp(new MoreauYosidaPenalty<Real>(obj,bnd,x,state->searchSize));
    myPen_->updateMultipliers(state->searchSize,x);
    // Initialize the algorithm state
    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;
    // Initialize additional storage
    x_ = x.clone();
    l_ = l.clone();
    // Update objective and constraint.
    Real zerotol = 0.0;
    myPen_->update(x,true,algo_state.iter);
    algo_state.value = myPen_->value(x, zerotol);
    algo_state.value = myPen_->getObjectiveValue();
    algo_state.nfval += myPen_->getNumberFunctionEvaluations();
    myPen_->gradient(*(state->gradientVec), x, zerotol);
    algo_state.ngrad += myPen_->getNumberGradientEvaluations();
    algo_state.gnorm = (state->gradientVec)->norm();
    state->constraintVec = c.clone();
    con.value(*(state->constraintVec),x,zerotol);
    algo_state.ncval++;
    algo_state.cnorm = (state->constraintVec)->norm();
  }

  /** \brief Compute step (equality and bound constraints).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, EqualityConstraint<Real> &con, 
                BoundConstraint<Real> &bnd, 
                AlgorithmState<Real> &algo_state ) {
    step_   = Teuchos::rcp(new CompositeStepSQP<Real>(*parlist_));
    status_ = Teuchos::rcp(new StatusTestSQP<Real>(eta_,omega_,1.e-6*eta_,maxit_));
    algo_   = Teuchos::rcp(new DefaultAlgorithm<Real>(*step_,*status_,false));
    x_->set(x); l_->set(l);
    algo_->run(*x_,*l_,*myPen_,con,print_);
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
    state->descentVec->set(s);

    x.plus(s);
    l.set(*l_);

    algo_state.iter++;
    con.update(x,true,algo_state.iter);
    myPen_->update(x,true,algo_state.iter);

    state->searchSize *= tau_;
    myPen_->updateMultipliers(state->searchSize,x);

    state->gradientVec->set(*((step_->getStepState())->gradientVec));
    state->constraintVec->set(*((step_->getStepState())->constraintVec));

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

    algo_state.nfval += myPen_->getNumberFunctionEvaluations() + ((algo_->getState())->nfval);
    algo_state.ngrad += myPen_->getNumberGradientEvaluations() + ((algo_->getState())->ngrad);
    algo_state.ncval += (algo_->getState())->ncval;
    algo_state.value = myPen_->getObjectiveValue();
    algo_state.cnorm = (state->constraintVec)->norm();
    algo_state.snorm = s.norm();
    algo_state.iterateVec->set(x);
    algo_state.lagmultVec->set(l);
  }

  /** \brief Print iterate header.
  */
  std::string printHeader( void ) const {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "fval";
    hist << std::setw(15) << std::left << "cnorm";
    hist << std::setw(15) << std::left << "gnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(15) << std::left << "penalty";
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
    hist << "\n" << " Moreau-Yosida Penalty solver";
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

}; // class MoreauYosidaPenaltyStep

} // namespace ROL

#endif
