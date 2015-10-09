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
  Teuchos::RCP<Algorithm<Real> > algo_;
  Teuchos::RCP<Vector<Real> > x_; 
  Teuchos::RCP<Vector<Real> > g_; 
  Teuchos::RCP<Vector<Real> > l_; 

  Real tau_;
  bool print_;

  Teuchos::RCP<Teuchos::ParameterList> parlist_;
  int subproblemIter_;

  void updateState(const Vector<Real> &x, const Vector<Real> &l,
                   Objective<Real> &obj,
                   EqualityConstraint<Real> &con, BoundConstraint<Real> &bnd,
                   AlgorithmState<Real> &algo_state) {
    Real zerotol = std::sqrt(ROL_EPSILON);
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    // Update objective and constraint.
    obj.update(x,true,algo_state.iter);
    con.update(x,true,algo_state.iter);
    myPen_->update(x,true,algo_state.iter);
    // Compute objective value, constraint value, & gradient of Lagrangian
    algo_state.value = myPen_->value(x, zerotol);
    con.value(*(state->constraintVec),x, zerotol);
    myPen_->gradient(*(state->gradientVec), x, zerotol);
    // Compute criticality measure
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
    algo_state.cnorm = (state->constraintVec)->norm();
    // Update state
    algo_state.nfval++;
    algo_state.ngrad++;
    algo_state.ncval++;
  }

public:
  ~MoreauYosidaPenaltyStep() {}

  MoreauYosidaPenaltyStep(Teuchos::ParameterList &parlist)
    : Step<Real>(), myPen_(Teuchos::null), algo_(Teuchos::null),
      x_(Teuchos::null), g_(Teuchos::null), l_(Teuchos::null),
      tau_(10.), print_(false), parlist_(Teuchos::null), subproblemIter_(0) {
    // Parse parameters
    Teuchos::ParameterList& steplist = parlist.sublist("Step").sublist("Moreau-Yosida Penalty");
    Step<Real>::getState()->searchSize = steplist.get("Initial Penalty Parameter",10.0);
    tau_   = steplist.get("Penalty Parameter Growth Factor",10.0);
    print_ = steplist.sublist("Subproblem").get("Print History",false);
    // Set parameters for step subproblem
    parlist_ = Teuchos::rcp(&parlist,false);
    Real gtol = steplist.sublist("Subproblem").get("Optimality Tolerance",1.e-8);
    Real ctol = steplist.sublist("Subproblem").get("Feasibility Tolerance",1.e-8);
    Real stol = 1.e-6*std::min(gtol,ctol);
    int maxit = steplist.sublist("Subproblem").get("Iteration Limit",1000);
    parlist_->sublist("Status Test").set("Gradient Tolerance",   gtol);
    parlist_->sublist("Status Test").set("Constraint Tolerance", ctol);
    parlist_->sublist("Status Test").set("Step Tolerance",       stol);
    parlist_->sublist("Status Test").set("Iteration Limit",      maxit);
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
    g_ = g.clone();
    l_ = l.clone();
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
    updateState(x,l,obj,con,bnd,algo_state);
  }

  /** \brief Compute step (equality and bound constraints).
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, EqualityConstraint<Real> &con, 
                BoundConstraint<Real> &bnd, 
                AlgorithmState<Real> &algo_state ) {
    algo_ = Teuchos::rcp(new Algorithm<Real>("Composite Step SQP",*parlist_,false));
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
    // Update iterate and Lagrange multiplier
    x.plus(s);
    l.set(*l_);
    // Update objective and constraint
    algo_state.iter++;
    con.update(x,true,algo_state.iter);
    myPen_->update(x,true,algo_state.iter);
    // Update multipliers
    state->searchSize *= tau_;
    myPen_->updateMultipliers(state->searchSize,x);
    // Update state
    updateState(x,l,obj,con,bnd,algo_state);
    algo_state.nfval += myPen_->getNumberFunctionEvaluations() + ((algo_->getState())->nfval);
    algo_state.ngrad += myPen_->getNumberGradientEvaluations() + ((algo_->getState())->ngrad);
    algo_state.ncval += (algo_->getState())->ncval;
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
