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
    \brief Implements the computation of optimization steps using Moreau-Yosida
           regularized bound constraints.

    To describe the generalized Moreau-Yosida penalty method, we consider the
    following abstract setting.  Suppose \f$\mathcal{X}\f$ is a Hilbert space
    of functions mapping \f$\Xi\f$ to \f$\mathbb{R}\f$.  For example, 
    \f$\Xi\subset\mathbb{R}^n\f$ and \f$\mathcal{X}=L^2(\Xi)\f$ or 
    \f$\Xi = \{1,\ldots,n\}\f$ and \f$\mathcal{X}=\mathbb{R}^n\f$. We assume
    \f$ f:\mathcal{X}\to\mathbb{R}\f$ is twice-continuously Fr&eacute;chet 
    differentiable and \f$a,\,b\in\mathcal{X}\f$ with \f$a\le b\f$ almost 
    everywhere in \f$\Xi\f$.  Note that the generalized Moreau-Yosida penalty
    method will also work with secant approximations of the Hessian. 

    The generalized Moreau-Yosida penalty method is a proveably convergent
    algorithm for convex optimization problems and may not converge for general
    nonlinear, nonconvex problems.  The algorithm solves
    \f[
       \min_x \quad f(x) \quad \text{s.t.} \quad c(x) = 0, \quad a \le x \le b.
    \f]
    We can respresent the bound constraints using the indicator function
    \f$\iota_{[a,b]}(x) = 0\f$ if \f$a \le x \le b\f$ and equals \f$\infty\f$
    otherwise.  Using this indicator function, we can write our optimization
    problem as the (nonsmooth) equality constrained program
    \f[
       \min_x \quad f(x) + \iota_{[a,b]}(x) \quad \text{s.t.}\quad c(x) = 0.
    \f]
    Since the indicator function is not continuously Fr&eacute;chet
    differentiable, we cannot apply our existing algorithms (such as, Composite
    Step SQP) to the above equality constrained problem.  To circumvent this
    issue, we smooth the indicator function using generalized Moreau-Yosida
    regularization, i.e., we replace \f$\iota_{[a,b]}\f$ in the objective
    function with
    \f[
       \varphi(x,\mu,c) = \inf_y\; \{\; \iota_{[a,b]}(x-y)
         + \langle \mu, y\rangle_{\mathcal{X}}
         + \frac{c}{2}\|y\|_{\mathcal{X}}^2 \;\}.
    \f]
    One can show that \f$\varphi(\cdot,\mu,c)\f$ for any \f$\mu\in\mathcal{X}\f$
    and \f$c > 0\f$ is continuously Fr&eacute;chet
    differentiable with respect to \f$x\f$.  Thus, using this penalty,
    Step::compute solves the following subproblem: given
    \f$c_k>0\f$ and \f$\mu_k\in\mathcal{X}\f$, determine \f$x_k\in\mathcal{X}\f$
    that solves
    \f[
      \min_{x} \quad f(x) + \varphi(x,\mu_k,c_k)\quad\text{s.t.}
         c(x) = 0.
    \f]
    The multipliers \f$\mu_k\f$ are then updated in Step::update as
    \f$\mu_{k+1} = \nabla_x\varphi(x_k,\mu_k,c_k)\f$ and \f$c_k\f$ is
    potentially increased (although this is not always necessary).

    For more information on this method see:
    \li D. P. Bertsekas. "Approximations Procedures Based on the Method of
    Multipliers." Journal of Optimization Theory and Applications,
    Vol. 23(4), 1977.
    \li K. Ito, K. Kunisch. "Augmented Lagrangian Methods for Nonsmooth,
    Convex, Optimization in Hilbert Space." Nonlinear Analysis, 2000.
*/


namespace ROL {

template <class Real>
class MoreauYosidaPenaltyStep : public Step<Real> {
private:
  Teuchos::RCP<Algorithm<Real> > algo_;
  Teuchos::RCP<Vector<Real> > x_; 
  Teuchos::RCP<Vector<Real> > g_; 
  Teuchos::RCP<Vector<Real> > l_; 

  Real tau_;
  bool print_;

  Teuchos::ParameterList parlist_;
  int subproblemIter_;

  void updateState(const Vector<Real> &x, const Vector<Real> &l,
                   Objective<Real> &obj,
                   EqualityConstraint<Real> &con, BoundConstraint<Real> &bnd,
                   AlgorithmState<Real> &algo_state) {
    MoreauYosidaPenalty<Real> &myPen
      = Teuchos::dyn_cast<MoreauYosidaPenalty<Real> >(obj);
    Real zerotol = std::sqrt(ROL_EPSILON<Real>()), one(1);
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    // Update objective and constraint.
    myPen.update(x,true,algo_state.iter);
    con.update(x,true,algo_state.iter);
    // Compute objective value, constraint value, & gradient of Lagrangian
    algo_state.value = myPen.value(x, zerotol);
    con.value(*(state->constraintVec),x, zerotol);
    myPen.gradient(*(state->gradientVec), x, zerotol);
    con.applyAdjointJacobian(*g_,l,x,zerotol);
    state->gradientVec->plus(*g_);
    // Compute criticality measure
    if (bnd.isActivated()) {
      x_->set(x);
      x_->axpy(-one,(state->gradientVec)->dual());
      bnd.project(*x_);
      x_->axpy(-one,x);
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

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

  ~MoreauYosidaPenaltyStep() {}

  MoreauYosidaPenaltyStep(Teuchos::ParameterList &parlist)
    : Step<Real>(), algo_(Teuchos::null),
      x_(Teuchos::null), g_(Teuchos::null), l_(Teuchos::null),
      tau_(10), print_(false), parlist_(parlist), subproblemIter_(0) {
    // Parse parameters
    Real ten(10), oem6(1.e-6), oem8(1.e-8);
    Teuchos::ParameterList& steplist = parlist.sublist("Step").sublist("Moreau-Yosida Penalty");
    Step<Real>::getState()->searchSize = steplist.get("Initial Penalty Parameter",ten);
    tau_   = steplist.get("Penalty Parameter Growth Factor",ten);
    print_ = steplist.sublist("Subproblem").get("Print History",false);
    // Set parameters for step subproblem
    Real gtol = steplist.sublist("Subproblem").get("Optimality Tolerance",oem8);
    Real ctol = steplist.sublist("Subproblem").get("Feasibility Tolerance",oem8);
    Real stol = oem6*std::min(gtol,ctol);
    int maxit = steplist.sublist("Subproblem").get("Iteration Limit",1000);
    parlist_.sublist("Status Test").set("Gradient Tolerance",   gtol);
    parlist_.sublist("Status Test").set("Constraint Tolerance", ctol);
    parlist_.sublist("Status Test").set("Step Tolerance",       stol);
    parlist_.sublist("Status Test").set("Iteration Limit",      maxit);
  }

  /** \brief Initialize step with equality constraint.
  */
  void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                   Objective<Real> &obj, EqualityConstraint<Real> &con, BoundConstraint<Real> &bnd,
                   AlgorithmState<Real> &algo_state ) {
    MoreauYosidaPenalty<Real> &myPen
      = Teuchos::dyn_cast<MoreauYosidaPenalty<Real> >(obj);
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
    myPen.updateMultipliers(state->searchSize,x);
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
    Real one(1);
    MoreauYosidaPenalty<Real> &myPen
      = Teuchos::dyn_cast<MoreauYosidaPenalty<Real> >(obj);
    algo_ = Teuchos::rcp(new Algorithm<Real>("Composite Step",parlist_,false));
    x_->set(x); l_->set(l);
    algo_->run(*x_,*l_,myPen,con,print_);
    s.set(*x_); s.axpy(-one,x);
    subproblemIter_ = (algo_->getState())->iter;
  }

  /** \brief Update step, if successful (equality and bound constraints).
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, EqualityConstraint<Real> &con,
               BoundConstraint<Real> &bnd,
               AlgorithmState<Real> &algo_state ) {
    MoreauYosidaPenalty<Real> &myPen
      = Teuchos::dyn_cast<MoreauYosidaPenalty<Real> >(obj);
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    state->descentVec->set(s);
    // Update iterate and Lagrange multiplier
    x.plus(s);
    l.set(*l_);
    // Update objective and constraint
    algo_state.iter++;
    con.update(x,true,algo_state.iter);
    myPen.update(x,true,algo_state.iter);
    // Update multipliers
    state->searchSize *= tau_;
    myPen.updateMultipliers(state->searchSize,x);
    // Update state
    updateState(x,l,obj,con,bnd,algo_state);
    algo_state.nfval += myPen.getNumberFunctionEvaluations() + ((algo_->getState())->nfval);
    algo_state.ngrad += myPen.getNumberGradientEvaluations() + ((algo_->getState())->ngrad);
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
    hist << std::setw(10) << std::left << "penalty";
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
      hist << std::scientific << std::setprecision(2);
      hist << std::setw(10) << std::left << Step<Real>::getStepState()->searchSize;
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
