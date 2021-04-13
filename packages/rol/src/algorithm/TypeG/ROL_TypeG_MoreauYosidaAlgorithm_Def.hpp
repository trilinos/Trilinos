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

#ifndef ROL_TYPEG_MOREAUYOSIDAALGORITHM_DEF_H
#define ROL_TYPEG_MOREAUYOSIDAALGORITHM_DEF_H

#include "ROL_TypeE_AlgorithmFactory.hpp"

namespace ROL {
namespace TypeG {

template<typename Real>
MoreauYosidaAlgorithm<Real>::MoreauYosidaAlgorithm(ParameterList &list)
  : TypeG::Algorithm<Real>::Algorithm(),
    tau_(10), print_(false), list_(list), subproblemIter_(0) {
  // Set status test
  status_->reset();
  status_->add(makePtr<ConstraintStatusTest<Real>>(list));

  // Parse parameters
  Real ten(10), oem6(1.e-6), oem8(1.e-8), oe8(1e8);
  ParameterList& steplist = list.sublist("Step").sublist("Moreau-Yosida Penalty");
  state_->searchSize = steplist.get("Initial Penalty Parameter",           ten);
  maxPenalty_        = steplist.get("Maximum Penalty Parameter",           oe8);
  tau_               = steplist.get("Penalty Parameter Growth Factor",     ten);
  updatePenalty_     = steplist.get("Update Penalty",                      true);
  updateMultiplier_  = steplist.get("Update Multiplier",                   true);
  print_             = steplist.sublist("Subproblem").get("Print History", false);
  // Set parameters for step subproblem
  Real gtol = steplist.sublist("Subproblem").get("Optimality Tolerance",  oem8);
  Real ctol = steplist.sublist("Subproblem").get("Feasibility Tolerance", oem8);
  int maxit = steplist.sublist("Subproblem").get("Iteration Limit",       1000);
  Real stol = oem6*std::min(gtol,ctol);
  list_.sublist("Status Test").set("Gradient Tolerance",   gtol);
  list_.sublist("Status Test").set("Constraint Tolerance", ctol);
  list_.sublist("Status Test").set("Step Tolerance",       stol);
  list_.sublist("Status Test").set("Iteration Limit",      maxit);
  // Get step name from parameterlist
  stepname_ = steplist.sublist("Subproblem").get("Step Type","Augmented Lagrangian");

  // Output settings
  verbosity_    = list.sublist("General").get("Output Level", 0);
  printHeader_  = verbosity_ > 2;
  print_              = (verbosity_ > 2 ? true : print_);
  list_.sublist("General").set("Output Level",(print_ ? verbosity_ : 0));
}

template<typename Real>
void MoreauYosidaAlgorithm<Real>::initialize(Vector<Real>                &x,
                                             const Vector<Real>          &g,
                                             const Vector<Real>          &l,
                                             const Vector<Real>          &c,
                                             MoreauYosidaObjective<Real> &myobj,
                                             BoundConstraint<Real>       &bnd,
                                             Constraint<Real>            &con,
                                             Vector<Real>                &pwa,
                                             Vector<Real>                &dwa,
                                             std::ostream                &outStream) {
  hasPolyProj_ = true;
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
    hasPolyProj_ = false;
  }
  proj_->project(x,outStream);
  // Initialize data
  TypeG::Algorithm<Real>::initialize(x,g,l,c);
  // Initialize the algorithm state
  state_->nfval = 0;
  state_->ngrad = 0;
  state_->ncval = 0;
  updateState(x,l,myobj,bnd,con,pwa,dwa,outStream);
}


template<typename Real>
void MoreauYosidaAlgorithm<Real>::updateState(const Vector<Real>          &x,
                                              const Vector<Real>          &l,
                                              MoreauYosidaObjective<Real> &myobj,
                                              BoundConstraint<Real>       &bnd,
                                              Constraint<Real>            &con,
                                              Vector<Real>                &pwa,
                                              Vector<Real>                &dwa,
                                              std::ostream                &outStream) {
  const Real one(1);
  Real zerotol = std::sqrt(ROL_EPSILON<Real>());
  // Update objective and constraint
  if (state_->iter == 0) {
    myobj.update(x,UpdateType::Initial,state_->iter);
    con.update(x,UpdateType::Initial,state_->iter);
  }
  //else {
  //  myobj.update(x,UpdateType::Accept,state_->iter);
  //  con.update(x,UpdateType::Accept,state_->iter);
  //}
  // Compute norm of the gradient of the Lagrangian
  state_->value = myobj.getObjectiveValue(x, zerotol);
  myobj.getObjectiveGradient(*state_->gradientVec, x, zerotol);
  //myobj.gradient(*state_->gradientVec, x, zerotol);
  con.applyAdjointJacobian(dwa, l, x, zerotol);
  state_->gradientVec->plus(dwa);
  //gnorm_ = state_->gradientVec->norm();
  pwa.set(x);
  pwa.axpy(-one,state_->gradientVec->dual());
  proj_->project(pwa,outStream);
  pwa.axpy(-one,x);
  gnorm_ = pwa.norm();
  // Compute constraint violation
  con.value(*state_->constraintVec, x, zerotol);
  state_->cnorm = state_->constraintVec->norm();
  compViolation_ = myobj.testComplementarity(x);
  state_->gnorm = std::max(gnorm_,compViolation_);
  // Update state
  state_->nfval++;
  state_->ngrad++;
  state_->ncval++;
}

template<typename Real>
void MoreauYosidaAlgorithm<Real>::run( Vector<Real>          &x,
                                       const Vector<Real>    &g, 
                                       Objective<Real>       &obj,
                                       BoundConstraint<Real> &bnd,
                                       Constraint<Real>      &econ,
                                       Vector<Real>          &emul,
                                       const Vector<Real>    &eres,
                                       std::ostream          &outStream ) {
  const Real one(1);
  Ptr<Vector<Real>> pwa = x.clone(), dwa = g.clone();
  // Initialize Moreau-Yosida data
  MoreauYosidaObjective<Real> myobj(makePtrFromRef(obj),makePtrFromRef(bnd),
                                  x,g,state_->searchSize,updateMultiplier_,
                                  updatePenalty_);
  initialize(x,g,emul,eres,myobj,bnd,econ,*pwa,*dwa,outStream);
  Ptr<TypeE::Algorithm<Real>> algo;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Solve augmented Lagrangian subproblem
    algo = TypeE::AlgorithmFactory<Real>(list_);
    emul.zero();
    if (hasPolyProj_) algo->run(x,g,myobj,econ,emul,eres,
                                *proj_->getLinearConstraint(),
                                *proj_->getMultiplier(),
                                *proj_->getResidual(),outStream);
    else              algo->run(x,g,myobj,econ,emul,eres,outStream);
    subproblemIter_ = algo->getState()->iter;
    state_->nfval += algo->getState()->nfval;
    state_->ngrad += algo->getState()->ngrad;
    state_->ncval += algo->getState()->ncval;

    // Compute step
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();
    state_->lagmultVec->axpy(-one,emul);
    state_->snorm += state_->lagmultVec->norm();

    // Update iterate and Lagrange multiplier
    state_->iterateVec->set(x);
    state_->lagmultVec->set(emul);

    // Update objective and constraint
    state_->iter++;

    // Update state
    updateState(x,emul,myobj,bnd,econ,*pwa,*dwa);

    // Update multipliers
    if (updatePenalty_) {
      state_->searchSize *= tau_;
      state_->searchSize  = std::min(state_->searchSize,maxPenalty_);
    }
    myobj.updateMultipliers(state_->searchSize,x);

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,printHeader_);
  }
  if (verbosity_ > 0) TypeG::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void MoreauYosidaAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(109,'-') << std::endl;
    hist << "Moreau-Yosida Penalty Solver";
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  fval     - Objective function value" << std::endl;
    hist << "  cnorm    - Norm of the constraint" << std::endl;
    hist << "  gLnorm   - Norm of the gradient of the Lagrangian" << std::endl;
    hist << "  ifeas    - Infeasibility metric" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  penalty  - Penalty parameter for bound constraints" << std::endl;
    hist << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << "  #cval    - Cumulative number of times the constraint was evaluated" << std::endl;
    hist << "  subiter  - Number of subproblem iterations" << std::endl;
    hist << std::string(109,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "fval";
  hist << std::setw(15) << std::left << "cnorm";
  hist << std::setw(15) << std::left << "gLnorm";
  hist << std::setw(15) << std::left << "ifeas";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(10) << std::left << "penalty";
  hist << std::setw(8) << std::left << "#fval";
  hist << std::setw(8) << std::left << "#grad";
  hist << std::setw(8) << std::left << "#cval";
  hist << std::setw(8) << std::left << "subIter";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void MoreauYosidaAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "Moreau-Yosida Penalty Solver (Type G, General Constraints)";
  hist << std::endl;
  hist << "Subproblem Solver: " << stepname_ << std::endl;
  os << hist.str();
}

template<typename Real>
void MoreauYosidaAlgorithm<Real>::writeOutput( std::ostream& os, const bool print_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( print_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->cnorm;
    hist << std::setw(15) << std::left << gnorm_;
    hist << std::setw(15) << std::left << compViolation_;
    hist << std::setw(15) << std::left << "---";
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << state_->searchSize;
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::setw(8) << std::left << state_->ncval;
    hist << std::setw(8) << std::left << "---";
    hist << std::endl;
  }
  else {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->cnorm;
    hist << std::setw(15) << std::left << gnorm_;
    hist << std::setw(15) << std::left << compViolation_;
    hist << std::setw(15) << std::left << state_->snorm;
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << state_->searchSize;
    hist << std::scientific << std::setprecision(6);
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::setw(8) << std::left << state_->ncval;
    hist << std::setw(8) << std::left << subproblemIter_;
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeG
} // namespace ROL

#endif
