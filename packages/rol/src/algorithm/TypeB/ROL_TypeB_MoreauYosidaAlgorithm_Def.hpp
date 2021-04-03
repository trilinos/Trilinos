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

#ifndef ROL_TYPEB_MOREAUYOSIDAALGORITHM_DEF_HPP
#define ROL_TYPEB_MOREAUYOSIDAALGORITHM_DEF_HPP

#include "ROL_TypeU_AlgorithmFactory.hpp"

namespace ROL {
namespace TypeB {

template<typename Real>
MoreauYosidaAlgorithm<Real>::MoreauYosidaAlgorithm(ParameterList &list)
  : TypeB::Algorithm<Real>::Algorithm(),
    tau_(10), print_(false), list_(list), subproblemIter_(0) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

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
  stepname_ = steplist.sublist("Subproblem").get("Step Type","Trust Region");

  // Output settings
  verbosity_    = list.sublist("General").get("Output Level", 0);
  writeHeader_  = verbosity_ > 2;
  print_              = (verbosity_ > 2 ? true : print_);
  list_.sublist("General").set("Output Level",(print_ ? verbosity_ : 0));
}

template<typename Real>
void MoreauYosidaAlgorithm<Real>::initialize(Vector<Real>                &x,
                                             const Vector<Real>          &g,
                                             MoreauYosidaObjective<Real> &myobj,
                                             BoundConstraint<Real>       &bnd,
                                             Vector<Real>                &pwa,
                                             std::ostream                &outStream) {
  hasEcon_ = true;
  if (proj_ == nullPtr) {
    proj_ = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
    hasEcon_ = false;
  }
  // Initialize data
  TypeB::Algorithm<Real>::initialize(x,g);
  // Initialize the algorithm state
  state_->nfval = 0;
  state_->ngrad = 0;
  updateState(x,myobj,bnd,pwa,outStream);
}


template<typename Real>
void MoreauYosidaAlgorithm<Real>::updateState(const Vector<Real>          &x,
                                              MoreauYosidaObjective<Real> &myobj,
                                              BoundConstraint<Real>       &bnd,
                                              Vector<Real>                &pwa,
                                              std::ostream                &outStream) {
  const Real one(1);
  Real zerotol = std::sqrt(ROL_EPSILON<Real>());
  // Update objective and constraint.
  if (state_->iter == 0) {
    myobj.update(x,UpdateType::Initial,state_->iter);
  }
  //else {
  //  myobj.update(x,UpdateType::Accept,state_->iter);
  //}
  // Compute norm of the gradient of the Lagrangian
  state_->value = myobj.getObjectiveValue(x, zerotol);
  myobj.getObjectiveGradient(*state_->gradientVec, x, zerotol);
  //myobj.gradient(*state_->gradientVec, x, zerotol);
  //gnorm_ = state_->gradientVec->norm();
  pwa.set(x);
  pwa.axpy(-one,state_->gradientVec->dual());
  proj_->project(pwa,outStream);
  pwa.axpy(-one,x);
  gnorm_ = pwa.norm();
  // Compute constraint violation
  compViolation_ = myobj.testComplementarity(x);
  state_->gnorm = std::max(gnorm_,compViolation_);
  // Update state
  state_->nfval++;
  state_->ngrad++;
}

template<typename Real>
void MoreauYosidaAlgorithm<Real>::run( Vector<Real>          &x,
                                       const Vector<Real>    &g,
                                       Objective<Real>       &obj,
                                       BoundConstraint<Real> &bnd,
                                       std::ostream          &outStream ) {
  const Real one(1);
  Ptr<Vector<Real>> pwa = x.clone();
  // Initialize Moreau-Yosida data
  MoreauYosidaObjective<Real> myobj(makePtrFromRef(obj),makePtrFromRef(bnd),
                                    x,g,state_->searchSize,updateMultiplier_,
                                    updatePenalty_);
  initialize(x,g,myobj,bnd,*pwa,outStream);
  Ptr<TypeU::Algorithm<Real>> algo;

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  while (status_->check(*state_)) {
    // Solve augmented Lagrangian subproblem
    algo = TypeU::AlgorithmFactory<Real>(list_);
    if (hasEcon_) algo->run(x,g,myobj,*proj_->getLinearConstraint(),
                            *proj_->getMultiplier(),*proj_->getResidual(),
                            outStream);
    else          algo->run(x,g,myobj,outStream);
    subproblemIter_ = algo->getState()->iter;

    // Compute step
    state_->stepVec->set(x);
    state_->stepVec->axpy(-one,*state_->iterateVec);
    state_->snorm = state_->stepVec->norm();

    // Update iterate and Lagrange multiplier
    state_->iterateVec->set(x);

    // Update objective and constraint
    state_->iter++;

    // Update state
    updateState(x,myobj,bnd,*pwa,outStream);

    // Update multipliers
    if (updatePenalty_) {
      state_->searchSize *= tau_;
      state_->searchSize  = std::min(state_->searchSize,maxPenalty_);
    }
    myobj.updateMultipliers(state_->searchSize,x);

    state_->nfval += myobj.getNumberFunctionEvaluations() + algo->getState()->nfval;
    state_->ngrad += myobj.getNumberGradientEvaluations() + algo->getState()->ngrad;

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeB::Algorithm<Real>::writeExitStatus(outStream);
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
    hist << "  gnorm    - Norm of the gradient" << std::endl;
    hist << "  ifeas    - Infeasibility metric" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  penalty  - Penalty parameter for bound constraints" << std::endl;
    hist << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << "  subiter  - Number of subproblem iterations" << std::endl;
    hist << std::string(109,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "fval";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "ifeas";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(10) << std::left << "penalty";
  hist << std::setw(8) << std::left << "#fval";
  hist << std::setw(8) << std::left << "#grad";
  hist << std::setw(8) << std::left << "subIter";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void MoreauYosidaAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << " Moreau-Yosida Penalty Solver";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void MoreauYosidaAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( write_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << gnorm_;
    hist << std::setw(15) << std::left << compViolation_;
    hist << std::setw(15) << std::left << "---";
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << state_->searchSize;
    hist << std::scientific << std::setprecision(6);
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::setw(8) << std::left << "---";
    hist << std::endl;
  }
  else {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << gnorm_;
    hist << std::setw(15) << std::left << compViolation_;
    hist << std::setw(15) << std::left << state_->snorm;
    hist << std::scientific << std::setprecision(2);
    hist << std::setw(10) << std::left << state_->searchSize;
    hist << std::scientific << std::setprecision(6);
    hist << std::setw(8) << std::left << state_->nfval;
    hist << std::setw(8) << std::left << state_->ngrad;
    hist << std::setw(8) << std::left << subproblemIter_;
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeB
} // namespace ROL

#endif
