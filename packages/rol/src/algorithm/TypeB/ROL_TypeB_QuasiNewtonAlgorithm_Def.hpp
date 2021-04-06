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

#ifndef ROL_TYPEB_QUASINEWTONALGORITHM_DEF_HPP
#define ROL_TYPEB_QUASINEWTONALGORITHM_DEF_HPP

#include "ROL_TypeB_GradientAlgorithm.hpp"
#include "ROL_TypeB_SpectralGradientAlgorithm.hpp"
#include "ROL_TypeB_LinMoreAlgorithm.hpp"
#include "ROL_TypeB_PrimalDualActiveSetAlgorithm.hpp"
#include "ROL_TypeB_MoreauYosidaAlgorithm.hpp"
#include "ROL_TypeB_InteriorPointAlgorithm.hpp"
#include "ROL_PQNObjective.hpp"

namespace ROL {
namespace TypeB {

template<typename Real>
QuasiNewtonAlgorithm<Real>::QuasiNewtonAlgorithm(ParameterList           &list,
                                                 const Ptr<Secant<Real>> &secant)
  : secant_(secant), esec_(SECANT_USERDEFINED), list_(list), hasLEC_(true) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  // Parse parameter list
  ParameterList &lslist = list.sublist("Step").sublist("Line Search");
  maxit_        = lslist.get("Function Evaluation Limit",                        20);
  c1_           = lslist.get("Sufficient Decrease Tolerance",                  1e-4);
  rhodec_       = lslist.sublist("Line-Search Method").get("Backtracking Rate", 0.5);
  sigma1_       = lslist.sublist("PQN").get("Lower Step Size Safeguard",        0.1);
  sigma2_       = lslist.sublist("PQN").get("Upper Step Size Safeguard",        0.9);
  algoName_     = lslist.sublist("PQN").get("Subproblem Solver","Spectral Gradient");
  int sp_maxit  = lslist.sublist("PQN").get("Subproblem Iteration Limit",      1000);
  sp_tol1_      = lslist.sublist("PQN").get("Subproblem Absolute Tolerance",   1e-4);
  sp_tol2_      = lslist.sublist("PQN").get("Subproblem Relative Tolerance",   1e-2);
  Real opt_tol  = lslist.sublist("Status Test").get("Gradient Tolerance",      1e-8);
  sp_tol_min_   = static_cast<Real>(1e-2)*opt_tol;
  verbosity_    = list.sublist("General").get("Output Level",                     0);
  writeHeader_  = verbosity_ > 2;

  list_.sublist("Status Test").set("Iteration Limit",    sp_maxit);
  list_.sublist("General").set("Output Level", verbosity_>0 ? verbosity_-1 : 0);

  if ( secant_ == nullPtr ) {
    secantName_ = list.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS");
    esec_ = StringToESecant(secantName_);
    secant_ = SecantFactory<Real>(list);
  }
  else {
    secantName_ = list.sublist("General").sublist("Secant").get("User Defined Secant Name",
                                                                "Unspecified User Defined Secant Method");
  }
}


template<typename Real>
void QuasiNewtonAlgorithm<Real>::initialize(Vector<Real>          &x,
                                            const Vector<Real>    &g,
                                            Objective<Real>       &obj,
                                            BoundConstraint<Real> &bnd,
                                            std::ostream &outStream) {
  const Real one(1);
  if (proj_ == nullPtr) {
    proj_   = makePtr<PolyhedralProjection<Real>>(makePtrFromRef(bnd));
    hasLEC_ = false;
  }
  // Initialize data
  TypeB::Algorithm<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  proj_->project(x,outStream); state_->nproj++;
  state_->iterateVec->set(x);
  obj.update(x,UpdateType::Initial,state_->iter);    
  state_->value = obj.value(x,ftol); state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
  state_->stepVec->set(x);
  state_->stepVec->axpy(-one,state_->gradientVec->dual());
  proj_->project(*state_->stepVec,outStream); state_->nproj++;
  state_->stepVec->axpy(-one,x);
  state_->gnorm = state_->stepVec->norm();
  state_->snorm = ROL_INF<Real>();
}

template<typename Real>
void QuasiNewtonAlgorithm<Real>::run( Vector<Real>          &x,
                                      const Vector<Real>    &g, 
                                      Objective<Real>       &obj,
                                      BoundConstraint<Real> &bnd,
                                      std::ostream          &outStream ) {
  const Real half(0.5), one(1);
  // Initialize trust-region data
  initialize(x,g,obj,bnd,outStream);
  Ptr<Vector<Real>> s = x.clone(), gp = x.clone(), gold = g.clone(), xs = x.clone();
  Real ftrial(0), gs(0), alphaTmp(0), tol(std::sqrt(ROL_EPSILON<Real>())), gtol(1);

  Ptr<TypeB::Algorithm<Real>> algo;
  Ptr<PQNObjective<Real>> qobj = makePtr<PQNObjective<Real>>(secant_,x,g);
  Ptr<Problem<Real>> problem = makePtr<Problem<Real>>(qobj,xs);
  problem->addBoundConstraint(makePtrFromRef(bnd));
  if (hasLEC_) {
    problem->addLinearConstraint("LEC",proj_->getLinearConstraint(),
                                       proj_->getMultiplier(),
                                       proj_->getResidual());
    problem->setProjectionAlgorithm(list_);
  }
  problem->finalize(false,verbosity_>2,outStream);

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  // Compute steepest descent step
  gp->set(state_->gradientVec->dual());
  while (status_->check(*state_)) {
    // Compute step
    qobj->setAnchor(x,*state_->gradientVec);
    xs->set(x); xs->axpy(-one,*gp); proj_->project(*xs,outStream); state_->nproj++;
    gtol = std::max(sp_tol_min_,std::min(sp_tol1_,sp_tol2_*state_->gnorm));
    list_.sublist("Status Test").set("Gradient Tolerance",gtol);
    if (algoName_ == "Trust Region")                algo = makePtr<TypeB::LinMoreAlgorithm<Real>>(list_);
    else if (algoName_ == "Line Search")            algo = makePtr<TypeB::GradientAlgorithm<Real>>(list_);
    else if (algoName_ == "Primal Dual Active Set") algo = makePtr<TypeB::PrimalDualActiveSetAlgorithm<Real>>(list_);
    else if (algoName_ == "Moreau-Yosida")          algo = makePtr<TypeB::MoreauYosidaAlgorithm<Real>>(list_);
    else if (algoName_ == "Interior Point")         algo = makePtr<TypeB::InteriorPointAlgorithm<Real>>(list_);
    else                                            algo = makePtr<TypeB::SpectralGradientAlgorithm<Real>>(list_);
    algo->run(*problem,outStream);
    s->set(*xs); s->axpy(-one,x);
    spgIter_ = algo->getState()->iter;
    state_->nproj += staticPtrCast<const TypeB::AlgorithmState<Real>>(algo->getState())->nproj;

    // Perform backtracking line search 
    state_->searchSize = one;
    x.set(*state_->iterateVec);
    x.axpy(state_->searchSize,*s);
    obj.update(x,UpdateType::Trial);
    ftrial = obj.value(x,tol); ls_nfval_ = 1;
    gs = state_->gradientVec->apply(*s);
    if (verbosity_ > 1) {
      outStream << "  In TypeB::QuasiNewtonAlgorithm: Line Search"  << std::endl;
      outStream << "    Step size:                        " << state_->searchSize         << std::endl;
      outStream << "    Trial objective value:            " << ftrial                     << std::endl;
      outStream << "    Computed reduction:               " << state_->value-ftrial       << std::endl;
      outStream << "    Dot product of gradient and step: " << gs                         << std::endl;
      outStream << "    Sufficient decrease bound:        " << -gs*state_->searchSize*c1_ << std::endl;
      outStream << "    Number of function evaluations:   " << ls_nfval_                  << std::endl;
    }
    while ( ftrial > state_->value + c1_*state_->searchSize*gs && ls_nfval_ < maxit_ ) {
      alphaTmp = -half*state_->searchSize*state_->searchSize*gs
                 / (ftrial-state_->value-state_->searchSize*gs);
      state_->searchSize = (sigma1_*state_->searchSize <= alphaTmp && alphaTmp <= sigma2_*state_->searchSize)
                            ? alphaTmp : rhodec_*state_->searchSize;
      //state_->searchSize *= rhodec_;
      x.set(*state_->iterateVec);
      x.axpy(state_->searchSize,*s);
      obj.update(x,UpdateType::Trial);
      ftrial = obj.value(x,tol); ls_nfval_++;
      if (verbosity_ > 1) {
        outStream << std::endl;
        outStream << "    Step size:                        " << state_->searchSize         << std::endl;
        outStream << "    Trial objective value:            " << ftrial                     << std::endl;
        outStream << "    Computed reduction:               " << state_->value-ftrial       << std::endl;
        outStream << "    Dot product of gradient and step: " << gs                         << std::endl;
        outStream << "    Sufficient decrease bound:        " << -gs*state_->searchSize*c1_ << std::endl;
        outStream << "    Number of function evaluations:   " << ls_nfval_                  << std::endl;
      }
    }
    state_->nfval += ls_nfval_;

    // Compute norm of step
    state_->stepVec->set(*s);
    state_->stepVec->scale(state_->searchSize);
    state_->snorm = state_->stepVec->norm();

    // Update iterate
    state_->iterateVec->set(x);

    // Compute new value and gradient
    state_->iter++;
    state_->value = ftrial;
    obj.update(x,UpdateType::Accept,state_->iter);
    gold->set(*state_->gradientVec);
    obj.gradient(*state_->gradientVec,x,tol); state_->ngrad++;
    gp->set(state_->gradientVec->dual());

    // Compute projected gradient norm
    s->set(x); s->axpy(-one,*gp);
    proj_->project(*s,outStream); state_->nproj++;
    s->axpy(-one,x);
    state_->gnorm = s->norm();

    // Update secant
    secant_->updateStorage(x,*state_->gradientVec,*gold,*state_->stepVec,state_->snorm,state_->iter);

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeB::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void QuasiNewtonAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(114,'-') << std::endl;
    hist << "Line-Search Projected Quasi-Newton with " << secantName_ << " Hessian approximation";
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  value    - Objective function value" << std::endl;
    hist << "  gnorm    - Norm of the gradient" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  alpha    - Line search step length" << std::endl;
    hist << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << "  #proj    - Cumulative number of times the projection was computed" << std::endl;
    hist << "  ls_#fval - Number of the times the objective function was evaluated during the line search" << std::endl;
    hist << "  sp_iter  - Number iterations to compute quasi-Newton step" << std::endl;
    hist << std::string(114,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(15) << std::left << "alpha";
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::setw(10) << std::left << "#proj";
  hist << std::setw(10) << std::left << "#ls_fval";
  hist << std::setw(10) << std::left << "sp_iter";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void QuasiNewtonAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "Line-Search Projected Quasi-Newton (Type B, Bound Constraints)" << std::endl;
  os << hist.str();
}

template<typename Real>
void QuasiNewtonAlgorithm<Real>::writeOutput( std::ostream& os, bool write_header ) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) writeName(os);
  if ( write_header )      writeHeader(os);
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << "---";
    hist << std::setw(15) << std::left << "---";
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << state_->nproj;
    hist << std::setw(10) << std::left << "---";
    hist << std::setw(10) << std::left << "---";
    hist << std::endl;
  }
  else {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << state_->snorm;
    hist << std::setw(15) << std::left << state_->searchSize;
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << state_->nproj;
    hist << std::setw(10) << std::left << ls_nfval_;
    hist << std::setw(10) << std::left << spgIter_;
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeB
} // namespace ROL

#endif
