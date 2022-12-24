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

#ifndef ROL_TYPEP_QUASINEWTONALGORITHM_DEF_HPP
#define ROL_TYPEP_QUASINEWTONALGORITHM_DEF_HPP

#include "ROL_TypeP_ProxGradientAlgorithm.hpp"
#include "ROL_TypeP_SpectralGradientAlgorithm.hpp"
#include "ROL_TypeP_iPianoAlgorithm.hpp"
#include "ROL_PQNObjective.hpp"

namespace ROL {
namespace TypeP {

template<typename Real>
QuasiNewtonAlgorithm<Real>::QuasiNewtonAlgorithm(ParameterList           &list,
                                                 const Ptr<Secant<Real>> &secant)
  : secant_(secant), esec_(SECANT_USERDEFINED), list_(list), hasLEC_(true) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(list));

  // Parse parameter list
  ParameterList &lslist = list.sublist("Step").sublist("Line Search");
  t0_           = list.sublist("Status Test").get("Gradient Scale"            , 1.0);
  initProx_     = lslist.get("Apply Prox to Initial Guess",                   false);
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
                                            Objective<Real>       &sobj,
                                            Objective<Real>       &nobj,
                                            Vector<Real>          &dg,
                                            std::ostream &outStream) {
  const Real one(1);
  Real tol(std::sqrt(ROL_EPSILON<Real>()));
  // Initialize data
  TypeP::Algorithm<Real>::initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  if (initProx_) {
    state_->iterateVec->set(x);
    nobj.prox(x,*state_->iterateVec,one,tol); state_->nprox++;
  }
  sobj.update(x,UpdateType::Initial,state_->iter);    
  nobj.update(x,UpdateType::Initial,state_->iter);    
  state_->svalue = sobj.value(x,ftol); state_->nsval++;
  state_->nvalue = nobj.value(x,ftol); state_->nnval++;
  state_->value  = state_->svalue + state_->nvalue;
  sobj.gradient(*state_->gradientVec,x,ftol); state_->ngrad++;
  dg.set(state_->gradientVec->dual());
  pgstep(*state_->iterateVec,*state_->stepVec,nobj,x,dg,t0_,tol);
  state_->gnorm = state_->stepVec->norm() / t0_;
  state_->snorm = ROL_INF<Real>();
}

template<typename Real>
void QuasiNewtonAlgorithm<Real>::run( Vector<Real>          &x,
                                      const Vector<Real>    &g, 
                                      Objective<Real>       &sobj,
                                      Objective<Real>       &nobj,
                                      std::ostream          &outStream ) {
  const Real half(0.5), one(1);
  // Initialize trust-region data
  Ptr<Vector<Real>> s = x.clone(), gp = x.clone(), gold = g.clone(), xs = x.clone();
  initialize(x,g,sobj,nobj,*gp,outStream);
  Real strial(0), ntrial(0), ftrial(0), gs(0), Qk(0), rhoTmp(0);
  Real tol(std::sqrt(ROL_EPSILON<Real>())), gtol(1);

  Ptr<TypeP::Algorithm<Real>> algo;
  Ptr<PQNObjective<Real>> qobj = makePtr<PQNObjective<Real>>(secant_,x,g);
  //Ptr<Problem<Real>> problem = makePtr<Problem<Real>>(qobj,xs);

  // Output
  if (verbosity_ > 0) writeOutput(outStream,true);

  // Compute steepest descent step
  xs->set(*state_->iterateVec);
  state_->iterateVec->set(x);
  while (status_->check(*state_)) {
    // Compute step
    qobj->setAnchor(x,*state_->gradientVec);
    gtol = std::max(sp_tol_min_,std::min(sp_tol1_,sp_tol2_*state_->gnorm));
    list_.sublist("Status Test").set("Gradient Tolerance",gtol);
    if (algoName_ == "Line Search") algo = makePtr<TypeP::ProxGradientAlgorithm<Real>>(list_);
    else if (algoName_ == "iPiano") algo = makePtr<TypeP::iPianoAlgorithm<Real>>(list_);
    else                            algo = makePtr<TypeP::SpectralGradientAlgorithm<Real>>(list_);
    algo->run(*xs,*qobj,nobj,outStream);
    s->set(*xs); s->axpy(-one,x);
    spgIter_ = algo->getState()->iter;
    state_->nprox += staticPtrCast<const TypeP::AlgorithmState<Real>>(algo->getState())->nprox;

    // Perform backtracking line search 
    state_->searchSize = one;
    x.set(*state_->iterateVec);
    x.axpy(state_->searchSize,*s);
    sobj.update(x,UpdateType::Trial);
    nobj.update(x,UpdateType::Trial);
    strial = sobj.value(x,tol);
    ntrial = nobj.value(x,tol);
    ftrial = strial + ntrial;
    ls_nfval_ = 1;
    gs = state_->gradientVec->apply(*s);
    Qk = gs + ntrial - state_->nvalue;
    if (verbosity_ > 1) {
      outStream << "  In TypeP::QuasiNewtonAlgorithm: Line Search"  << std::endl;
      outStream << "    Step size:                        " << state_->searchSize   << std::endl;
      outStream << "    Trial objective value:            " << ftrial               << std::endl;
      outStream << "    Computed reduction:               " << state_->value-ftrial << std::endl;
      outStream << "    Dot product of gradient and step: " << gs                   << std::endl;
      outStream << "    Sufficient decrease bound:        " << -Qk*c1_              << std::endl;
      outStream << "    Number of function evaluations:   " << ls_nfval_            << std::endl;
    }
    while ( ftrial > state_->value + c1_*Qk && ls_nfval_ < maxit_ ) {
      rhoTmp = -half * Qk / (strial-state_->svalue-state_->searchSize*gs);
      state_->searchSize = ((sigma1_ <= rhoTmp && rhoTmp <= sigma2_) ? rhoTmp : rhodec_) * state_->searchSize;
      x.set(*state_->iterateVec);
      x.axpy(state_->searchSize,*s);
      sobj.update(x,UpdateType::Trial);
      nobj.update(x,UpdateType::Trial);
      strial = sobj.value(x,tol);
      ntrial = nobj.value(x,tol);
      ftrial = strial + ntrial;
      Qk     = state_->searchSize * gs + ntrial - state_->nvalue;
      ls_nfval_++;
      if (verbosity_ > 1) {
        outStream << std::endl;
        outStream << "    Step size:                        " << state_->searchSize   << std::endl;
        outStream << "    Trial objective value:            " << ftrial               << std::endl;
        outStream << "    Computed reduction:               " << state_->value-ftrial << std::endl;
        outStream << "    Dot product of gradient and step: " << gs                   << std::endl;
        outStream << "    Sufficient decrease bound:        " << -Qk*c1_              << std::endl;
        outStream << "    Number of function evaluations:   " << ls_nfval_            << std::endl;
      }
    }
    state_->nsval += ls_nfval_;
    state_->nnval += ls_nfval_;

    // Compute norm of step
    state_->stepVec->set(*s);
    state_->stepVec->scale(state_->searchSize);
    state_->snorm = state_->stepVec->norm();

    // Update iterate
    state_->iterateVec->set(x);

    // Compute new value and gradient
    state_->iter++;
    state_->value  = ftrial;
    state_->svalue = strial;
    state_->nvalue = ntrial;
    sobj.update(x,UpdateType::Accept,state_->iter);
    nobj.update(x,UpdateType::Accept,state_->iter);
    gold->set(*state_->gradientVec);
    sobj.gradient(*state_->gradientVec,x,tol); state_->ngrad++;
    gp->set(state_->gradientVec->dual());

    // Compute projected gradient norm
    pgstep(*xs,*s,nobj,x,*gp,t0_,tol);
    state_->gnorm = s->norm() / t0_;

    // Update secant
    secant_->updateStorage(x,*state_->gradientVec,*gold,*state_->stepVec,state_->snorm,state_->iter);

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,writeHeader_);
  }
  if (verbosity_ > 0) TypeP::Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void QuasiNewtonAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(114,'-') << std::endl;
    hist << "Line-Search Proximal Quasi-Newton with " << secantName_ << " Hessian approximation";
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  value    - Objective function value" << std::endl;
    hist << "  gnorm    - Norm of the gradient" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  alpha    - Line search step length" << std::endl;
    hist << "  #sval    - Cumulative number of times the smooth objective function was evaluated" << std::endl;
    hist << "  #nval    - Cumulative number of times the nonsmooth objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << "  #prox    - Cumulative number of times the projection was computed" << std::endl;
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
  hist << std::setw(10) << std::left << "#sval";
  hist << std::setw(10) << std::left << "#nval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::setw(10) << std::left << "#prox";
  hist << std::setw(10) << std::left << "#ls_fval";
  hist << std::setw(10) << std::left << "sp_iter";
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void QuasiNewtonAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << "Line-Search Proximal Quasi-Newton (Type P)" << std::endl;
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
    hist << std::setw(10) << std::left << state_->nsval;
    hist << std::setw(10) << std::left << state_->nnval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << state_->nprox;
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
    hist << std::setw(10) << std::left << state_->nsval;
    hist << std::setw(10) << std::left << state_->nnval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << state_->nprox;
    hist << std::setw(10) << std::left << ls_nfval_;
    hist << std::setw(10) << std::left << spgIter_;
    hist << std::endl;
  }
  os << hist.str();
}

} // namespace TypeP
} // namespace ROL

#endif
