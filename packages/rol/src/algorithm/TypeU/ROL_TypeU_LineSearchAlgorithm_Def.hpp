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

#ifndef ROL_TYPEU_LINESEARCHALGORITHM_DEF_H
#define ROL_TYPEU_LINESEARCHALGORITHM_DEF_H

#include "ROL_LineSearch_U_Factory.hpp"
#include "ROL_DescentDirection_U_Factory.hpp"

namespace ROL {
namespace TypeU {

template<typename Real>
LineSearchAlgorithm<Real>::LineSearchAlgorithm( ParameterList &parlist,
                              const Ptr<DescentDirection_U<Real>> &descent,
                              const Ptr<LineSearch_U<Real>> &lineSearch )
  : Algorithm<Real>(), desc_(descent), lineSearch_(lineSearch),
    edesc_(DESCENT_U_USERDEFINED), els_(LINESEARCH_U_USERDEFINED),
    econd_(CURVATURECONDITION_U_WOLFE), verbosity_(0) {
  // Set status test
  status_->reset();
  status_->add(makePtr<StatusTest<Real>>(parlist));

  // Parse parameter list
  ParameterList& Llist = parlist.sublist("Step").sublist("Line Search");
  ParameterList& Glist = parlist.sublist("General");
  econd_ = StringToECurvatureConditionU(Llist.sublist("Curvature Condition").get("Type","Strong Wolfe Conditions") );
  acceptLastAlpha_ = Llist.get("Accept Last Alpha", false); 
  verbosity_ = Glist.get("Output Level",0);
  printHeader_ = verbosity_ > 2;
  // Initialize Line Search
  if (lineSearch_ == nullPtr) {
    lineSearchName_ = Llist.sublist("Line-Search Method").get("Type","Cubic Interpolation"); 
    els_ = StringToELineSearchU(lineSearchName_);
    lineSearch_ = LineSearchUFactory<Real>(parlist);
  } 
  else { // User-defined linesearch provided
    lineSearchName_ = Llist.sublist("Line-Search Method").get("User Defined Line Search Name","Unspecified User Defined Line Search");
  }
  if (desc_ == nullPtr) {
    ParameterList& dlist = Llist.sublist("Descent Method");
    descentName_ = dlist.get("Type","Quasi-Newton Method");
    edesc_ = StringToEDescentU(descentName_);
    desc_  = DescentDirectionUFactory<Real>(parlist);
  }
  else {
    descentName_ = Llist.sublist("Descent Method").get("User Defined Descent Direction Name","Unspecified User Defined Descent Direction");
  }
}

template<typename Real>
void LineSearchAlgorithm<Real>::initialize(const Vector<Real> &x,
                                           const Vector<Real> &g,
                                           Objective<Real>    &obj,
                                           std::ostream &outStream) {
  // Initialize data
  Algorithm<Real>::initialize(x,g);
  lineSearch_->initialize(x,g);
  desc_->initialize(x,g);
  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  obj.update(x,UpdateType::Initial,state_->iter);    
  state_->value = obj.value(x,ftol); 
  state_->nfval++;
  obj.gradient(*state_->gradientVec,x,ftol);
  state_->ngrad++;
  state_->gnorm = state_->gradientVec->norm();
  state_->snorm = ROL_INF<Real>();
}

template<typename Real>
void LineSearchAlgorithm<Real>::run( Vector<Real>       &x,
                                     const Vector<Real> &g,
                                     Objective<Real>    &obj,
                                     std::ostream       &outStream ) {
  const Real one(1);
  // Initialize trust-region data
  Real ftrial(0), gs(0), tol(std::sqrt(ROL_EPSILON<Real>()));
  initialize(x,g,obj,outStream);
  state_->searchSize = one;
  Ptr<Vector<Real>> gprev = g.clone();

  // Output
  if (verbosity_ > 0) writeOutput(outStream, true);

  while (status_->check(*state_)) {
    // Compute descent direction
    desc_->compute(*state_->stepVec,state_->snorm,gs,SPiter_,SPflag_,x,
                   *state_->gradientVec,obj);

    // Perform line search
    ftrial = state_->value;
    ls_nfval_ = 0; ls_ngrad_ = 0;
    lineSearch_->run(state_->searchSize,ftrial,ls_nfval_,ls_ngrad_,gs,*state_->stepVec,x,obj);

    // Make correction if maximum function evaluations reached
    if(!acceptLastAlpha_) {
      lineSearch_->setMaxitUpdate(state_->searchSize,ftrial,state_->value);
    }

    // Compute scaled descent direction
    state_->stepVec->scale(state_->searchSize);
    state_->snorm *= state_->searchSize;

    // Update iterate
    x.plus(*state_->stepVec);

    // Compute new value and gradient
    obj.update(x,UpdateType::Accept,state_->iter);
    state_->value = obj.value(x,tol);
    state_->nfval++;
    gprev->set(*state_->gradientVec);
    obj.gradient(*state_->gradientVec,x,tol);
    state_->ngrad++;
    state_->gnorm = state_->gradientVec->norm();

    // Update algorithm state
    state_->iterateVec->set(x);
    state_->nfval += ls_nfval_;
    state_->ngrad += ls_ngrad_;
    desc_->update(x,*state_->stepVec,*gprev,*state_->gradientVec,state_->snorm,state_->iter);
    state_->iter++;

    // Update Output
    if (verbosity_ > 0) writeOutput(outStream,printHeader_);
  }
  if (verbosity_ > 0) Algorithm<Real>::writeExitStatus(outStream);
}

template<typename Real>
void LineSearchAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  std::stringstream hist;
  if (verbosity_ > 1) {
    hist << std::string(109,'-') << std::endl;
    hist << descentName_;
    hist << " status output definitions" << std::endl << std::endl;
    hist << "  iter     - Number of iterates (steps taken)" << std::endl;
    hist << "  value    - Objective function value" << std::endl;
    hist << "  gnorm    - Norm of the gradient" << std::endl;
    hist << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    hist << "  alpha    - Line search step length" << std::endl;
    hist << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    hist << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    hist << "  ls_#fval - Number of the times the objective function was evaluated during the line search" << std::endl;
    hist << "  ls_#grad - Number of the times the gradient was evaluated during the line search" << std::endl;
    if (edesc_ == DESCENT_U_NEWTONKRYLOV) {
      hist << "  iterCG   - Number of Krylov iterations used to compute search direction" << std::endl;
      hist << "  flagCG   - Krylov solver flag" << std::endl;
    }
    hist << std::string(109,'-') << std::endl;
  }

  hist << "  ";
  hist << std::setw(6)  << std::left << "iter";
  hist << std::setw(15) << std::left << "value";
  hist << std::setw(15) << std::left << "gnorm";
  hist << std::setw(15) << std::left << "snorm";
  hist << std::setw(15) << std::left << "alpha";
  hist << std::setw(10) << std::left << "#fval";
  hist << std::setw(10) << std::left << "#grad";
  hist << std::setw(10) << std::left << "ls_#fval";
  hist << std::setw(10) << std::left << "ls_#grad";
  if (edesc_ == DESCENT_U_NEWTONKRYLOV) {
    hist << std::setw(10) << std::left << "iterCG";
    hist << std::setw(10) << std::left << "flagCG";
  }
  hist << std::endl;
  os << hist.str();
}

template<typename Real>
void LineSearchAlgorithm<Real>::writeName( std::ostream& os ) const {
  std::stringstream hist;
  hist << std::endl << desc_->printName();
  hist << std::endl;
  hist << "Line Search: " << lineSearchName_;
  hist << " satisfying " << ECurvatureConditionUToString(econd_) << std::endl;
  os << hist.str();
}

template<typename Real>
void LineSearchAlgorithm<Real>::writeOutput( std::ostream& os, bool print_header) const {
  std::stringstream hist;
  hist << std::scientific << std::setprecision(6);
  if ( state_->iter == 0 ) {
    writeName(os);
  }
  if ( print_header ) {
    writeHeader(os);
  }
  if ( state_->iter == 0 ) {
    hist << "  ";
    hist << std::setw(6)  << std::left << state_->iter;
    hist << std::setw(15) << std::left << state_->value;
    hist << std::setw(15) << std::left << state_->gnorm;
    hist << std::setw(15) << std::left << "---";
    hist << std::setw(15) << std::left << "---";
    hist << std::setw(10) << std::left << state_->nfval;
    hist << std::setw(10) << std::left << state_->ngrad;
    hist << std::setw(10) << std::left << "---";
    hist << std::setw(10) << std::left << "---";
    if (edesc_ == DESCENT_U_NEWTONKRYLOV) {
      hist << std::setw(10) << std::left << "---";
      hist << std::setw(10) << std::left << "---";
    }
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
    hist << std::setw(10) << std::left << ls_nfval_;
    hist << std::setw(10) << std::left << ls_ngrad_;
    if (edesc_ == DESCENT_U_NEWTONKRYLOV) {
      hist << std::setw(10) << std::left << SPiter_;
      hist << std::setw(10) << std::left << SPflag_;
    }
    hist << std::endl;
  }
  os << hist.str();
}
} // namespace TypeU
} // namespace ROL

#endif
