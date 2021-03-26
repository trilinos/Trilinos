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

#ifndef ROL2_TYPEU_LINESEARCHALGORITHM_DEF_HPP
#define ROL2_TYPEU_LINESEARCHALGORITHM_DEF_HPP

namespace ROL2 {
namespace TypeU {

template<typename Real>
LineSearchAlgorithm<Real>::LineSearchAlgorithm(       ParameterList&               parlist,
                                                const Ptr<DescentDirection<Real>>& desc,
                                                const Ptr<LineSearch<Real>>&       ls )
  : Algorithm<Real>(), desc_(desc), lineSearch_(ls), verbosity_(0) {

  // Set status test
  auto& status = Algorithm<Real>::getStatus();
  status.reset();
  status.add(makePtr<StatusTest<Real>>(parlist));

  // Parse parameter list
  auto& Llist  = parlist.sublist("Step").sublist("Line Search");
  auto& LMlist = Llist.sublist("Line-Search Method");
  auto& Glist  = parlist.sublist("General");
  auto& CClist = Llist.sublist("Curvature Condition");

  econd_ = stringToEnum(CClist.get("Type","Strong Wolfe Conditions"), CurvatureCond{} );
  acceptLastAlpha_ = Llist.get("Accept Last Alpha", false); 
  verbosity_ = Glist.get("Output Level",0);
  printHeader_ = verbosity_ > 2;

  // Initialize Line Search
  if( lineSearch_ == nullPtr) {
    lineSearchName_ = LMlist.get("Type","Cubic Interpolation"); 
    els_ = stringToEnum(lineSearchName_, Type{});
    lineSearch_ = LineSearchUFactory<Real>(parlist); // TODO
  } 
  else { // User-defined linesearch provided
    lineSearchName_ = LMlist.get("User Defined Line Search Name","Unspecified User Defined Line Search");
  }
  if (desc_ == nullPtr) {
    auto& dlist = Llist.sublist("Descent Method");
    descentName_ = dlist.get("Type","Quasi-Newton Method");
    edesc_ = stringToEnum(descentName_,DescentType{});
    desc_  = DescentDirectionUFactory<Real>(parlist); // TODO
  }
  else {
    descentName_ = LMlist.get("User Defined Descent Direction Name","Unspecified User Defined Descent Direction");
  }
}

template<typename Real>
void LineSearchAlgorithm<Real>::initialize( const Vector<Real>&    x,
                                            const Vector<Real>&    g,
                                                  Objective<Real>& obj,
                                                  std::ostream&    outStream ) {
  auto& state = Algorithm<Real>::getState();

  // Initialize data
  Algorithm<Real>::initialize(x,g);
  lineSearch_->initialize(x,g);
  desc_->initialize(x,g);

  // Update approximate gradient and approximate objective function.
  Real ftol = std::sqrt(ROL_EPSILON<Real>());
  obj.update(x,UPDATE_INITIAL,state.iter);    
  state.value = obj.value(x,ftol); 
  state.nfval++;
  obj.gradient(state.gradientVec,x,ftol);
  state.ngrad++;
  state.gnorm = state.gradientVec->norm();
  state.snorm = ROL_INF<Real>;
}

template<typename Real>
void LineSearchAlgorithm<Real>::run(       Vector<Real>&    x,
                                     const Vector<Real>&    g,
                                           Objective<Real>& obj,
                                           std::ostream&    outStream ) {
  auto& state = Algorithm<Real>::getState();
  const Real one(1);

  // Initialize trust-region data
  Real ftrial(0), gs(0), tol(std::sqrt(ROL_EPSILON<Real>));
  initialize(x,g,obj,outStream);
  state.searchSize = one;
  auto gprev = g.clone();

  // Output
  output.push_back(print(true));
  if (verbosity_ > 0) outStream << print(true);

  while( status.check(state) ) {

    // Compute descent direction
    desc_->compute(state.stepVec,state.snorm,gs,SPiter_,SPflag_,x,
                   state.gradientVec,obj);

    // Perform line search
    ftrial = state.value;
    ls_nfval_ = 0; 
    ls_ngrad_ = 0;
    lineSearch_->run(state.searchSize,ftrial,ls_nfval_,ls_ngrad_,gs,state.stepVec,x,obj);

    // Make correction if maximum function evaluations reached
    if(!acceptLastAlpha_) {
      lineSearch_->setMaxitUpdate(state.searchSize,ftrial,state.value);
    }

    // Compute scaled descent direction
    state.stepVec->scale(state.searchSize);
    state.snorm *= state.searchSize;

    // Update iterate
    x.plus(state.stepVec);

    // Compute new value and gradient
    obj.update(x,UPDATE_ACCEPT,state.iter);
    state.value = obj.value(x,tol);
    state.nfval++;
    gprev->set(state.gradientVec);
    obj.gradient(state.gradientVec,x,tol);
    state.ngrad++;
    state.gnorm = state.gradientVec->norm();

    // Update algorithm state
    state.iterateVec->set(x);
    state.nfval += ls_nfval_;
    state.ngrad += ls_ngrad_;
    desc_->update(x,state.stepVec,*gprev,state.gradientVec,state.snorm,state.iter);
    state.iter++;

    // Update Output
    if (verbosity_ > 0) writeHeader(outStream);
  }
  if (verbosity_ > 0) writeOutput( Algorithm<Real>::printExitStatus() );
}

template<typename Real>
void LineSearchAlgorithm<Real>::writeHeader( std::ostream& os ) const {
  if( verbosity_ > 1 ) {
    os << std::string(109,'-') << std::endl;
    os << descentName_;
    os << " status output definitions" << std::endl << std::endl;
    os << "  iter     - Number of iterates (steps taken)" << std::endl;
    os << "  value    - Objective function value" << std::endl;
    os << "  gnorm    - Norm of the gradient" << std::endl;
    os << "  snorm    - Norm of the step (update to optimization vector)" << std::endl;
    os << "  alpha    - Line search step length" << std::endl;
    os << "  #fval    - Cumulative number of times the objective function was evaluated" << std::endl;
    os << "  #grad    - Cumulative number of times the gradient was computed" << std::endl;
    os << "  ls_#fval - Number of the times the objective function was evaluated during the line search" << std::endl;
    os << "  ls_#grad - Number of the times the gradient was evaluated during the line search" << std::endl;
    if( edesc_ == DescentType::NewtonKrylov ) {
      os << "  iterCG   - Number of Krylov iterations used to compute search direction" << std::endl;
      os << "  flagCG   - Krylov solver flag" << std::endl;
    }
    os << std::string(109,'-') << std::endl;
  }

  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "gnorm";
  os << std::setw(15) << std::left << "snorm";
  os << std::setw(15) << std::left << "alpha";
  os << std::setw(10) << std::left << "#fval";
  os << std::setw(10) << std::left << "#grad";
  os << std::setw(10) << std::left << "ls_#fval";
  os << std::setw(10) << std::left << "ls_#grad";
  if (edesc_ == DESCENT_NEWTONKRYLOV) {
    os << std::setw(10) << std::left << "iterCG";
    os << std::setw(10) << std::left << "flagCG";
  }
  os << std::endl;
}

template<typename Real>
void LineSearchAlgorithm<Real>::writeName( std::ostream& os ) const {
  desc_->writeName(os);
  os << "\n\n";
  os << "Line Search: " << lineSearchName_;
  os << " satisfying " << typeToString(econd_) << std::endl;
}

template<typename Real>
std::string LineSearchAlgorithm<Real>::print( print_header ) const {
  os << std::scientific << std::setprecision(6);
  if( state.iter == 0 ) writeName(os);
  if( print_header )    writeHeader(os);
  if( state.iter == 0 ) {
    os << "  ";
    os << std::setw(6)  << std::left << state.iter;
    os << std::setw(15) << std::left << state.value;
    os << std::setw(15) << std::left << state.gnorm;
    os << std::setw(15) << std::left << "---";
    os << std::setw(15) << std::left << "---";
    os << std::setw(10) << std::left << state.nfval;
    os << std::setw(10) << std::left << state.ngrad;
    os << std::setw(10) << std::left << "---";
    os << std::setw(10) << std::left << "---";
    if( edesc_ == DescentType::Newton) {
      os << std::setw(10) << std::left << "---";
      os << std::setw(10) << std::left << "---";
    }
    os << std::endl;
  }
  else {
    os << "  ";
    os << std::setw(6)  << std::left << state.iter;
    os << std::setw(15) << std::left << state.value;
    os << std::setw(15) << std::left << state.gnorm;
    os << std::setw(15) << std::left << state.snorm;
    os << std::setw(15) << std::left << state.searchSize;
    os << std::setw(10) << std::left << state.nfval;
    os << std::setw(10) << std::left << state.ngrad;
    os << std::setw(10) << std::left << ls_nfval_;
    os << std::setw(10) << std::left << ls_ngrad_;
    if(edesc_ == DESCENT_NEWTONKRYLOV) {
      os << std::setw(10) << std::left << SPiter_;
      os << std::setw(10) << std::left << SPflag_;
    }
    os << std::endl;
  }
  return os.str();
}

} // namespace TypeU
} // namespace ROL

#endif // ROL2_TYPEU_LINESEARCHALGORITHM_DEF_HPP
