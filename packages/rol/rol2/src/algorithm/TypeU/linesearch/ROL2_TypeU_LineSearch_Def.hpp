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

#pragma once
#ifndef ROL2_TYPEU_LINESEARCH_DEF_H
#define ROL2_TYPEU_LINESEARCH_DEF_H

/** \class ROL2::TypeU::LineSearch
    \brief Provides interface for and implements line searches.
*/

namespace ROL2 {
namespace TypeU {

template<class Real>
LineSearch<Real>::LineSearch( ParameterList& parlist ) {
  Real one(1), p9(0.9), p6(0.6), p4(0.4), oem4(1.e-4), zero(0);
  auto& ls = parlist.sublist("Step").sublist("Line Search");
  auto& dm = ls.sublist("Descent Method");
  auto& cc = ls.sublist("Curvature Condition");
  edesc_ = DescentDirection<Real>::type_dict[dm.get("Type","Quasi-Newton Method")];
  econd_ = curvature_dict[cc.get("Type","Strong Wolfe Conditions")];

  // Linesearch Parameters
  alpha0_       = ls.get("Initial Step Size",one);
  alpha0bnd_    = ls.get("Lower Bound for Initial Step Size",one);
  useralpha_    = ls.get("User Defined Initial Step Size",false);
  usePrevAlpha_ = ls.get("Use Previous Step Length as Initial Guess",false);
  acceptMin_    = ls.get("Accept Linesearch Minimizer",false);
  maxit_        = ls.get("Function Evaluation Limit",20);
  c1_           = ls.get("Sufficient Decrease Tolerance",oem4);
  c2_           = cc.get("General Parameter",p9);
  c3_           = cc.get("Generalized Wolfe Parameter",p6);
 
  fmin_       = ROL_MAX<Real>;
  alphaMin_   = 0; 
  itcond_     = false;
  FDdirDeriv_ = ls.get("Finite Difference Directional Derivative",false);

  c1_ = ((c1_ < zero) ? oem4 : c1_);
  c2_ = ((c2_ < zero) ? p9   : c2_);
  c3_ = ((c3_ < zero) ? p9   : c3_);
  if ( c2_ <= c1_ ) {
    c1_ = oem4;
    c2_ = p9;
  }
  if ( edesc_ == DescentType::NonlinearCG ) {
    c2_ = p4;
    c3_ = std::min(one-c2_,c3_);
  }
}

template<class Real>
void LineSearch<Real>::initialize( const Vector<Real>& x, 
                                   const Vector<Real>& g ) {
  xtst_ = x.clone();
}

// use this function to modify alpha and fval if the maximum number of iterations
// are reached
template<class Real>
void LineSearch<Real>::setMaxitUpdate(       Real& alpha, 
                                             Real& fnew, 
                                       const Real& fold ) {
  // Use local minimizer
  if( itcond_ && acceptMin_ ) {
    alpha = alphaMin_;
    fnew = fmin_;
  }
  // Take no step
  else if(itcond_ && !acceptMin_) {
    alpha = 0;
    fnew = fold;
  }
  setNextInitialAlpha(alpha);
}

template<class Real>
bool LineSearch<Real>::status( LineSearch<Real>::Type type, 
                                     int&             ls_neval, 
                                     int&             ls_ngrad, 
                                     Real             alpha, 
                                     Real             fold, 
                                     Real             sgold, 
                                     Real             fnew, 
                               const Vector<Real>&    x, 
                               const Vector<Real>&    s, 
                                     Objective<Real>& obj ) { 
  const Real one(1), two(2);

  // Check Armijo Condition
  bool armijo =  fnew <= (fold + c1_*alpha*sgold);

  // Check Maximum Iteration
  bool itcond_ = ( ls_neval >= maxit_ );

  // Check Curvature Condition
  bool curvcond = false;
  if ( armijo && ((type != Type::BackTracking && type != Type::CubicInterp) ||
                  (edesc_ == DescentType::NonlinearCG )) ) {
    if ( econd_ == CurvatureCond::Goldstein )
      curvcond = fnew >= (fold + (one-c1_)*alpha*sgold);
    else if (econd_ == CurvatureCond::Null)
      curvcond = true;
    else {
      Real sgnew = dirDeriv(x,s,alpha,fnew,obj); //ls_ngrad++;
      if (    ((econd_ == CurvatureCond::Wolfe)
                   && (sgnew >= c2_*sgold))
           || ((econd_ == CurvatureCond::StrongWolfe)
                   && (std::abs(sgnew) <= c2_*std::abs(sgold)))
           || ((econd_ == CurvatureCond::GeneralizedWolfe)
                   && (c2_*sgold <= sgnew && sgnew <= -c3_*sgold))
           || ((econd_ == CurvatureCond::ApproximateWolfe)
                   && (c2_*sgold <= sgnew && sgnew <= (two*c1_ - one)*sgold)) ) {
        curvcond = true;
      }
    }
  }

  if( fnew < fmin_ ) {
    fmin_ = fnew;
    alphaMin_ = alpha;
  }

  if (type == Type::BackTracking || type == Type::CubicInterp) {
    if (edesc_ == DescentType::NonlinearCG)  
      return ((armijo && curvcond) || itcond_);
    else return (armijo || itcond_);
  }
  else return ((armijo && curvcond) || itcond_);
}

template<class Real>
Real LineSearch<Real>::getInitialAlpha(       int&             ls_neval, 
                                              int&             ls_ngrad, 
                                              Real             fval, 
                                              Real             gs, 
                                        const Vector<Real>&    x, 
                                        const Vector<Real>&    s, 
                                              Objective<Real>& obj) {
  Real val(1);
  
  if (useralpha_ || usePrevAlpha_ )  val = alpha0_;
  else {
    const Real one(1), half(0.5);
    if (edesc_ == DescentType::Steepest || edesc_ == DescentType::NonlinearCG ) {
      Real tol = default_tolerance<Real>();
    
      // Evaluate objective at x + s
      xtst_->set(x);
      xtst_->plus(s);
      obj.update(*xtst_,UpdateType::Trial);
      Real fnew = obj.value(*xtst_,tol);
      ls_neval++;
    
      // Minimize quadratic interpolate to compute new alpha
      Real denom = (fnew - fval - gs);
      Real alpha = ((denom > ROL_EPSILON<Real>) ? -half*gs/denom : one);
      val = ((alpha > alpha0bnd_) ? alpha : one);
    }
    else val = one;
  }
  return val;
}

template<class Real>
void LineSearch<Real>:: setNextInitialAlpha( Real alpha ) {
  if( usePrevAlpha_ ) {
    alpha0_  = alpha; 
  }
}

template<class Real>
bool LineSearch<Real>::useLocalMinimizer() {
  return itcond_ && acceptMin_;
}

template<class Real>
bool LineSearch<Real>::takeNoStep() {
  return itcond_ && !acceptMin_;
}

template<class Real>
Real LineSearch<Real>::dirDeriv( const Vector<Real>&    x,      // current iterate
                                 const Vector<Real>&    s,      // current trial step
                                       Real             alpha,  // current step length
                                       Real             fnew,   // f(x+alpha*s)
                                       Objective<Real>& obj) {
  Real tol = default_tolerance<Real>();
  Real val(0);
  xtst_->set(x); 
  xtst_->axpy(alpha,s);
  if (FDdirDeriv_) {
    Real snorm = s.norm();
    if (snorm > 0) {
      Real xnorm   = xtst_->norm();
      Real cbrteps = std::cbrt(ROL_EPSILON<Real>);
      Real h       = cbrteps*std::max(xnorm/snorm,1);
      xtst_->axpy(h,s);
      obj.update(*xtst_,UpdateType::Trial);
      Real ftrial = obj.value(*xtst_,tol);
      val = (ftrial - fnew) / h;
    }
  }
  else val = obj.dirDeriv(*xtst_,s,tol);
  return val;
}


template<class Real>
Ptr<LineSearch<Real>> 
LineSearch<Real>::create( ParameterList& parlist ) {
  auto lslist = parlist.sublist("Step").sublist("Line Search");
  auto lsm = lslist.sublist("Line-Search Method").get("Type","Cubic Interpolation");
  auto lstype = type_dict[lsm];
  switch( lstype ) {
    case Type::IterationScaling:     return makePtr<IterationScaling<Real>>(parlist);
    case Type::PathBasedTargetLevel: return makePtr<PathBasedTargetLevel<Real>>(parlist);
    case Type::Backtracking:         return makePtr<BackTracking<Real>>(parlist);
    case Type::CubicInterp:          return makePtr<CubicInterp<Real>>(parlist);
    case Type::Bisection:            
    case Type::GoldenSection:
    case Type::Brents:              // return makePtr<ScalarMinimizationLineSearch<Real>>(parlist); 
    default:  return nullPtr; // TODO should we throw an exception here?
  };
} // LineSearch<Real>::create


} // namespace TypeU
} // namespace ROL

#endif // ROL2_TYPEU_LINESEARCH_DEF_H
