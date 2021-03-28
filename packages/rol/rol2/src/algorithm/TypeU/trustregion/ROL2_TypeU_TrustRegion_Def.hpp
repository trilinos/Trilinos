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
#ifndef ROL2_TYPEU_TRUSTREGION_DEF_H
#define ROL2_TYPEU_TRUSTREGION_DEF_H

namespace ROL2 {
namespace TypeU {


template<class Real>
Real TrustRegion<Real>::initialRadius(       int&                    nfval,
                                       const Vector<Real>&           x,
                                       const Vector<Real>&           g,
                                             Vector<Real>&           Bg,
                                             Real                    fx,
                                             Real                    gnorm,
                                             Objective<Real>&        obj,
                                             TrustRegionModel<Real>& model,
                                             Real                    delMax,
                                             std::ostream&           os,
                                             bool                    print ) {
  const Real zero{0}, half{0.5}, one{1}, two{2}, three{3}, six{6};
  const Real eps{ ROL_EPSILON<Real> };
  Real del{ ROL_INF<Real> };
  auto xcp = x.clone();
  model.setData(obj,x,g);
  Real htol = std::sqrt(eps);
  model.hessVec(Bg,g.dual(),x,htol);
  Real gBg = Bg.dot(g);
  Real alpha = one;
  if ( gBg > eps )  alpha = gnorm*gnorm/gBg;
   
  // Evaluate the objective function at the Cauchy point
  xcp->set(g.dual());
  xcp->scale(-alpha);

  //Real gs = xcp->dot(g.dual());
  Real gs = xcp->apply(g);
  xcp->plus(x);
  obj.update(*xcp,UpdateType::Temp);
  Real ftol = static_cast<Real>(0.1)*ROL_OVERFLOW<Real>; 
  Real fnew = obj.value(*xcp,ftol); // MUST DO SOMETHING HERE WITH FTOL
  nfval++;

  // Perform cubic interpolation to determine initial trust region radius
  Real a = fnew - fx - gs - half*alpha*alpha*gBg;
  if ( std::abs(a) < eps ) { 
    // a = 0 implies the objective is quadratic in the negative gradient direction
    del = std::min(alpha*gnorm,delMax);
  }
  else {
    Real b = half*alpha*alpha*gBg;
    Real c = gs;
    Real discriminant = b*b-three*a*c;
    if ( discriminant > eps ) {
      // There is at least one critical point
      Real t1 = (-b-std::sqrt(discriminant))/(three*a);
      Real t2 = (-b+std::sqrt(discriminant))/(three*a);
      if ( six*a*t1 + two*b > zero ) {
        // t1 is the minimizer
        del = std::min(t1*alpha*gnorm,delMax);
      }
      else {
        // t2 is the minimizer
        del = std::min(t2*alpha*gnorm,delMax);
      }
    }
    else del = std::min(alpha*gnorm,delMax);
  }
  if (del <= eps*gnorm)  del = one;

  obj.update(x,UpdateType::Revert);
  if ( print ) {
    os << "  In TrustRegionUtilities::initialRadius"      << std::endl;
    os << "    Initial radius:                          " << del << std::endl;
  }
  return del;
}

template<class Real>
void TrustRegion<Real>::analyzeRatio( Real&                    rho,
                                      TrustRegion<Real>::Flag& flag,
                                      Real                     fold,
                                      Real                     ftrial,
                                      Real                     pRed,
                                      Real                     epsi,
                                      std::ostream&            os,
                                      bool                     print ) {
  const Real zero{0}, one{1};
  Real eps       = epsi*std::max(one,fold);
  Real aRed      = fold - ftrial;
  Real aRed_safe = aRed + eps, pRed_safe = pRed + eps;
  if (((std::abs(aRed_safe) < epsi) && (std::abs(pRed_safe) < epsi)) || aRed == pRed) {
    rho  = one;
    flag = Flag::Success; 
  }
  else if ( std::isnan(aRed_safe) || std::isnan(pRed_safe) ) {
    rho  = -one;
    flag = Flag::TRNaN;
  }
  else {
    rho = aRed_safe/pRed_safe;
    if (pRed_safe < zero && aRed_safe > zero)       flag = Flag::PosPredNeg;
    else if (aRed_safe <= zero && pRed_safe > zero) flag = Flag::NPosPredPos;
    else if (aRed_safe <= zero && pRed_safe < zero) flag = Flag::NPosPredNeg;
    else                                            flag = Flag::Success;
  }

  if ( print ) {
    os << "  In TrustRegionUtilities::analyzeRatio"       << std::endl;
    os << "    Current objective function value:        " << fold      << std::endl;
    os << "    New objective function value:            " << ftrial    << std::endl;
    os << "    Actual reduction:                        " << aRed      << std::endl;
    os << "    Predicted reduction:                     " << pRed      << std::endl;
    os << "    Safeguard:                               " << epsi      << std::endl;
    os << "    Actual reduction with safeguard:         " << aRed_safe << std::endl;
    os << "    Predicted reduction with safeguard:      " << pRed_safe << std::endl;
    os << "    Ratio of actual and predicted reduction: " << rho       << std::endl;
    os << "    Trust-region flag:                       " << flag      << std::endl;
    os << std::endl;
  }
}

template<class Real>
Real TrustRegion<Real>::interpolateRadius( const Vector<Real>& g,
                                           const Vector<Real>& s,
                                                 Real          snorm,
                                                 Real          pRed,
                                                 Real          fold,
                                                 Real          ftrial,
                                                 Real          del,
                                                 Real          gamma0,
                                                 Real          gamma1,
                                                 Real          eta2,
                                                 std::ostream& os,
                                                 bool          print ) {
  const Real one(1);

  Real gs = g.apply(s);
  Real modelVal = fold - pRed;
  Real theta = (one-eta2)*gs/((one-eta2)*(fold+gs)+eta2*modelVal-ftrial);

  if ( print ) {
    os << "  In TrustRegionUtilities::interpolateRadius"  << std::endl;
    os << "    Interpolation model value:               " << modelVal << std::endl;
    os << "    Interpolation step length:               " << theta    << std::endl;
    os << std::endl;
  }
  return std::min(gamma1*std::min(snorm,del),std::max(gamma0,theta)*del);
}

template<class Real>
Ptr<TrustRegion<Real>> 
TrustRegion<Real>::create( ParameterList& parlist ) {
  return nullPtr;
//  auto& trlist = parlist.sublist("Step").sublist("Trust Region");
//  auto soltype = type_dict[trlist.get("Subproblem Solver","DogLeg")];
//  switch( soltype ) {
//    case: Type::CauchyPoint:  return makePtr<CauchyPoint<Real>>();
//    case: Type::DogLeg:       return makePtr<DogLeg<Real>>();
//    case: Type::DoubleDogLeg: return makePtr<<Real>>();
//    case: Type::TruncateCG:   return makePtr<<Real>>(parlist);
//    case: Type::SPG:          return makePtr<<Real>>(parlist);
//    default: return nullPtr; // TODO: Should we not throw an exception?
//  }
} // TrustRegion<Real>::create


} // namespace TypeU
} // namespace ROL2


#endif // ROL2_TYPEU_TRUSTREGION_DEF_H
