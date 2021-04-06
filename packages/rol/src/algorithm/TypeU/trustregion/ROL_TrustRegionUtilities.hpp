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

#ifndef ROL_TRUSTREGIONUTILITIES_U_H
#define ROL_TRUSTREGIONUTILITIES_U_H

#include "ROL_TrustRegionModel_U.hpp"

namespace ROL {
namespace TRUtils {

/** \enum  ROL::TRUtils::ETRFlag 
    \brief Enumation of flags used by trust-region solvers.

    \arg SUCCESS        Actual and predicted reductions are positive 
    \arg POSPREDNEG     Reduction is positive, predicted negative (impossible)
    \arg NPOSPREDPOS    Reduction is nonpositive, predicted positive
    \arg NPOSPREDNEG    Reduction is nonpositive, predicted negative (impossible)
    \arg TRNAN          Actual and/or predicted reduction is NaN

*/
enum ETRFlag {
  SUCCESS = 0,
  POSPREDNEG,
  NPOSPREDPOS,
  NPOSPREDNEG,
  TRNAN,
  QMINSUFDEC,
  UNDEFINED 
};

inline std::string ETRFlagToString(ETRFlag trf) {
  std::string retString;
  switch(trf) {
    case SUCCESS:  
      retString = "Both actual and predicted reductions are positive (success)";
      break;
    case POSPREDNEG: 
      retString = "Actual reduction is positive and predicted reduction is negative (impossible)";
      break;
    case NPOSPREDPOS: 
      retString = "Actual reduction is nonpositive and predicted reduction is positive";
      break;
    case NPOSPREDNEG:
      retString = "Actual reduction is nonpositive and predicted reduction is negative (impossible)";
      break;
    case TRNAN:
      retString = "Actual and/or predicted reduction is a NaN";
      break;
    case QMINSUFDEC:
      retString = "Subproblem solution did not produce sufficient decrease";
      break;
    default:
      retString = "INVALID ETRFlag";       
  }
  return retString;
}

template<typename Real>
inline Real initialRadius(int &nfval,
                          const Vector<Real> &x,
                          const Vector<Real> &g,
                          Vector<Real> &Bg,
                          const Real fx,
                          const Real gnorm,
                          Objective<Real> &obj,
                          TrustRegionModel_U<Real> &model,
                          const Real delMax,
                          std::ostream &outStream,
                          const bool print = false) {
  const Real zero(0), half(0.5), one(1), two(2), three(3), six(6);
  const Real eps(ROL_EPSILON<Real>());
  Real del(ROL_INF<Real>());
  Ptr<Vector<Real>> xcp = x.clone();
  model.setData(obj,x,g);
  Real htol = std::sqrt(eps);
  model.hessVec(Bg,g.dual(),x,htol);
  Real gBg = Bg.dot(g);
  Real alpha = one;
  if ( gBg > eps ) {
    alpha = gnorm*gnorm/gBg;
  }
  // Evaluate the objective function at the Cauchy point
  xcp->set(g.dual());
  xcp->scale(-alpha);
  //Real gs = xcp->dot(g.dual());
  Real gs = xcp->apply(g);
  xcp->plus(x);
  obj.update(*xcp,UpdateType::Temp);
  Real ftol = static_cast<Real>(0.1)*ROL_OVERFLOW<Real>(); 
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
    if ( b*b-three*a*c > eps ) {
      // There is at least one critical point
      Real t1 = (-b-std::sqrt(b*b-three*a*c))/(three*a);
      Real t2 = (-b+std::sqrt(b*b-three*a*c))/(three*a);
      if ( six*a*t1 + two*b > zero ) {
        // t1 is the minimizer
        del = std::min(t1*alpha*gnorm,delMax);
      }
      else {
        // t2 is the minimizer
        del = std::min(t2*alpha*gnorm,delMax);
      }
    }
    else {
      del = std::min(alpha*gnorm,delMax);
    }
  }
  if (del <= eps*gnorm) {
    del = one;
  }
  obj.update(x,UpdateType::Revert);
  if ( print ) {
    outStream << "  In TrustRegionUtilities::initialRadius"      << std::endl;
    outStream << "    Initial radius:                          " << del << std::endl;
  }
  return del;
}

template<typename Real>
inline void analyzeRatio(Real &rho,
                         ETRFlag &flag,
                   const Real fold,
                   const Real ftrial,
                   const Real pRed,
                   const Real epsi,
                         std::ostream &outStream = std::cout,
                   const bool print = false) {
  const Real zero(0), one(1);
  Real eps       = epsi*std::max(one,fold);
  Real aRed      = fold - ftrial;
  Real aRed_safe = aRed + eps, pRed_safe = pRed + eps;
  if (((std::abs(aRed_safe) < epsi) && (std::abs(pRed_safe) < epsi)) || aRed == pRed) {
    rho  = one;
    flag = SUCCESS;
  }
  else if ( std::isnan(aRed_safe) || std::isnan(pRed_safe) ) {
    rho  = -one;
    flag = TRNAN;
  }
  else {
    rho = aRed_safe/pRed_safe;
    if (pRed_safe < zero && aRed_safe > zero) {
      flag = POSPREDNEG;
    }
    else if (aRed_safe <= zero && pRed_safe > zero) {
      flag = NPOSPREDPOS;
    }
    else if (aRed_safe <= zero && pRed_safe < zero) {
      flag = NPOSPREDNEG;
    }
    else {
      flag = SUCCESS;
    }
  }
  if ( print ) {
    outStream << "  In TrustRegionUtilities::analyzeRatio"       << std::endl;
    outStream << "    Current objective function value:        " << fold      << std::endl;
    outStream << "    New objective function value:            " << ftrial    << std::endl;
    outStream << "    Actual reduction:                        " << aRed      << std::endl;
    outStream << "    Predicted reduction:                     " << pRed      << std::endl;
    outStream << "    Safeguard:                               " << epsi      << std::endl;
    outStream << "    Actual reduction with safeguard:         " << aRed_safe << std::endl;
    outStream << "    Predicted reduction with safeguard:      " << pRed_safe << std::endl;
    outStream << "    Ratio of actual and predicted reduction: " << rho       << std::endl;
    outStream << "    Trust-region flag:                       " << flag      << std::endl;
    outStream << std::endl;
  }
}

template<typename Real>
inline Real interpolateRadius(const Vector<Real> &g,
                              const Vector<Real> &s,
                              const Real snorm,
                              const Real pRed,
                              const Real fold,
                              const Real ftrial,
                              const Real del,
                              const Real gamma0,
                              const Real gamma1,
                              const Real eta2,
                              std::ostream &outStream = std::cout,
                              const bool print = false) {
  const Real one(1);
  //Real gs = g.dot(s.dual());
  Real gs = g.apply(s);
  Real modelVal = fold - pRed;
  Real theta = (one-eta2)*gs/((one-eta2)*(fold+gs)+eta2*modelVal-ftrial);
  if ( print ) {
    outStream << "  In TrustRegionUtilities::interpolateRadius"  << std::endl;
    outStream << "    Interpolation model value:               " << modelVal << std::endl;
    outStream << "    Interpolation step length:               " << theta    << std::endl;
    outStream << std::endl;
  }
  return std::min(gamma1*std::min(snorm,del),std::max(gamma0,theta)*del);
}

} // namespace TrustRegion
} // namespace ROL


#endif
