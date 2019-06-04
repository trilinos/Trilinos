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

#ifndef ROL_CHEBYSHEVKUSUOKA_HPP
#define ROL_CHEBYSHEVKUSUOKA_HPP

#include "ROL_SpectralRisk.hpp"
#include "ROL_GaussChebyshev1Quadrature.hpp"
#include "ROL_GaussChebyshev2Quadrature.hpp"
#include "ROL_GaussChebyshev3Quadrature.hpp"

/** @ingroup risk_group
    \class ROL::ChebyshevSpectral
    \brief Provides an interface for the Chebyshev-Spectral risk measure.

    The Chebyshev-Spectral risk measure is defined as
    \f[
       \mathcal{R}(X) = \int_{\alpha_0}^{\alpha_1} w(\alpha)
          \mathrm{CVaR}_{\alpha}(X) \,\mathrm{d}\alpha
    \f]
    where \f$0\le \alpha_0 < \alpha_1 < 1\f$ and the conditional value-at-risk
    (CVaR) with confidence level \f$0\le \alpha < 1\f$ is
    \f[
       \mathrm{CVaR}_\alpha(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \frac{1}{1-\alpha} \mathbb{E}\left[(X-t)_+\right]
         \right\}, \quad (x)_+ = \max\{0,x\}.
    \f]
    There are three choices of weight functions \f$w\f$: (i) the first weight
    function generates the Chebyshev polynomials of the first kind and has
    the specific form
    \f[
       w(x) = \frac{1}{\sqrt{(x-\alpha_0)(\alpha_1-x)}};
    \f]
    (ii) the second weight function generates the Chebyshev polynomials of the
    second kind and has the specific form
    \f[
       w(x) = \sqrt{(x-\alpha_0)(\alpha_1-x)};
    \f]
    and (iii) the third weight function is related again to the Chebyshev
    polynomials of the first kind and has the specific form
    \f[
       w(x) = \sqrt{\frac{x-\alpha_0}{\alpha_1-x}}.
    \f]
    As defined, \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.

    ROL implements \f$\mathcal{R}\f$ by approximating the integral with
    the appropriate Gauss-Chebyshev quadrature rule.  The corresponding
    quadrature points and weights are then used to construct a
    ROL::MixedQuantileQuadrangle risk measure.
    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class ChebyshevSpectral : public SpectralRisk<Real> {
private:
  ROL::Ptr<PlusFunction<Real> > plusFunction_;

  Real lower_, upper_;
  int nQuad_; 
  int wType_;

  std::vector<Real> wts_;
  std::vector<Real> pts_;

  void checkInputs(void) const {
    ROL_TEST_FOR_EXCEPTION(lower_ > upper_, std::invalid_argument,
      ">>> ERROR (ROL::ChebyshevSpectral): Lower bound exceeds upper!");
    ROL_TEST_FOR_EXCEPTION(lower_ < static_cast<Real>(0), std::invalid_argument,
      ">>> ERROR (ROL::ChebyshevSpectral): Lower bound is less than zero!");
    ROL_TEST_FOR_EXCEPTION(static_cast<Real>(1) < upper_, std::invalid_argument,
      ">>> ERROR (ROL::ChebyshevSpectral): Upper bound is greater than one!");
    ROL_TEST_FOR_EXCEPTION((wType_ < 1 || wType_ > 3), std::invalid_argument,
      ">>> ERROR (ROL::ChebyshevSpectral): Weight must be 1, 2 or 3!");
    ROL_TEST_FOR_EXCEPTION(plusFunction_ == ROL::nullPtr, std::invalid_argument,
      ">>> ERROR (ROL::ChebyshevSpectral): PlusFunction pointer is null!");
  }

  void initializeQuad(void) {
     ROL::Ptr<Quadrature1D<Real> > quad;
     if ( wType_ == 1 ) {
       quad = ROL::makePtr<GaussChebyshev1Quadrature<Real>>(nQuad_);
     }
     else if ( wType_ == 2 ) {
       quad = ROL::makePtr<GaussChebyshev2Quadrature<Real>>(nQuad_);
     }
     else if ( wType_ == 3 ) {
       quad = ROL::makePtr<GaussChebyshev3Quadrature<Real>>(nQuad_);
     }
     // quad->test();
     quad->get(pts_,wts_);
     Real sum(0), half(0.5), one(1);
     for (int i = 0; i < nQuad_; ++i) {
       sum += wts_[i];
     }
     for (int i = 0; i < nQuad_; ++i) {
       wts_[i] /= sum;
       pts_[i] = lower_ + (upper_-lower_)*half*(pts_[i] + one);
     }
     SpectralRisk<Real>::buildMixedQuantile(pts_,wts_,plusFunction_);
  }

public:
  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Chebyshev-Spectral"
      and the "Chebyshev-Spectral" sublist should have the following parameters
      \li "Lower Bound" (between 0 and 1)
      \li "Upper Bound" (between 0 and 1, greater than "Lower Bound")
      \li "Weight Type" (either 1, 2, or 3)
      \li "Number of Quadrature Points"
      \li A sublist for plus function information.
  */
  ChebyshevSpectral( ROL::ParameterList &parlist )
    : SpectralRisk<Real>() {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Chebyshev Spectral Risk");
    // Grab confidence level and quadrature order
    lower_ = list.get("Lower Bound",0.0);
    upper_ = list.get("Upper Bound",1.0);
    nQuad_ = list.get("Number of Quadrature Points",5);
    wType_ = list.get("Weight Type",1);
    plusFunction_ = ROL::makePtr<PlusFunction<Real>>(list);
    // Check inputs
    checkInputs();
    initializeQuad();
  }

  /** \brief Constructor.

      @param[in]     lower   is the lower confidence level (between 0 and 1)
      @param[in]     upper   is the upper confidence level (between 0 and 1, greater than lower)
      @param[in]     nQuad   is the number of quadrature points
      @param[in]     wType   is the weight type (either 1, 2, or 3)
      @param[in]     pf      is the plus function or an approximation
  */
  ChebyshevSpectral(const Real lower, const Real upper,
                   const int nQuad, const int wType,
                   const ROL::Ptr<PlusFunction<Real> > &pf)
    : SpectralRisk<Real>(), plusFunction_(pf),
      lower_(lower), upper_(upper), nQuad_(nQuad), wType_(wType) {
    // Check inputs
    checkInputs();
    initializeQuad();
  }
};

}

#endif
