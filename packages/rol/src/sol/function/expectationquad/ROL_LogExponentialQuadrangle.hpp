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

#ifndef ROL_LOGEXPONENTIALQUAD_HPP
#define ROL_LOGEXPONENTIALQUAD_HPP

#include "ROL_ExpectationQuad.hpp"

/** @ingroup risk_group
    \class ROL::LogExponentialQuadrangle
    \brief Provides an interface for the entropic risk using the expectation
           risk quadrangle.

    The entropic risk measure (also called the exponential utility and the
    log-exponential risk measure) is
    \f[
       \mathcal{R}(X) = \lambda
       \log\mathbb{E}\left[\exp\left(\frac{X}{\lambda}\right)\right]
    \f]
    for \f$\lambda > 0\f$.  The entropic risk is convex, translation
    equivariant and monotonic.

    This class defines the entropic risk measure using the framework of the
    expectation risk quadrangle.  In this case, the scalar regret function
    is
    \f[
       v(x) = \lambda(\exp\left(\frac{x}{\lambda}\right)-1).
    \f]
    The entropic risk measure is then implemented as
    \f[
       \mathcal{R}(X) = \inf_{t\in\mathbb{R}}\left\{
           t + \mathbb{E}[v(X-t)] \right\}.
    \f]
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$t\f$, then minimizes jointly for \f$(x_0,t)\f$.
*/

namespace ROL {

template<class Real>
class LogExponentialQuadrangle : public ExpectationQuad<Real> {
private:
  Real coeff_;

  void parseParameterList(Teuchos::ParameterList &parlist) {
    std::string type = parlist.sublist("SOL").get("Stochastic Component Type","Risk Averse");
    Teuchos::ParameterList list;
    if (type == "Risk Averse") {
      list = parlist.sublist("SOL").sublist("Risk Measure").sublist("Entropic Risk");
    }
    else if (type == "Error") {
      list = parlist.sublist("SOL").sublist("Error Measure").sublist("Exponential");
    }
    else if (type == "Deviation") {
      list = parlist.sublist("SOL").sublist("Deviation Measure").sublist("Entropic");
    }
    else if (type == "Regret") {
      list = parlist.sublist("SOL").sublist("Regret Measure").sublist("Exponential");
    }
    coeff_ = list.get<Real>("Rate");
  }

  void checkInputs(void) const {
    Real zero(0);
    TEUCHOS_TEST_FOR_EXCEPTION((coeff_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::LogExponentialQuadrangle): Rate must be positive!");
  }

public:
  /** \brief Constructor.

      @param[in]     coeff    is the scale parameter \f$\lambda\f$
  */
  LogExponentialQuadrangle(const Real coeff = 1)
    : ExpectationQuad<Real>(), coeff_(coeff) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measures"->"Log-Exponential Quadrangle"
      and withing the "Log-Exponential Quadrangle" sublist should have
      \li "Rate" (greater than 0). 
  */
  LogExponentialQuadrangle(Teuchos::ParameterList &parlist)
    : ExpectationQuad<Real>() {
    parseParameterList(parlist);
    checkInputs();
  }

  Real error(Real x, int deriv = 0) {
    Real err(0), one(1), cx = coeff_*x;
    if (deriv==0) {
      err = (std::exp(cx) - cx - one)/coeff_;
    }
    else if (deriv==1) {
      err = std::exp(cx) - one;
    }
    else {
      err = coeff_*std::exp(cx);
    }
    return err;
  }

  Real regret(Real x, int deriv = 0) {
    Real zero(0), one(1);
    Real X = ((deriv==0) ? x : ((deriv==1) ? one : zero));
    Real reg = error(x,deriv) + X;
    return reg;
  }

};

}
#endif
