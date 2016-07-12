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

#ifndef ROL_RISKMEASUREFACTORY_HPP
#define ROL_RISKMEASUREFACTORY_HPP

#include "Teuchos_ParameterList.hpp"

#include "ROL_Types.hpp"

// Standard Risk Measure Implementations
#include "ROL_CVaR.hpp"
#include "ROL_CoherentExpUtility.hpp"
#include "ROL_ExpUtility.hpp"
#include "ROL_HMCR.hpp"
#include "ROL_MeanDeviationFromTarget.hpp"
#include "ROL_MeanDeviation.hpp"
#include "ROL_MeanVarianceFromTarget.hpp"
#include "ROL_MeanVariance.hpp"
#include "ROL_MoreauYosidaCVaR.hpp"

// Risk Quadrangle Risk Measure Implementations
#include "ROL_LogExponentialQuadrangle.hpp"
#include "ROL_LogQuantileQuadrangle.hpp"
#include "ROL_MeanVarianceQuadrangle.hpp"
#include "ROL_MixedQuantileQuadrangle.hpp"
#include "ROL_SuperQuantileQuadrangle.hpp"
#include "ROL_Chebyshev1Kusuoka.hpp"
#include "ROL_Chebyshev2Kusuoka.hpp"
#include "ROL_Chebyshev3Kusuoka.hpp"
#include "ROL_QuantileQuadrangle.hpp"
#include "ROL_QuantileRadiusQuadrangle.hpp"
#include "ROL_SmoothedWorstCaseQuadrangle.hpp"
#include "ROL_TruncatedMeanQuadrangle.hpp"

// F-Divergence Distributionally Robust Risk Measure Implementations
#include "ROL_Chi2Divergence.hpp"
#include "ROL_KLDivergence.hpp"

namespace ROL {

  enum ERiskMeasure {
    RISKMEASURE_CVAR = 0,
    RISKMEASURE_COHERENTEXPUTILITY,
    RISKMEASURE_EXPUTILITY,
    RISKMEASURE_HMCR,
    RISKMEASURE_MEANDEVIATIONFROMTARGET, 
    RISKMEASURE_MEANDEVIATION,
    RISKMEASURE_MEANVARIANCEFROMTARGET,
    RISKMEASURE_MEANVARIANCE,
    RISKMEASURE_MOREAUYOSIDACVAR,
    RISKMEASURE_LOGEXPONENTIALQUADRANGLE,
    RISKMEASURE_LOGQUANTILEQUADRANGLE,
    RISKMEASURE_MEANVARIANCEQUADRANGLE,
    RISKMEASURE_MIXEDQUANTILEQUADRANGLE,
    RISKMEASURE_QUANTILEQUADRANGLE,
    RISKMEASURE_QUANTILERADIUSQUADRANGLE,
    RISKMEASURE_SMOOTHEDWORSTCASEQUADRANGLE,
    RISKMEASURE_SUPERQUANTILEQUADRANGLE,
    RISKMEASURE_CHEBYSHEV1KUSUOKA,
    RISKMEASURE_CHEBYSHEV2KUSUOKA,
    RISKMEASURE_CHEBYSHEV3KUSUOKA,
    RISKMEASURE_TRUNCATEDMEANQUADRANGLE,
    RISKMEASURE_CHI2DIVERGENCE,
    RISKMEASURE_KLDIVERGENCE,
    RISKMEASURE_LAST
  };

  inline std::string ERiskMeasureToString(ERiskMeasure ed) {
    std::string retString;
    switch(ed) {
      case RISKMEASURE_CVAR:
             retString = "CVaR";                                    break;
      case RISKMEASURE_COHERENTEXPUTILITY:
             retString = "Coherent Exponential Utility";            break;
      case RISKMEASURE_EXPUTILITY:
             retString = "Exponential Utility";                     break;
      case RISKMEASURE_HMCR:
             retString = "HMCR";                                    break;
      case RISKMEASURE_MEANDEVIATIONFROMTARGET:
             retString = "Mean Plus Deviation From Target";         break;
      case RISKMEASURE_MEANDEVIATION:
             retString = "Mean Plus Deviation";                     break;
      case RISKMEASURE_MEANVARIANCEFROMTARGET:
             retString = "Mean Plus Variance From Target";          break;
      case RISKMEASURE_MEANVARIANCE:
             retString = "Mean Plus Variance";                      break;
      case RISKMEASURE_MOREAUYOSIDACVAR:
             retString = "Moreau-Yosida CVaR";                      break;
      case RISKMEASURE_LOGEXPONENTIALQUADRANGLE:
             retString = "Log-Exponential Quadrangle";              break;
      case RISKMEASURE_LOGQUANTILEQUADRANGLE:
             retString = "Log-Quantile Quadrangle";                 break;
      case RISKMEASURE_MEANVARIANCEQUADRANGLE:
             retString = "Mean-Variance Quadrangle";                break;
      case RISKMEASURE_MIXEDQUANTILEQUADRANGLE:
             retString = "Mixed-Quantile Quadrangle";               break;
      case RISKMEASURE_SUPERQUANTILEQUADRANGLE:
             retString = "Super Quantile Quadrangle";               break;
      case RISKMEASURE_CHEBYSHEV1KUSUOKA:
             retString = "Chebyshev 1 Kusuoka";                     break;
      case RISKMEASURE_CHEBYSHEV2KUSUOKA:
             retString = "Chebyshev 2 Kusuoka";                     break;
      case RISKMEASURE_CHEBYSHEV3KUSUOKA:
             retString = "Chebyshev 3 Kusuoka";                     break;
      case RISKMEASURE_QUANTILEQUADRANGLE:
             retString = "Quantile-Based Quadrangle";               break;
      case RISKMEASURE_QUANTILERADIUSQUADRANGLE:
             retString = "Quantile-Radius Quadrangle";              break;
      case RISKMEASURE_SMOOTHEDWORSTCASEQUADRANGLE:
             retString = "Smoothed Worst-Case Quadrangle";          break;
      case RISKMEASURE_TRUNCATEDMEANQUADRANGLE:
             retString = "Truncated Mean Quadrangle";               break;
      case RISKMEASURE_CHI2DIVERGENCE:
             retString = "Chi-Squared Divergence";                  break;
      case RISKMEASURE_KLDIVERGENCE:
             retString = "KL Divergence";                           break;
      case RISKMEASURE_LAST:
             retString = "Last Type (Dummy)";                       break;
      default:
             retString = "INVALID ERiskMeasure";                    break;
    }
    return retString;
  }

  inline int isValidRiskMeasure(ERiskMeasure ed) {
    return( (ed == RISKMEASURE_CVAR) ||
            (ed == RISKMEASURE_COHERENTEXPUTILITY) ||
            (ed == RISKMEASURE_EXPUTILITY) ||
            (ed == RISKMEASURE_HMCR) ||
            (ed == RISKMEASURE_MEANDEVIATIONFROMTARGET) ||
            (ed == RISKMEASURE_MEANDEVIATION) ||
            (ed == RISKMEASURE_MEANVARIANCEFROMTARGET) ||
            (ed == RISKMEASURE_MEANVARIANCE) ||
            (ed == RISKMEASURE_MOREAUYOSIDACVAR) ||
            (ed == RISKMEASURE_LOGEXPONENTIALQUADRANGLE) ||
            (ed == RISKMEASURE_LOGQUANTILEQUADRANGLE) ||
            (ed == RISKMEASURE_MEANVARIANCEQUADRANGLE) ||
            (ed == RISKMEASURE_MIXEDQUANTILEQUADRANGLE) ||
            (ed == RISKMEASURE_SUPERQUANTILEQUADRANGLE) ||
            (ed == RISKMEASURE_CHEBYSHEV1KUSUOKA) ||
            (ed == RISKMEASURE_CHEBYSHEV2KUSUOKA) ||
            (ed == RISKMEASURE_CHEBYSHEV3KUSUOKA) ||
            (ed == RISKMEASURE_QUANTILEQUADRANGLE) ||
            (ed == RISKMEASURE_QUANTILERADIUSQUADRANGLE) ||
            (ed == RISKMEASURE_SMOOTHEDWORSTCASEQUADRANGLE) ||
            (ed == RISKMEASURE_TRUNCATEDMEANQUADRANGLE) ||
            (ed == RISKMEASURE_CHI2DIVERGENCE) ||
            (ed == RISKMEASURE_KLDIVERGENCE) );
  }

  inline ERiskMeasure & operator++(ERiskMeasure &type) {
    return type = static_cast<ERiskMeasure>(type+1);
  }

  inline ERiskMeasure operator++(ERiskMeasure &type, int) {
    ERiskMeasure oldval = type;
    ++type;
    return oldval;
  }

  inline ERiskMeasure & operator--(ERiskMeasure &type) {
    return type = static_cast<ERiskMeasure>(type-1);
  }

  inline ERiskMeasure operator--(ERiskMeasure &type, int) {
    ERiskMeasure oldval = type;
    --type;
    return oldval;
  }

  inline ERiskMeasure StringToERiskMeasure(std::string s) {
    s = removeStringFormat(s);
    for ( ERiskMeasure tr = RISKMEASURE_CVAR; tr < RISKMEASURE_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(ERiskMeasureToString(tr))) ) {
        return tr;
      }
    }
    return RISKMEASURE_LAST;
  }

  template<class Real>
  inline Teuchos::RCP<RiskMeasure<Real> > RiskMeasureFactory(Teuchos::ParameterList &parlist) {
    std::string risk = parlist.sublist("SOL").sublist("Risk Measure").get("Name","CVaR");
    ERiskMeasure ed = StringToERiskMeasure(risk);
    switch(ed) {
      case RISKMEASURE_CVAR:
             return Teuchos::rcp(new CVaR<Real>(parlist));
      case RISKMEASURE_COHERENTEXPUTILITY:
             return Teuchos::rcp(new CoherentExpUtility<Real>());
      case RISKMEASURE_EXPUTILITY:
             return Teuchos::rcp(new ExpUtility<Real>(parlist));
      case RISKMEASURE_HMCR:
             return Teuchos::rcp(new HMCR<Real>(parlist));
      case RISKMEASURE_MEANDEVIATIONFROMTARGET:
             return Teuchos::rcp(new MeanDeviationFromTarget<Real>(parlist));
      case RISKMEASURE_MEANDEVIATION:
             return Teuchos::rcp(new MeanDeviation<Real>(parlist));
      case RISKMEASURE_MEANVARIANCEFROMTARGET:
             return Teuchos::rcp(new MeanVarianceFromTarget<Real>(parlist));
      case RISKMEASURE_MEANVARIANCE:
             return Teuchos::rcp(new MeanVariance<Real>(parlist));
      case RISKMEASURE_MOREAUYOSIDACVAR:
             return Teuchos::rcp(new MoreauYosidaCVaR<Real>(parlist));
      case RISKMEASURE_LOGEXPONENTIALQUADRANGLE:
             return Teuchos::rcp(new LogExponentialQuadrangle<Real>(parlist));
      case RISKMEASURE_LOGQUANTILEQUADRANGLE:
             return Teuchos::rcp(new LogQuantileQuadrangle<Real>(parlist));
      case RISKMEASURE_MEANVARIANCEQUADRANGLE:
             return Teuchos::rcp(new MeanVarianceQuadrangle<Real>(parlist));
      case RISKMEASURE_MIXEDQUANTILEQUADRANGLE:
             return Teuchos::rcp(new MixedQuantileQuadrangle<Real>(parlist));
      case RISKMEASURE_SUPERQUANTILEQUADRANGLE:
             return Teuchos::rcp(new SuperQuantileQuadrangle<Real>(parlist));
      case RISKMEASURE_CHEBYSHEV1KUSUOKA:
             return Teuchos::rcp(new Chebyshev1Kusuoka<Real>(parlist));
      case RISKMEASURE_CHEBYSHEV2KUSUOKA:
             return Teuchos::rcp(new Chebyshev2Kusuoka<Real>(parlist));
      case RISKMEASURE_CHEBYSHEV3KUSUOKA:
             return Teuchos::rcp(new Chebyshev3Kusuoka<Real>(parlist));
      case RISKMEASURE_QUANTILEQUADRANGLE:
             return Teuchos::rcp(new QuantileQuadrangle<Real>(parlist));
      case RISKMEASURE_QUANTILERADIUSQUADRANGLE:
             return Teuchos::rcp(new QuantileRadiusQuadrangle<Real>(parlist));
      case RISKMEASURE_SMOOTHEDWORSTCASEQUADRANGLE:
             return Teuchos::rcp(new SmoothedWorstCaseQuadrangle<Real>(parlist));
      case RISKMEASURE_TRUNCATEDMEANQUADRANGLE:
             return Teuchos::rcp(new TruncatedMeanQuadrangle<Real>(parlist));
      case RISKMEASURE_CHI2DIVERGENCE:
             return Teuchos::rcp(new Chi2Divergence<Real>(parlist));
      case RISKMEASURE_KLDIVERGENCE:
             return Teuchos::rcp(new KLDivergence<Real>(parlist));
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                   "Invalid risk measure type " << risk << "!");
    }
  }
}
#endif
