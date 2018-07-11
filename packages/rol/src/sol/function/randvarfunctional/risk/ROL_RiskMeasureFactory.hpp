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

#include "ROL_ParameterList.hpp"

#include "ROL_Types.hpp"

// Standard Risk Measure Implementations
#include "ROL_CVaR.hpp"
#include "ROL_MixedCVaR.hpp"
#include "ROL_SecondOrderCVaR.hpp"
#include "ROL_ChebyshevSpectral.hpp"
#include "ROL_SpectralRisk.hpp"
#include "ROL_QuantileRadius.hpp"
#include "ROL_HMCR.hpp"
#include "ROL_EntropicRisk.hpp"
#include "ROL_CoherentEntropicRisk.hpp"
#include "ROL_MeanDeviationFromTarget.hpp"
#include "ROL_MeanDeviation.hpp"
#include "ROL_MeanVarianceFromTarget.hpp"
#include "ROL_MeanVariance.hpp"

// Risk Quadrangle Risk Measure Implementations
#include "ROL_ExpectationQuadRisk.hpp"
#include "ROL_LogQuantileQuadrangle.hpp"
#include "ROL_SmoothedWorstCaseQuadrangle.hpp"
#include "ROL_TruncatedMeanQuadrangle.hpp"
#include "ROL_MoreauYosidaCVaR.hpp"
#include "ROL_GenMoreauYosidaCVaR.hpp"
#include "ROL_LogExponentialQuadrangle.hpp"
#include "ROL_MeanVarianceQuadrangle.hpp"

// F-Divergence Distributionally Robust Risk Measure Implementations
#include "ROL_Chi2Divergence.hpp"
#include "ROL_KLDivergence.hpp"

namespace ROL {

  enum ERiskMeasure {
    RISKMEASURE_CVAR = 0,
    RISKMEASURE_MOREAUYOSIDACVAR,
    RISKMEASURE_GENMOREAUYOSIDACVAR,
    RISKMEASURE_MIXEDCVAR,
    RISKMEASURE_SPECTRALRISK,
    RISKMEASURE_SECONDORDERCVAR,
    RISKMEASURE_CHEBYSHEVSPECTRAL,
    RISKMEASURE_QUANTILERADIUS,
    RISKMEASURE_HMCR,
    RISKMEASURE_ENTROPICRISK,
    RISKMEASURE_COHERENTENTROPICRISK,
    RISKMEASURE_MEANDEVIATIONFROMTARGET, 
    RISKMEASURE_MEANDEVIATION,
    RISKMEASURE_MEANVARIANCEFROMTARGET,
    RISKMEASURE_MEANVARIANCE,
    RISKMEASURE_TRUNCATEDMEAN,
    RISKMEASURE_LOGQUANTILE,
    RISKMEASURE_SMOOTHEDWORSTCASE,
    RISKMEASURE_LOGEXPONENTIAL,
    RISKMEASURE_SAFETYMARGIN,
    RISKMEASURE_CHI2DIVERGENCE,
    RISKMEASURE_KLDIVERGENCE,
    RISKMEASURE_LAST
  };

  inline std::string ERiskMeasureToString(ERiskMeasure ed) {
    std::string retString;
    switch(ed) {
      case RISKMEASURE_CVAR:
             retString = "CVaR";                                    break;
      case RISKMEASURE_MOREAUYOSIDACVAR:
             retString = "Moreau-Yosida CVaR";                      break;
      case RISKMEASURE_GENMOREAUYOSIDACVAR:
             retString = "Generalized Moreau-Yosida CVaR";          break;
      case RISKMEASURE_MIXEDCVAR:
             retString = "Mixed CVaR";                              break;
      case RISKMEASURE_SPECTRALRISK:
             retString = "Spectral Risk";                           break;
      case RISKMEASURE_SECONDORDERCVAR:
             retString = "Second Order CVaR";                       break;
      case RISKMEASURE_CHEBYSHEVSPECTRAL:
             retString = "Chebyshev Spectral Risk";                 break;
      case RISKMEASURE_QUANTILERADIUS:
             retString = "Quantile Radius";                         break;
      case RISKMEASURE_HMCR:
             retString = "HMCR";                                    break;
      case RISKMEASURE_ENTROPICRISK:
             retString = "Entropic Risk";                           break;
      case RISKMEASURE_COHERENTENTROPICRISK:
             retString = "Coherent Entropic Risk";                  break;
      case RISKMEASURE_MEANDEVIATIONFROMTARGET:
             retString = "Mean Plus Deviation From Target";         break;
      case RISKMEASURE_MEANDEVIATION:
             retString = "Mean Plus Deviation";                     break;
      case RISKMEASURE_MEANVARIANCEFROMTARGET:
             retString = "Mean Plus Variance From Target";          break;
      case RISKMEASURE_MEANVARIANCE:
             retString = "Mean Plus Variance";                      break;
      case RISKMEASURE_TRUNCATEDMEAN:
             retString = "Truncated Mean";                          break;
      case RISKMEASURE_LOGQUANTILE:
             retString = "Log Quantile";                            break;
      case RISKMEASURE_SMOOTHEDWORSTCASE:
             retString = "Smoothed Worst Case";                     break;
      case RISKMEASURE_LOGEXPONENTIAL:
             retString = "Log Exponential";                         break;
      case RISKMEASURE_SAFETYMARGIN:
             retString = "Safety Margin";                           break;
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
    return( (ed == RISKMEASURE_CVAR)                    ||
            (ed == RISKMEASURE_MOREAUYOSIDACVAR)        ||
            (ed == RISKMEASURE_GENMOREAUYOSIDACVAR)     ||
            (ed == RISKMEASURE_MIXEDCVAR)               ||
            (ed == RISKMEASURE_SPECTRALRISK)            ||
            (ed == RISKMEASURE_SECONDORDERCVAR)         ||
            (ed == RISKMEASURE_CHEBYSHEVSPECTRAL)       ||
            (ed == RISKMEASURE_QUANTILERADIUS)          ||
            (ed == RISKMEASURE_HMCR)                    ||
            (ed == RISKMEASURE_ENTROPICRISK)            ||
            (ed == RISKMEASURE_COHERENTENTROPICRISK)    ||
            (ed == RISKMEASURE_MEANDEVIATIONFROMTARGET) ||
            (ed == RISKMEASURE_MEANDEVIATION)           ||
            (ed == RISKMEASURE_MEANVARIANCEFROMTARGET)  ||
            (ed == RISKMEASURE_MEANVARIANCE)            ||
            (ed == RISKMEASURE_TRUNCATEDMEAN)           ||
            (ed == RISKMEASURE_LOGQUANTILE)             ||
            (ed == RISKMEASURE_SMOOTHEDWORSTCASE)       ||
            (ed == RISKMEASURE_LOGEXPONENTIAL)          ||
            (ed == RISKMEASURE_SAFETYMARGIN)            ||
            (ed == RISKMEASURE_CHI2DIVERGENCE)          ||
            (ed == RISKMEASURE_KLDIVERGENCE));
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
  inline Ptr<RandVarFunctional<Real> > RiskMeasureFactory(ROL::ParameterList &parlist) {
    std::string risk = parlist.sublist("SOL").sublist("Risk Measure").get("Name","CVaR");
    ERiskMeasure ed = StringToERiskMeasure(risk);
    switch(ed) {
      case RISKMEASURE_CVAR:
             return makePtr<CVaR<Real>>(parlist);
      case RISKMEASURE_MOREAUYOSIDACVAR:
             return makePtr<ExpectationQuadRisk<Real>>(makePtr<MoreauYosidaCVaR<Real>>(parlist));
      case RISKMEASURE_GENMOREAUYOSIDACVAR:
             return makePtr<ExpectationQuadRisk<Real>>(makePtr<GenMoreauYosidaCVaR<Real>>(parlist));
      case RISKMEASURE_MIXEDCVAR:
             return makePtr<MixedCVaR<Real>>(parlist);
      case RISKMEASURE_SPECTRALRISK:
             return makePtr<SpectralRisk<Real>>(parlist);
      case RISKMEASURE_SECONDORDERCVAR:
             return makePtr<SecondOrderCVaR<Real>>(parlist);
      case RISKMEASURE_CHEBYSHEVSPECTRAL:
             return makePtr<ChebyshevSpectral<Real>>(parlist);
      case RISKMEASURE_QUANTILERADIUS:
             return makePtr<QuantileRadius<Real>>(parlist);
      case RISKMEASURE_HMCR:
             return makePtr<HMCR<Real>>(parlist);
      case RISKMEASURE_ENTROPICRISK:
             return makePtr<EntropicRisk<Real>>(parlist);
      case RISKMEASURE_COHERENTENTROPICRISK:
             return makePtr<CoherentEntropicRisk<Real>>();
      case RISKMEASURE_MEANDEVIATIONFROMTARGET:
             return makePtr<MeanDeviationFromTarget<Real>>(parlist);
      case RISKMEASURE_MEANDEVIATION:
             return makePtr<MeanDeviation<Real>>(parlist);
      case RISKMEASURE_MEANVARIANCEFROMTARGET:
             return makePtr<MeanVarianceFromTarget<Real>>(parlist);
      case RISKMEASURE_MEANVARIANCE:
             return makePtr<MeanVariance<Real>>(parlist);
      case RISKMEASURE_TRUNCATEDMEAN:
             return makePtr<ExpectationQuadRisk<Real>>(makePtr<TruncatedMeanQuadrangle<Real>>(parlist));
      case RISKMEASURE_LOGQUANTILE:
             return makePtr<ExpectationQuadRisk<Real>>(makePtr<LogQuantileQuadrangle<Real>>(parlist));
      case RISKMEASURE_SMOOTHEDWORSTCASE:
             return makePtr<ExpectationQuadRisk<Real>>(makePtr<SmoothedWorstCaseQuadrangle<Real>>(parlist));
      case RISKMEASURE_LOGEXPONENTIAL:
             return makePtr<ExpectationQuadRisk<Real>>(makePtr<LogExponentialQuadrangle<Real>>(parlist));
      case RISKMEASURE_SAFETYMARGIN:
             return makePtr<ExpectationQuadRisk<Real>>(makePtr<MeanVarianceQuadrangle<Real>>(parlist));
      case RISKMEASURE_CHI2DIVERGENCE:
             return makePtr<Chi2Divergence<Real>>(parlist);
      case RISKMEASURE_KLDIVERGENCE:
             return makePtr<KLDivergence<Real>>(parlist);
      default:
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                   "Invalid risk measure type " << risk << "!");
    }
  }
}
#endif
