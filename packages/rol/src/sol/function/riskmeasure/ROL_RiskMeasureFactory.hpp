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
#include "ROL_ExpUtility.hpp"
#include "ROL_HMCR.hpp"
#include "ROL_MeanDeviationFromTarget.hpp"
#include "ROL_MeanDeviation.hpp"
#include "ROL_MeanVarianceFromTarget.hpp"
#include "ROL_MeanVariance.hpp"

// Risk Quadrangle Risk Measure Implementations
#include "ROL_MoreauYosidaQuantileQuadrangle.hpp"
#include "ROL_LogExponentialQuadrangle.hpp"
#include "ROL_LogQuantileQuadrangle.hpp"
#include "ROL_QuantileQuadrangle.hpp"
#include "ROL_TruncatedMeanQuadrangle.hpp"

namespace ROL {

  enum ERiskMeasure {
    RISKMEASURE_CVAR = 0,
    RISKMEASURE_EXPUTILITY,
    RISKMEASURE_HMCR,
    RISKMEASURE_MEANDEVIATIONFROMTARGET, 
    RISKMEASURE_MEANDEVIATION,
    RISKMEASURE_MEANVARIANCEFROMTARGET,
    RISKMEASURE_MEANVARIANCE,
    RISKMEASURE_MOREAUYOSIDAQUANTILEQUADRANGLE,
    RISKMEASURE_LOGEXPONENTIALQUADRANGLE,
    RISKMEASURE_LOGQUANTILEQUADRANGLE,
    RISKMEASURE_QUANTILEQUADRANGLE,
    RISKMEASURE_TRUNCATEDMEANQUADRANGLE,
    RISKMEASURE_LAST
  };

  inline std::string ERiskMeasureToString(ERiskMeasure ed) {
    std::string retString;
    switch(ed) {
      case RISKMEASURE_CVAR:
             retString = "CVaR";                                    break;
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
      case RISKMEASURE_MOREAUYOSIDAQUANTILEQUADRANGLE:
             retString = "Moreau-Yosida Quantile-Based Quadrangle"; break;
      case RISKMEASURE_LOGEXPONENTIALQUADRANGLE:
             retString = "Log-Exponential Quadrangle";              break;
      case RISKMEASURE_LOGQUANTILEQUADRANGLE:
             retString = "Log-Quantile Quadrangle";                 break;
      case RISKMEASURE_QUANTILEQUADRANGLE:
             retString = "Quantile-Based Quadrangle";               break;
      case RISKMEASURE_TRUNCATEDMEANQUADRANGLE:
             retString = "Truncated Mean Quadrangle";               break;
      case RISKMEASURE_LAST:
             retString = "Last Type (Dummy)";                       break;
      default:
             retString = "INVALID ERiskMeasure";                    break;
    }
    return retString;
  }

  inline int isValidRiskMeasure(ERiskMeasure ed) {
    return( (ed == RISKMEASURE_CVAR) ||
            (ed == RISKMEASURE_EXPUTILITY) ||
            (ed == RISKMEASURE_MEANDEVIATIONFROMTARGET) ||
            (ed == RISKMEASURE_MEANDEVIATION) ||
            (ed == RISKMEASURE_MEANVARIANCEFROMTARGET) ||
            (ed == RISKMEASURE_MEANVARIANCE) ||
            (ed == RISKMEASURE_MOREAUYOSIDAQUANTILEQUADRANGLE) ||
            (ed == RISKMEASURE_LOGEXPONENTIALQUADRANGLE) ||
            (ed == RISKMEASURE_LOGQUANTILEQUADRANGLE) ||
            (ed == RISKMEASURE_QUANTILEQUADRANGLE) ||
            (ed == RISKMEASURE_TRUNCATEDMEANQUADRANGLE) );
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
    std::string dist = parlist.sublist("SOL").sublist("Risk Measure").get("Name","CVaR");
    ERiskMeasure ed = StringToERiskMeasure(dist);
    switch(ed) {
      case RISKMEASURE_CVAR:
             return Teuchos::rcp(new CVaR<Real>(parlist));
      case RISKMEASURE_EXPUTILITY:
             return Teuchos::rcp(new ExpUtility<Real>);
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
      case RISKMEASURE_MOREAUYOSIDAQUANTILEQUADRANGLE:
             return Teuchos::rcp(new MoreauYosidaQuantileQuadrangle<Real>(parlist));
      case RISKMEASURE_LOGEXPONENTIALQUADRANGLE:
             return Teuchos::rcp(new LogExponentialQuadrangle<Real>);
      case RISKMEASURE_LOGQUANTILEQUADRANGLE:
             return Teuchos::rcp(new LogQuantileQuadrangle<Real>(parlist));
      case RISKMEASURE_QUANTILEQUADRANGLE:
             return Teuchos::rcp(new QuantileQuadrangle<Real>(parlist));
      case RISKMEASURE_TRUNCATEDMEANQUADRANGLE:
             return Teuchos::rcp(new TruncatedMeanQuadrangle<Real>(parlist));
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                   "Invalid risk measure type" << dist);
    }
  }
}
#endif
