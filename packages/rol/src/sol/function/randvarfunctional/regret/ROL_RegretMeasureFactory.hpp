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

#ifndef ROL_REGRETMEASUREFACTORY_HPP
#define ROL_REGRETMEASUREFACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "ROL_Types.hpp"
#include "ROL_ExpectationQuadRegret.hpp"
#include "ROL_QuantileQuadrangle.hpp"
#include "ROL_MoreauYosidaCVaR.hpp"
#include "ROL_GenMoreauYosidaCVaR.hpp"
#include "ROL_MeanVarianceQuadrangle.hpp"
#include "ROL_TruncatedMeanQuadrangle.hpp"
#include "ROL_LogExponentialQuadrangle.hpp"
#include "ROL_LogQuantileQuadrangle.hpp"
#include "ROL_SmoothedWorstCaseQuadrangle.hpp"

namespace ROL {

  enum ERegretMeasure {
    REGRETMEASURE_MEANABSOLUTELOSS = 0,
    REGRETMEASURE_MOREAUYOSIDAMEANABSOLUTELOSS,
    REGRETMEASURE_GENMOREAUYOSIDAMEANABSOLUTELOSS,
    REGRETMEASURE_EXPONENTIAL,
    REGRETMEASURE_MEANL2,
    REGRETMEASURE_TRUNCATEDMEAN,
    REGRETMEASURE_LOGQUANTILE,
    REGRETMEASURE_SMOOTHEDWORSTCASE,
    REGRETMEASURE_LAST
  };

  inline std::string ERegretMeasureToString(ERegretMeasure ed) {
    std::string retString;
    switch(ed) {
      case REGRETMEASURE_MEANABSOLUTELOSS:
             retString = "Mean Absolute Loss";                           break;
      case REGRETMEASURE_MOREAUYOSIDAMEANABSOLUTELOSS:
             retString = "Moreau-Yosida Mean Absolute Loss";             break;
      case REGRETMEASURE_GENMOREAUYOSIDAMEANABSOLUTELOSS:
             retString = "Generalized Moreau-Yosida Mean Absolute Loss"; break;
      case REGRETMEASURE_EXPONENTIAL:
             retString = "Exponential";                                  break;
      case REGRETMEASURE_MEANL2:
             retString = "Mean L2";                                      break;
      case REGRETMEASURE_TRUNCATEDMEAN:
             retString = "Truncated Mean";                               break;
      case REGRETMEASURE_LOGQUANTILE:
             retString = "Log Quantile";                                 break;
      case REGRETMEASURE_SMOOTHEDWORSTCASE:
             retString = "Smoothed Worst Case";                          break;
      case REGRETMEASURE_LAST:
             retString = "Last Type (Dummy)";                            break;
      default:
             retString = "INVALID ERegretMeasure";                       break;
    }
    return retString;
  }

  inline int isValidRegretMeasure(ERegretMeasure ed) {
    return( (ed == REGRETMEASURE_MEANABSOLUTELOSS)                ||
            (ed == REGRETMEASURE_MOREAUYOSIDAMEANABSOLUTELOSS)    ||
            (ed == REGRETMEASURE_GENMOREAUYOSIDAMEANABSOLUTELOSS) ||
            (ed == REGRETMEASURE_EXPONENTIAL)                     ||
            (ed == REGRETMEASURE_MEANL2)                          ||
            (ed == REGRETMEASURE_TRUNCATEDMEAN)                   ||
            (ed == REGRETMEASURE_LOGQUANTILE)                     ||
            (ed == REGRETMEASURE_SMOOTHEDWORSTCASE) );
  }

  inline ERegretMeasure & operator++(ERegretMeasure &type) {
    return type = static_cast<ERegretMeasure>(type+1);
  }

  inline ERegretMeasure operator++(ERegretMeasure &type, int) {
    ERegretMeasure oldval = type;
    ++type;
    return oldval;
  }

  inline ERegretMeasure & operator--(ERegretMeasure &type) {
    return type = static_cast<ERegretMeasure>(type-1);
  }

  inline ERegretMeasure operator--(ERegretMeasure &type, int) {
    ERegretMeasure oldval = type;
    --type;
    return oldval;
  }

  inline ERegretMeasure StringToERegretMeasure(std::string s) {
    s = removeStringFormat(s);
    for ( ERegretMeasure tr = REGRETMEASURE_MEANABSOLUTELOSS; tr < REGRETMEASURE_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(ERegretMeasureToString(tr))) ) {
        return tr;
      }
    }
    return REGRETMEASURE_LAST;
  }

  template<class Real>
  inline Ptr<RandVarFunctional<Real> > RegretMeasureFactory(Teuchos::ParameterList &parlist) {
    std::string regret = parlist.sublist("SOL").sublist("Regret Measure").get("Name","Mean Absolute Loss");
    ERegretMeasure ed = StringToERegretMeasure(regret);
    switch(ed) {
      case REGRETMEASURE_MEANABSOLUTELOSS:
             return makePtr<ExpectationQuadRegret<Real>>(makePtr<QuantileQuadrangle<Real>>(parlist));
      case REGRETMEASURE_MOREAUYOSIDAMEANABSOLUTELOSS:
             return makePtr<ExpectationQuadRegret<Real>>(makePtr<MoreauYosidaCVaR<Real>>(parlist));
      case REGRETMEASURE_GENMOREAUYOSIDAMEANABSOLUTELOSS:
             return makePtr<ExpectationQuadRegret<Real>>(makePtr<GenMoreauYosidaCVaR<Real>>(parlist));
      case REGRETMEASURE_EXPONENTIAL:
             return makePtr<ExpectationQuadRegret<Real>>(makePtr<LogExponentialQuadrangle<Real>>(parlist));
      case REGRETMEASURE_MEANL2:
             return makePtr<ExpectationQuadRegret<Real>>(makePtr<MeanVarianceQuadrangle<Real>>(parlist));
      case REGRETMEASURE_TRUNCATEDMEAN:
             return makePtr<ExpectationQuadRegret<Real>>(makePtr<TruncatedMeanQuadrangle<Real>>(parlist));
      case REGRETMEASURE_LOGQUANTILE:
             return makePtr<ExpectationQuadRegret<Real>>(makePtr<LogQuantileQuadrangle<Real>>(parlist));
      case REGRETMEASURE_SMOOTHEDWORSTCASE:
             return makePtr<ExpectationQuadRegret<Real>>(makePtr<SmoothedWorstCaseQuadrangle<Real>>(parlist));
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                   "Invalid regret measure type " << regret << "!");
    }
  }
}
#endif
