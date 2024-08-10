// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REGRETMEASUREFACTORY_HPP
#define ROL_REGRETMEASUREFACTORY_HPP

#include "ROL_ParameterList.hpp"
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
  inline Ptr<RandVarFunctional<Real>> RegretMeasureFactory(ParameterList &parlist) {
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
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                   "Invalid regret measure type " << regret << "!");
    }
  }
}
#endif
