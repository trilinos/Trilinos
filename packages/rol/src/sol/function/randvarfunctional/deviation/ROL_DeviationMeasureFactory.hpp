// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DEVIATIONMEASUREFACTORY_HPP
#define ROL_DEVIATIONMEASUREFACTORY_HPP

#include "ROL_ParameterList.hpp"

#include "ROL_Types.hpp"

// Risk Quadrangle Implementations
#include "ROL_ExpectationQuadDeviation.hpp"
#include "ROL_GenMoreauYosidaCVaR.hpp"
#include "ROL_MoreauYosidaCVaR.hpp"
#include "ROL_LogExponentialQuadrangle.hpp"
#include "ROL_LogQuantileQuadrangle.hpp"
#include "ROL_MeanVarianceQuadrangle.hpp"
#include "ROL_QuantileQuadrangle.hpp"
#include "ROL_SmoothedWorstCaseQuadrangle.hpp"
#include "ROL_TruncatedMeanQuadrangle.hpp"

namespace ROL {

  enum EDeviationMeasure {
    DEVIATIONMEASURE_MEANVARIANCEQUADRANGLE = 0,
    DEVIATIONMEASURE_TRUNCATEDMEANQUADRANGLE,
    DEVIATIONMEASURE_QUANTILEQUADRANGLE,
    DEVIATIONMEASURE_MOREAUYOSIDACVAR,
    DEVIATIONMEASURE_GENMOREAUYOSIDACVAR,
    DEVIATIONMEASURE_LOGEXPONENTIALQUADRANGLE,
    DEVIATIONMEASURE_LOGQUANTILEQUADRANGLE,
    DEVIATIONMEASURE_SMOOTHEDWORSTCASEQUADRANGLE,
    DEVIATIONMEASURE_LAST
  };

  inline std::string EDeviationMeasureToString(EDeviationMeasure ed) {
    std::string retString;
    switch(ed) {
      case DEVIATIONMEASURE_MEANVARIANCEQUADRANGLE:
             retString = "Variance";                                  break;
      case DEVIATIONMEASURE_TRUNCATEDMEANQUADRANGLE:
             retString = "Truncated Mean";                            break;
      case DEVIATIONMEASURE_QUANTILEQUADRANGLE:
             retString = "CVaR";                                      break;
      case DEVIATIONMEASURE_MOREAUYOSIDACVAR:
             retString = "Moreau-Yosida CVaR";                        break;
      case DEVIATIONMEASURE_GENMOREAUYOSIDACVAR:
             retString = "Generalized Moreau-Yosida CVaR";            break;
      case DEVIATIONMEASURE_LOGEXPONENTIALQUADRANGLE:
             retString = "Entropic";                                  break;
      case DEVIATIONMEASURE_LOGQUANTILEQUADRANGLE:
             retString = "Log Quantile";                              break;
      case DEVIATIONMEASURE_SMOOTHEDWORSTCASEQUADRANGLE:
             retString = "Smoothed Upper Range";                      break;
      case DEVIATIONMEASURE_LAST:
             retString = "Last Type (Dummy)";                         break;
      default:
             retString = "INVALID EDeviationMeasure";                 break;
    }
    return retString;
  }

  inline int isValidDeviationMeasure(EDeviationMeasure ed) {
    return( (ed == DEVIATIONMEASURE_MEANVARIANCEQUADRANGLE)      ||
            (ed == DEVIATIONMEASURE_TRUNCATEDMEANQUADRANGLE)     || 
            (ed == DEVIATIONMEASURE_QUANTILEQUADRANGLE)          ||
            (ed == DEVIATIONMEASURE_MOREAUYOSIDACVAR)            ||
            (ed == DEVIATIONMEASURE_GENMOREAUYOSIDACVAR)         ||
            (ed == DEVIATIONMEASURE_LOGEXPONENTIALQUADRANGLE)    ||
            (ed == DEVIATIONMEASURE_LOGQUANTILEQUADRANGLE)       ||
            (ed == DEVIATIONMEASURE_SMOOTHEDWORSTCASEQUADRANGLE) ); 
  }

  inline EDeviationMeasure & operator++(EDeviationMeasure &type) {
    return type = static_cast<EDeviationMeasure>(type+1);
  }

  inline EDeviationMeasure operator++(EDeviationMeasure &type, int) {
    EDeviationMeasure oldval = type;
    ++type;
    return oldval;
  }

  inline EDeviationMeasure & operator--(EDeviationMeasure &type) {
    return type = static_cast<EDeviationMeasure>(type-1);
  }

  inline EDeviationMeasure operator--(EDeviationMeasure &type, int) {
    EDeviationMeasure oldval = type;
    --type;
    return oldval;
  }

  inline EDeviationMeasure StringToEDeviationMeasure(std::string s) {
    s = removeStringFormat(s);
    for ( EDeviationMeasure tr = DEVIATIONMEASURE_MEANVARIANCEQUADRANGLE; tr < DEVIATIONMEASURE_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(EDeviationMeasureToString(tr))) ) {
        return tr;
      }
    }
    return DEVIATIONMEASURE_LAST;
  }

  template<class Real>
  inline Ptr<RandVarFunctional<Real>> DeviationMeasureFactory(ParameterList &parlist) {
    std::string deviation = parlist.sublist("SOL").sublist("Deviation Measure").get("Name","Variance");
    EDeviationMeasure ed = StringToEDeviationMeasure(deviation);
    switch(ed) {
      case DEVIATIONMEASURE_MEANVARIANCEQUADRANGLE:
             return makePtr<ExpectationQuadDeviation<Real>>(makePtr<MeanVarianceQuadrangle<Real>>(parlist));
      case DEVIATIONMEASURE_TRUNCATEDMEANQUADRANGLE:
             return makePtr<ExpectationQuadDeviation<Real>>(makePtr<TruncatedMeanQuadrangle<Real>>(parlist));
      case DEVIATIONMEASURE_QUANTILEQUADRANGLE:
             return makePtr<ExpectationQuadDeviation<Real>>(makePtr<QuantileQuadrangle<Real>>(parlist));
      case DEVIATIONMEASURE_MOREAUYOSIDACVAR:
             return makePtr<ExpectationQuadDeviation<Real>>(makePtr<MoreauYosidaCVaR<Real>>(parlist));
      case DEVIATIONMEASURE_GENMOREAUYOSIDACVAR:
             return makePtr<ExpectationQuadDeviation<Real>>(makePtr<GenMoreauYosidaCVaR<Real>>(parlist));
      case DEVIATIONMEASURE_LOGEXPONENTIALQUADRANGLE:
             return makePtr<ExpectationQuadDeviation<Real>>(makePtr<LogExponentialQuadrangle<Real>>(parlist));
      case DEVIATIONMEASURE_LOGQUANTILEQUADRANGLE:
             return makePtr<ExpectationQuadDeviation<Real>>(makePtr<LogQuantileQuadrangle<Real>>(parlist));
      case DEVIATIONMEASURE_SMOOTHEDWORSTCASEQUADRANGLE:
             return makePtr<ExpectationQuadDeviation<Real>>(makePtr<SmoothedWorstCaseQuadrangle<Real>>(parlist));
      default:
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                   "Invalid deviation measure type " << deviation << "!");
    }
  }
}
#endif
