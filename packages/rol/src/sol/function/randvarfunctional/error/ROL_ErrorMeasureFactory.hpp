// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ERRORMEASUREFACTORY_HPP
#define ROL_ERRORMEASUREFACTORY_HPP

#include "ROL_ParameterList.hpp"

#include "ROL_Types.hpp"

// Risk Quadrangle Implementations
#include "ROL_ExpectationQuadError.hpp"
#include "ROL_GenMoreauYosidaCVaR.hpp"
#include "ROL_MoreauYosidaCVaR.hpp"
#include "ROL_LogExponentialQuadrangle.hpp"
#include "ROL_LogQuantileQuadrangle.hpp"
#include "ROL_MeanVarianceQuadrangle.hpp"
#include "ROL_QuantileQuadrangle.hpp"
#include "ROL_SmoothedWorstCaseQuadrangle.hpp"
#include "ROL_TruncatedMeanQuadrangle.hpp"

namespace ROL {

  enum EErrorMeasure {
    ERRORMEASURE_MEANVARIANCEQUADRANGLE = 0,
    ERRORMEASURE_TRUNCATEDMEANQUADRANGLE,
    ERRORMEASURE_QUANTILEQUADRANGLE,
    ERRORMEASURE_MOREAUYOSIDACVAR,
    ERRORMEASURE_GENMOREAUYOSIDACVAR,
    ERRORMEASURE_LOGEXPONENTIALQUADRANGLE,
    ERRORMEASURE_LOGQUANTILEQUADRANGLE,
    ERRORMEASURE_SMOOTHEDWORSTCASEQUADRANGLE,
    ERRORMEASURE_LAST
  };

  inline std::string EErrorMeasureToString(EErrorMeasure ed) {
    std::string retString;
    switch(ed) {
      case ERRORMEASURE_MEANVARIANCEQUADRANGLE:
             retString = "Least Squares";                             break;
      case ERRORMEASURE_TRUNCATEDMEANQUADRANGLE:
             retString = "Huber";                                     break;
      case ERRORMEASURE_QUANTILEQUADRANGLE:
             retString = "Koenker-Bassett";                           break;
      case ERRORMEASURE_MOREAUYOSIDACVAR:
             retString = "Moreau-Yosida-Koenker-Bassett";             break;
      case ERRORMEASURE_GENMOREAUYOSIDACVAR:
             retString = "Generalized Moreau-Yosida-Koenker-Bassett"; break;
      case ERRORMEASURE_LOGEXPONENTIALQUADRANGLE:
             retString = "Exponential";                               break;
      case ERRORMEASURE_LOGQUANTILEQUADRANGLE:
             retString = "Log Quantile";                              break;
      case ERRORMEASURE_SMOOTHEDWORSTCASEQUADRANGLE:
             retString = "Smoothed Worst Case";                       break;
      case ERRORMEASURE_LAST:
             retString = "Last Type (Dummy)";                         break;
      default:
             retString = "INVALID EErrorMeasure";                     break;
    }
    return retString;
  }

  inline int isValidErrorMeasure(EErrorMeasure ed) {
    return( (ed == ERRORMEASURE_MEANVARIANCEQUADRANGLE)      ||
            (ed == ERRORMEASURE_TRUNCATEDMEANQUADRANGLE)     || 
            (ed == ERRORMEASURE_QUANTILEQUADRANGLE)          ||
            (ed == ERRORMEASURE_MOREAUYOSIDACVAR)            ||
            (ed == ERRORMEASURE_GENMOREAUYOSIDACVAR)         ||
            (ed == ERRORMEASURE_LOGEXPONENTIALQUADRANGLE)    ||
            (ed == ERRORMEASURE_LOGQUANTILEQUADRANGLE)       ||
            (ed == ERRORMEASURE_SMOOTHEDWORSTCASEQUADRANGLE) ); 
  }

  inline EErrorMeasure & operator++(EErrorMeasure &type) {
    return type = static_cast<EErrorMeasure>(type+1);
  }

  inline EErrorMeasure operator++(EErrorMeasure &type, int) {
    EErrorMeasure oldval = type;
    ++type;
    return oldval;
  }

  inline EErrorMeasure & operator--(EErrorMeasure &type) {
    return type = static_cast<EErrorMeasure>(type-1);
  }

  inline EErrorMeasure operator--(EErrorMeasure &type, int) {
    EErrorMeasure oldval = type;
    --type;
    return oldval;
  }

  inline EErrorMeasure StringToEErrorMeasure(std::string s) {
    s = removeStringFormat(s);
    for ( EErrorMeasure tr = ERRORMEASURE_MEANVARIANCEQUADRANGLE; tr < ERRORMEASURE_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(EErrorMeasureToString(tr))) ) {
        return tr;
      }
    }
    return ERRORMEASURE_LAST;
  }

  template<class Real>
  inline Ptr<RandVarFunctional<Real>> ErrorMeasureFactory(ParameterList &parlist) {
    std::string error = parlist.sublist("SOL").sublist("Error Measure").get("Name","Least Squares");
    EErrorMeasure ed = StringToEErrorMeasure(error);
    switch(ed) {
      case ERRORMEASURE_MEANVARIANCEQUADRANGLE:
             return makePtr<ExpectationQuadError<Real>>(makePtr<MeanVarianceQuadrangle<Real>>(parlist));
      case ERRORMEASURE_TRUNCATEDMEANQUADRANGLE:
             return makePtr<ExpectationQuadError<Real>>(makePtr<TruncatedMeanQuadrangle<Real>>(parlist));
      case ERRORMEASURE_QUANTILEQUADRANGLE:
             return makePtr<ExpectationQuadError<Real>>(makePtr<QuantileQuadrangle<Real>>(parlist));
      case ERRORMEASURE_MOREAUYOSIDACVAR:
             return makePtr<ExpectationQuadError<Real>>(makePtr<MoreauYosidaCVaR<Real>>(parlist));
      case ERRORMEASURE_GENMOREAUYOSIDACVAR:
             return makePtr<ExpectationQuadError<Real>>(makePtr<GenMoreauYosidaCVaR<Real>>(parlist));
      case ERRORMEASURE_LOGEXPONENTIALQUADRANGLE:
             return makePtr<ExpectationQuadError<Real>>(makePtr<LogExponentialQuadrangle<Real>>(parlist));
      case ERRORMEASURE_LOGQUANTILEQUADRANGLE:
             return makePtr<ExpectationQuadError<Real>>(makePtr<LogQuantileQuadrangle<Real>>(parlist));
      case ERRORMEASURE_SMOOTHEDWORSTCASEQUADRANGLE:
             return makePtr<ExpectationQuadError<Real>>(makePtr<SmoothedWorstCaseQuadrangle<Real>>(parlist));
      default:
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                   "Invalid error measure type " << error << "!");
    }
  }
}
#endif
