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
  inline Ptr<RandVarFunctional<Real> > DeviationMeasureFactory(ROL::ParameterList &parlist) {
    std::string risk = parlist.sublist("SOL").sublist("Deviation Measure").get("Name","Least Squares");
    EDeviationMeasure ed = StringToEDeviationMeasure(risk);
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
                                   "Invalid deviation measure type " << risk << "!");
    }
  }
}
#endif
