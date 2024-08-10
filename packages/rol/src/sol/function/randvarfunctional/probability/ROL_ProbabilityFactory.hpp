// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROBABILITYFACTORY_HPP
#define ROL_PROBABILITYFACTORY_HPP

#include "ROL_ParameterList.hpp"

#include "ROL_Types.hpp"
#include "ROL_BPOE.hpp"
#include "ROL_SmoothedPOE.hpp"

namespace ROL {

  enum EProbability {
    PROBABILITY_BPOE = 0,
    PROBABILITY_SMOOTHEDPOE,
    PROBABILITY_LAST
  };

  inline std::string EProbabilityToString(EProbability ed) {
    std::string retString;
    switch(ed) {
      case PROBABILITY_BPOE:
             retString = "bPOE";                                    break;
      case PROBABILITY_SMOOTHEDPOE:
             retString = "Smoothed POE";                            break;
      case PROBABILITY_LAST:
             retString = "Last Type (Dummy)";                       break;
      default:
             retString = "INVALID EProbability";                    break;
    }
    return retString;
  }

  inline int isValidProbability(EProbability ed) {
    return( (ed == PROBABILITY_BPOE)                    ||
            (ed == PROBABILITY_SMOOTHEDPOE));
  }

  inline EProbability & operator++(EProbability &type) {
    return type = static_cast<EProbability>(type+1);
  }

  inline EProbability operator++(EProbability &type, int) {
    EProbability oldval = type;
    ++type;
    return oldval;
  }

  inline EProbability & operator--(EProbability &type) {
    return type = static_cast<EProbability>(type-1);
  }

  inline EProbability operator--(EProbability &type, int) {
    EProbability oldval = type;
    --type;
    return oldval;
  }

  inline EProbability StringToEProbability(std::string s) {
    s = removeStringFormat(s);
    for ( EProbability tr = PROBABILITY_BPOE; tr < PROBABILITY_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(EProbabilityToString(tr))) ) {
        return tr;
      }
    }
    return PROBABILITY_LAST;
  }

  template<class Real>
  inline Ptr<RandVarFunctional<Real>> ProbabilityFactory(ParameterList &parlist) {
    std::string prob = parlist.sublist("SOL").sublist("Probability").get("Name","bPOE");
    EProbability ed = StringToEProbability(prob);
    switch(ed) {
      case PROBABILITY_BPOE:
             return makePtr<BPOE<Real>>(parlist);
      case PROBABILITY_SMOOTHEDPOE:
             return makePtr<SmoothedPOE<Real>>(parlist);
      default:
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                   "Invalid probability type " << prob << "!");
    }
  }
}
#endif
