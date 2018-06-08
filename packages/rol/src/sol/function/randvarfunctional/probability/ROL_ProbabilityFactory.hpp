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
  inline Ptr<RandVarFunctional<Real> > ProbabilityFactory(ROL::ParameterList &parlist) {
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
