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

#ifndef ROL_TYPEU_ALGORITHM_FACTORY_H
#define ROL_TYPEU_ALGORITHM_FACTORY_H

#include "ROL_TypeU_BundleAlgorithm.hpp"
#include "ROL_TypeU_LineSearchAlgorithm.hpp"
#include "ROL_TypeU_TrustRegionAlgorithm.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace TypeU {

/** \enum   ROL::EAlgorithmU
    \brief  Enumeration of unconstrained algorithm types.

    \arg    ALGORITHM_U_BUNDLE         describe
    \arg    ALGORITHM_U_LINESEARCH     describe
    \arg    ALGORITHM_U_TRUSTREGION    describe
 */
enum EAlgorithmU{
  ALGORITHM_U_BUNDLE = 0,
  ALGORITHM_U_LINESEARCH,
  ALGORITHM_U_TRUSTREGION,
  ALGORITHM_U_LAST
};

inline std::string EAlgorithmUToString(EAlgorithmU alg) {
  std::string retString;
  switch(alg) {
    case ALGORITHM_U_BUNDLE:      retString = "Bundle";            break;
    case ALGORITHM_U_LINESEARCH:  retString = "Line Search";       break;
    case ALGORITHM_U_TRUSTREGION: retString = "Trust Region";      break;
    case ALGORITHM_U_LAST:        retString = "Last Type (Dummy)"; break;
    default:                      retString = "INVALID EAlgorithmU";
  }
  return retString;
}

/** \brief  Verifies validity of a AlgorithmU enum.
  
    \param  ls  [in]  - enum of the AlgorithmU
    \return 1 if the argument is a valid AlgorithmU; 0 otherwise.
  */
inline int isValidAlgorithmU(EAlgorithmU alg){
  return( (alg == ALGORITHM_U_BUNDLE)      ||
          (alg == ALGORITHM_U_LINESEARCH)  ||
          (alg == ALGORITHM_U_TRUSTREGION) ||
          (alg == ALGORITHM_U_LAST)
        );
}

inline EAlgorithmU & operator++(EAlgorithmU &type) {
  return type = static_cast<EAlgorithmU>(type+1);
}

inline EAlgorithmU operator++(EAlgorithmU &type, int) {
  EAlgorithmU oldval = type;
  ++type;
  return oldval;
}

inline EAlgorithmU & operator--(EAlgorithmU &type) {
  return type = static_cast<EAlgorithmU>(type-1);
}

inline EAlgorithmU operator--(EAlgorithmU &type, int) {
  EAlgorithmU oldval = type;
  --type;
  return oldval;
}

inline EAlgorithmU StringToEAlgorithmU(std::string s) {
  s = removeStringFormat(s);
  for ( EAlgorithmU alg = ALGORITHM_U_BUNDLE; alg < ALGORITHM_U_LAST; alg++ ) {
    if ( !s.compare(removeStringFormat(EAlgorithmUToString(alg))) ) {
      return alg;
    }
  }
  return ALGORITHM_U_TRUSTREGION;
}

template<typename Real>
inline Ptr<Algorithm<Real>> AlgorithmFactory(ParameterList &parlist) {
  EAlgorithmU ealg = StringToEAlgorithmU(parlist.sublist("Step").get("Type","Trust Region"));
  switch(ealg) {
    case ALGORITHM_U_BUNDLE:      return makePtr<BundleAlgorithm<Real>>(parlist);
    case ALGORITHM_U_LINESEARCH:  return makePtr<LineSearchAlgorithm<Real>>(parlist);
    case ALGORITHM_U_TRUSTREGION: return makePtr<TrustRegionAlgorithm<Real>>(parlist);
    default:                      return nullPtr;
  }
}

} // namespace TypeU
} // namespace ROL

#endif
