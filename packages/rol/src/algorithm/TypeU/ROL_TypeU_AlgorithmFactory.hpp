// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
inline Ptr<Algorithm<Real>> AlgorithmFactory(ParameterList &parlist, const Ptr<Secant<Real>> &secant = nullPtr) {
  EAlgorithmU ealg = StringToEAlgorithmU(parlist.sublist("Step").get("Type","Trust Region"));
  switch(ealg) {
    case ALGORITHM_U_BUNDLE:      return makePtr<BundleAlgorithm<Real>>(parlist);
    case ALGORITHM_U_LINESEARCH:  return makePtr<LineSearchAlgorithm<Real>>(parlist,secant);
    case ALGORITHM_U_TRUSTREGION: return makePtr<TrustRegionAlgorithm<Real>>(parlist,secant);
    default:                      return nullPtr;
  }
}

} // namespace TypeU
} // namespace ROL

#endif
