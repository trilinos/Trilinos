// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEP_ALGORITHMFACTORY_H
#define ROL_TYPEP_ALGORITHMFACTORY_H

#include "ROL_TypeP_ProxGradientAlgorithm.hpp"
#include "ROL_TypeP_SpectralGradientAlgorithm.hpp"
#include "ROL_TypeP_iPianoAlgorithm.hpp"
#include "ROL_TypeP_QuasiNewtonAlgorithm.hpp"
#include "ROL_TypeP_TrustRegionAlgorithm.hpp"
#include "ROL_TypeP_InexactNewtonAlgorithm.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace TypeP {

/** \enum   ROL::EAlgorithmP
    \brief  Enumeration of bound constrained algorithm types.

    \arg    ALGORITHM_P_LINESEARCH          describe
    \arg    ALGORITHM_P_TRUSTREGION         describe
    \arg    ALGORITHM_P_SPECTRALGRADIENT    describe
    \arg    ALGORITHM_P_IPIANO              describe
 */
enum EAlgorithmP{
  ALGORITHM_P_LINESEARCH = 0,
  ALGORITHM_P_TRUSTREGION,
  ALGORITHM_P_SPECTRALGRADIENT,
  ALGORITHM_P_IPIANO,
  ALGORITHM_P_LAST
};

inline std::string EAlgorithmPToString(EAlgorithmP alg) {
  std::string retString;
  switch(alg) {
    case ALGORITHM_P_LINESEARCH:          retString = "Line Search";            break;
    case ALGORITHM_P_TRUSTREGION:         retString = "Trust Region";           break;
    case ALGORITHM_P_SPECTRALGRADIENT:    retString = "Spectral Gradient";      break;
    case ALGORITHM_P_IPIANO:              retString = "iPiano";                 break;
    case ALGORITHM_P_LAST:                retString = "Last Type (Dummy)";      break;
    default:                              retString = "INVALID EAlgorithmP";
  }
  return retString;
}

/** \brief  Verifies validity of a AlgorithmP enum.
  
    \param  ls  [in]  - enum of the AlgorithmP
    \return 1 if the argument is a valid AlgorithmP; 0 otherwise.
  */
inline int isValidAlgorithmP(EAlgorithmP alg){
  return( (alg == ALGORITHM_P_LINESEARCH)          ||
          (alg == ALGORITHM_P_TRUSTREGION)         ||
          (alg == ALGORITHM_P_SPECTRALGRADIENT)    ||
          (alg == ALGORITHM_P_IPIANO)              ||
          (alg == ALGORITHM_P_LAST)
        );
}

inline EAlgorithmP & operator++(EAlgorithmP &type) {
  return type = static_cast<EAlgorithmP>(type+1);
}

inline EAlgorithmP operator++(EAlgorithmP &type, int) {
  EAlgorithmP oldval = type;
  ++type;
  return oldval;
}

inline EAlgorithmP & operator--(EAlgorithmP &type) {
  return type = static_cast<EAlgorithmP>(type-1);
}

inline EAlgorithmP operator--(EAlgorithmP &type, int) {
  EAlgorithmP oldval = type;
  --type;
  return oldval;
}

inline EAlgorithmP StringToEAlgorithmP(std::string s) {
  s = removeStringFormat(s);
  for ( EAlgorithmP alg = ALGORITHM_P_LINESEARCH; alg < ALGORITHM_P_LAST; alg++ ) {
    if ( !s.compare(removeStringFormat(EAlgorithmPToString(alg))) ) {
      return alg;
    }
  }
  return ALGORITHM_P_TRUSTREGION;
}

template<typename Real>
inline Ptr<Algorithm<Real>> AlgorithmFactory(ParameterList &parlist) {
  EAlgorithmP ealg = StringToEAlgorithmP(parlist.sublist("Step").get("Type","Trust Region"));
  switch(ealg) {
    case ALGORITHM_P_LINESEARCH:
    {
      std::string desc
        = parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").get("Type","Newton-Krylov");
      if (desc=="Newton-Krylov" || desc=="Newton")
        return makePtr<InexactNewtonAlgorithm<Real>>(parlist);
      else if (desc=="Quasi-Newton Method" || desc = "Quasi-Newton")
        return makePtr<QuasiNewtonAlgorithm<Real>>(parlist);
      else
        return makePtr<ProxGradientAlgorithm<Real>>(parlist);
    }
    case ALGORITHM_P_TRUSTREGION:         return makePtr<TrustRegionAlgorithm<Real>>(parlist);
    case ALGORITHM_P_SPECTRALGRADIENT:    return makePtr<SpectralGradientAlgorithm<Real>>(parlist);
    case ALGORITHM_P_IPIANO:              return makePtr<iPianoAlgorithm<Real>>(parlist);
    default:                              return nullPtr;
  }
}
} // namespace TypeP
} // namespace ROL

#endif
