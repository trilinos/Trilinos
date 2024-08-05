// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_ALGORITHMFACTORY_H
#define ROL_TYPEB_ALGORITHMFACTORY_H

#include "ROL_TypeB_GradientAlgorithm.hpp"
#include "ROL_TypeB_NewtonKrylovAlgorithm.hpp"
#include "ROL_TypeB_LinMoreAlgorithm.hpp"
#include "ROL_TypeB_MoreauYosidaAlgorithm.hpp"
#include "ROL_TypeB_InteriorPointAlgorithm.hpp"
#include "ROL_TypeB_PrimalDualActiveSetAlgorithm.hpp"
#include "ROL_TypeB_KelleySachsAlgorithm.hpp"
#include "ROL_TypeB_SpectralGradientAlgorithm.hpp"
#include "ROL_TypeB_LSecantBAlgorithm.hpp"
#include "ROL_TypeB_QuasiNewtonAlgorithm.hpp"
#include "ROL_TypeB_TrustRegionSPGAlgorithm.hpp"
#include "ROL_TypeB_ColemanLiAlgorithm.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace TypeB {

/** \enum   ROL::EAlgorithmB
    \brief  Enumeration of bound constrained algorithm types.

    \arg    ALGORITHM_B_LINESEARCH          describe
    \arg    ALGORITHM_B_TRUSTREGION         describe
    \arg    ALGORITHM_B_MOREAUYOSIDA        describe
    \arg    ALGORITHM_B_PRIMALDUALACTIVESET describe
    \arg    ALGORITHM_B_INTERIORPOINT       describe
    \arg    ALGORITHM_B_SPECTRALGRADIENT    describe
 */
enum EAlgorithmB{
  ALGORITHM_B_LINESEARCH = 0,
  ALGORITHM_B_TRUSTREGION,
  ALGORITHM_B_MOREAUYOSIDA,
  ALGORITHM_B_PRIMALDUALACTIVESET,
  ALGORITHM_B_INTERIORPOINT,
  ALGORITHM_B_SPECTRALGRADIENT,
  ALGORITHM_B_LAST
};

inline std::string EAlgorithmBToString(EAlgorithmB alg) {
  std::string retString;
  switch(alg) {
    case ALGORITHM_B_LINESEARCH:          retString = "Line Search";            break;
    case ALGORITHM_B_TRUSTREGION:         retString = "Trust Region";           break;
    case ALGORITHM_B_MOREAUYOSIDA:        retString = "Moreau-Yosida";          break;
    case ALGORITHM_B_PRIMALDUALACTIVESET: retString = "Primal Dual Active Set"; break;
    case ALGORITHM_B_INTERIORPOINT:       retString = "Interior Point";         break;
    case ALGORITHM_B_SPECTRALGRADIENT:    retString = "Spectral Gradient";      break;
    case ALGORITHM_B_LAST:                retString = "Last Type (Dummy)";      break;
    default:                              retString = "INVALID EAlgorithmB";
  }
  return retString;
}

/** \brief  Verifies validity of a AlgorithmB enum.
  
    \param  ls  [in]  - enum of the AlgorithmB
    \return 1 if the argument is a valid AlgorithmB; 0 otherwise.
  */
inline int isValidAlgorithmB(EAlgorithmB alg){
  return( (alg == ALGORITHM_B_LINESEARCH)          ||
          (alg == ALGORITHM_B_TRUSTREGION)         ||
          (alg == ALGORITHM_B_MOREAUYOSIDA)        ||
          (alg == ALGORITHM_B_PRIMALDUALACTIVESET) ||
          (alg == ALGORITHM_B_INTERIORPOINT)       ||
          (alg == ALGORITHM_B_SPECTRALGRADIENT)    ||
          (alg == ALGORITHM_B_LAST)
        );
}

inline EAlgorithmB & operator++(EAlgorithmB &type) {
  return type = static_cast<EAlgorithmB>(type+1);
}

inline EAlgorithmB operator++(EAlgorithmB &type, int) {
  EAlgorithmB oldval = type;
  ++type;
  return oldval;
}

inline EAlgorithmB & operator--(EAlgorithmB &type) {
  return type = static_cast<EAlgorithmB>(type-1);
}

inline EAlgorithmB operator--(EAlgorithmB &type, int) {
  EAlgorithmB oldval = type;
  --type;
  return oldval;
}

inline EAlgorithmB StringToEAlgorithmB(std::string s) {
  s = removeStringFormat(s);
  for ( EAlgorithmB alg = ALGORITHM_B_LINESEARCH; alg < ALGORITHM_B_LAST; alg++ ) {
    if ( !s.compare(removeStringFormat(EAlgorithmBToString(alg))) ) {
      return alg;
    }
  }
  return ALGORITHM_B_TRUSTREGION;
}

template<typename Real>
inline Ptr<Algorithm<Real>> AlgorithmFactory(ParameterList &parlist, const Ptr<Secant<Real>> &secant = nullPtr) {
  EAlgorithmB ealg = StringToEAlgorithmB(parlist.sublist("Step").get("Type","Trust Region"));
  switch(ealg) {
    case ALGORITHM_B_LINESEARCH:
    {
      std::string desc
        = parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").get("Type","Newton-Krylov");
      if (desc=="Newton-Krylov" || desc=="Newton")
        return makePtr<NewtonKrylovAlgorithm<Real>>(parlist,secant);
      else if (desc=="Quasi-Newton Method" || desc=="Quasi-Newton") {
        std::string method
          = parlist.sublist("Step").sublist("Line Search").sublist("Quasi-Newton").get("Method","L-Secant-B");
        if (method == "L-Secant-B")
          return makePtr<LSecantBAlgorithm<Real>>(parlist,secant);    // Similar to L-BFGS-B
        else
          return makePtr<QuasiNewtonAlgorithm<Real>>(parlist,secant); // PQN
      }
      else {
        return makePtr<GradientAlgorithm<Real>>(parlist);
      }
    }
    case ALGORITHM_B_TRUSTREGION:
    {
      std::string trmod
        = parlist.sublist("Step").sublist("Trust Region").get("Subproblem Model","Lin-More");
      if (trmod=="Kelley-Sachs")
        return makePtr<KelleySachsAlgorithm<Real>>(parlist,secant);
      else if (trmod=="SPG")
        return makePtr<TrustRegionSPGAlgorithm<Real>>(parlist,secant);
      else if (trmod=="Coleman-Li")
        return makePtr<ColemanLiAlgorithm<Real>>(parlist,secant);
      else
        return makePtr<LinMoreAlgorithm<Real>>(parlist,secant);
    }
    case ALGORITHM_B_MOREAUYOSIDA:        return makePtr<MoreauYosidaAlgorithm<Real>>(parlist,secant);
    case ALGORITHM_B_PRIMALDUALACTIVESET: return makePtr<PrimalDualActiveSetAlgorithm<Real>>(parlist,secant);
    case ALGORITHM_B_INTERIORPOINT:       return makePtr<InteriorPointAlgorithm<Real>>(parlist,secant);
    case ALGORITHM_B_SPECTRALGRADIENT:    return makePtr<SpectralGradientAlgorithm<Real>>(parlist);
    default:                              return nullPtr;
  }
}
} // namespace TypeB
} // namespace ROL

#endif
