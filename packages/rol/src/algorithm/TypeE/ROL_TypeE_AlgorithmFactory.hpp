// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEE_ALGORITHMFACTORY_H
#define ROL_TYPEE_ALGORITHMFACTORY_H

#include "ROL_TypeE_AugmentedLagrangianAlgorithm.hpp"
#include "ROL_TypeE_FletcherAlgorithm.hpp"
#include "ROL_TypeE_CompositeStepAlgorithm.hpp"
#include "ROL_TypeE_StabilizedLCLAlgorithm.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace TypeE {

/** \enum   ROL::TypeE::AlgorithmFactory
    \brief  Enumeration of equality constrained algorithm types.
 */
enum EAlgorithmE{
  ALGORITHM_E_AUGMENTEDLAGRANGIAN = 0,
  ALGORITHM_E_FLETCHER,
  ALGORITHM_E_COMPOSITESTEP,
  ALGORITHM_E_STABILIZEDLCL,
  ALGORITHM_E_LAST
};

inline std::string EAlgorithmEToString(EAlgorithmE alg) {
  std::string retString;
  switch(alg) {
    case ALGORITHM_E_AUGMENTEDLAGRANGIAN: retString = "Augmented Lagrangian"; break;
    case ALGORITHM_E_FLETCHER:            retString = "Fletcher";             break;
    case ALGORITHM_E_COMPOSITESTEP:       retString = "Composite Step";       break;
    case ALGORITHM_E_STABILIZEDLCL:       retString = "Stabilized LCL";       break;
    case ALGORITHM_E_LAST:                retString = "Last Type (Dummy)";    break;
    default:                              retString = "INVALID EAlgorithmE";
  }
  return retString;
}

/** \brief  Verifies validity of a AlgorithmE enum.
  
    \param  ls  [in]  - enum of the AlgorithmE
    \return 1 if the argument is a valid AlgorithmE; 0 otherwise.
  */
inline int isValidAlgorithmE(EAlgorithmE alg){
  return( (alg == ALGORITHM_E_AUGMENTEDLAGRANGIAN) ||
          (alg == ALGORITHM_E_FLETCHER)            ||
          (alg == ALGORITHM_E_COMPOSITESTEP)       ||
          (alg == ALGORITHM_E_STABILIZEDLCL)       ||
          (alg == ALGORITHM_E_LAST)
        );
}

inline EAlgorithmE & operator++(EAlgorithmE &type) {
  return type = static_cast<EAlgorithmE>(type+1);
}

inline EAlgorithmE operator++(EAlgorithmE &type, int) {
  EAlgorithmE oldval = type;
  ++type;
  return oldval;
}

inline EAlgorithmE & operator--(EAlgorithmE &type) {
  return type = static_cast<EAlgorithmE>(type-1);
}

inline EAlgorithmE operator--(EAlgorithmE &type, int) {
  EAlgorithmE oldval = type;
  --type;
  return oldval;
}

inline EAlgorithmE StringToEAlgorithmE(std::string s) {
  s = removeStringFormat(s);
  for ( EAlgorithmE alg = ALGORITHM_E_AUGMENTEDLAGRANGIAN; alg < ALGORITHM_E_LAST; alg++ ) {
    if ( !s.compare(removeStringFormat(EAlgorithmEToString(alg))) ) {
      return alg;
    }
  }
  return ALGORITHM_E_AUGMENTEDLAGRANGIAN;
}

template<typename Real>
inline Ptr<TypeE::Algorithm<Real>> AlgorithmFactory(ParameterList &parlist, const Ptr<Secant<Real>> &secant = nullPtr) {
  EAlgorithmE ealg = StringToEAlgorithmE(parlist.sublist("Step").get("Type","Augmented Lagrangian"));
  switch(ealg) {
    case ALGORITHM_E_AUGMENTEDLAGRANGIAN: return makePtr<TypeE::AugmentedLagrangianAlgorithm<Real>>(parlist,secant);
    case ALGORITHM_E_FLETCHER:            return makePtr<TypeE::FletcherAlgorithm<Real>>(parlist,secant);
    case ALGORITHM_E_COMPOSITESTEP:       return makePtr<TypeE::CompositeStepAlgorithm<Real>>(parlist);
    case ALGORITHM_E_STABILIZEDLCL:       return makePtr<TypeE::StabilizedLCLAlgorithm<Real>>(parlist,secant);
    default:                              return nullPtr;
  }
}

} // namespace TypeE
} // namespace ROL

#endif
