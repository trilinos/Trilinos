// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEG_ALGORITHMFACTORY_H
#define ROL_TYPEG_ALGORITHMFACTORY_H

#include "ROL_TypeG_AugmentedLagrangianAlgorithm.hpp"
#include "ROL_TypeG_MoreauYosidaAlgorithm.hpp"
#include "ROL_TypeG_InteriorPointAlgorithm.hpp"
#include "ROL_TypeG_StabilizedLCLAlgorithm.hpp"
#include "ROL_Types.hpp"

namespace ROL {
namespace TypeG {

/** \enum   ROL::EAlgorithmG
    \brief  Enumeration of generally constrained algorithm types.

    \arg    ALGORITHM_G_AUGMENTEDLAGRANGIAN describe
    \arg    ALGORITHM_G_MOREAUYOSIDA        describe
    \arg    ALGORITHM_G_INTERIORPOINT       describe
    \arg    ALGORITHM_G_STABILIZEDLCL       describe
 */
enum EAlgorithmG{
  ALGORITHM_G_AUGMENTEDLAGRANGIAN = 0,
  ALGORITHM_G_MOREAUYOSIDA,
  ALGORITHM_G_INTERIORPOINT,
  ALGORITHM_G_STABILIZEDLCL,
  ALGORITHM_G_LAST
};

inline std::string EAlgorithmGToString(EAlgorithmG alg) {
  std::string retString;
  switch(alg) {
    case ALGORITHM_G_AUGMENTEDLAGRANGIAN: retString = "Augmented Lagrangian"; break;
    case ALGORITHM_G_MOREAUYOSIDA:        retString = "Moreau-Yosida";          break;
    case ALGORITHM_G_INTERIORPOINT:       retString = "Interior Point";         break;
    case ALGORITHM_G_STABILIZEDLCL:       retString = "Stabilized LCL";         break;
    case ALGORITHM_G_LAST:                retString = "Last Type (Dummy)";    break;
    default:                              retString = "INVALID EAlgorithmG";
  }
  return retString;
}

/** \brief  Verifies validity of a AlgorithmG enum.
  
    \param  ls  [in]  - enum of the AlgorithmG
    \return 1 if the argument is a valid AlgorithmG; 0 otherwise.
  */
inline int isValidAlgorithmG(EAlgorithmG alg){
  return( (alg == ALGORITHM_G_AUGMENTEDLAGRANGIAN) ||
          (alg == ALGORITHM_G_MOREAUYOSIDA)        ||
          (alg == ALGORITHM_G_INTERIORPOINT)       ||
          (alg == ALGORITHM_G_STABILIZEDLCL)       ||
          (alg == ALGORITHM_G_LAST)
        );
}

inline EAlgorithmG & operator++(EAlgorithmG &type) {
  return type = static_cast<EAlgorithmG>(type+1);
}

inline EAlgorithmG operator++(EAlgorithmG &type, int) {
  EAlgorithmG oldval = type;
  ++type;
  return oldval;
}

inline EAlgorithmG & operator--(EAlgorithmG &type) {
  return type = static_cast<EAlgorithmG>(type-1);
}

inline EAlgorithmG operator--(EAlgorithmG &type, int) {
  EAlgorithmG oldval = type;
  --type;
  return oldval;
}

inline EAlgorithmG StringToEAlgorithmG(std::string s) {
  s = removeStringFormat(s);
  for ( EAlgorithmG alg = ALGORITHM_G_AUGMENTEDLAGRANGIAN; alg < ALGORITHM_G_LAST; alg++ ) {
    if ( !s.compare(removeStringFormat(EAlgorithmGToString(alg))) ) {
      return alg;
    }
  }
  return ALGORITHM_G_AUGMENTEDLAGRANGIAN;
}

template<typename Real>
inline Ptr<TypeG::Algorithm<Real>> AlgorithmFactory(ParameterList &parlist, const Ptr<Secant<Real>> &secant = nullPtr) {
  EAlgorithmG ealg = StringToEAlgorithmG(parlist.sublist("Step").get("Type","Augmented Lagrangian"));
  switch(ealg) {
    case ALGORITHM_G_AUGMENTEDLAGRANGIAN: return makePtr<TypeG::AugmentedLagrangianAlgorithm<Real>>(parlist,secant);
    case ALGORITHM_G_MOREAUYOSIDA:        return makePtr<TypeG::MoreauYosidaAlgorithm<Real>>(parlist,secant);
    case ALGORITHM_G_INTERIORPOINT:       return makePtr<TypeG::InteriorPointAlgorithm<Real>>(parlist,secant);
    case ALGORITHM_G_STABILIZEDLCL:       return makePtr<TypeG::StabilizedLCLAlgorithm<Real>>(parlist,secant);
    default:                              return nullPtr;
  }
}
} // namespace TypeG
} // namespace ROL

#endif
