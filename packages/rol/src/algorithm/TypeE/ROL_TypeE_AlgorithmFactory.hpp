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
inline Ptr<TypeE::Algorithm<Real>> AlgorithmFactory(ParameterList &parlist) {
  EAlgorithmE ealg = StringToEAlgorithmE(parlist.sublist("Step").get("Type","Augmented Lagrangian"));
  switch(ealg) {
    case ALGORITHM_E_AUGMENTEDLAGRANGIAN: return makePtr<TypeE::AugmentedLagrangianAlgorithm<Real>>(parlist);
    case ALGORITHM_E_FLETCHER:            return makePtr<TypeE::FletcherAlgorithm<Real>>(parlist);
    case ALGORITHM_E_COMPOSITESTEP:       return makePtr<TypeE::CompositeStepAlgorithm<Real>>(parlist);
    case ALGORITHM_E_STABILIZEDLCL:       return makePtr<TypeE::StabilizedLCLAlgorithm<Real>>(parlist);
    default:                              return nullPtr;
  }
}

} // namespace TypeE
} // namespace ROL

#endif
