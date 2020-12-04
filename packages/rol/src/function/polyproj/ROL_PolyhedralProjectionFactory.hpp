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


#ifndef ROL_POLYHEDRALPROJECTIONFACTORY_H
#define ROL_POLYHEDRALPROJECTIONFACTORY_H

#include "ROL_DaiFletcherProjection.hpp"
#include "ROL_DykstraProjection.hpp"
#include "ROL_SemismoothNewtonProjection.hpp"
#include "ROL_RiddersProjection.hpp"

namespace ROL {

/** \enum   ROL::EPolyProjAlgo
    \brief  Enumeration of polyhedral projecdtion algorithm types.

    \arg    PPA_DAIFLETCHER describe
    \arg    PPA_DYKSTRA     describe
    \arg    PPA_NEWTON      describe
    \arg    PPA_RIDDERS     describe
 */
enum EPolyProjAlgo{
  PPA_DAIFLETCHER = 0,
  PPA_DYKSTRA,
  PPA_NEWTON,
  PPA_RIDDERS,
  PPA_LAST
};

inline std::string EPolyProjAlgoToString(EPolyProjAlgo alg) {
  std::string retString;
  switch(alg) {
    case PPA_DAIFLETCHER: retString = "Dai-Fletcher";      break;
    case PPA_DYKSTRA:     retString = "Dysktra";           break;
    case PPA_NEWTON:      retString = "Semismooth Newton"; break;
    case PPA_RIDDERS:     retString = "Ridders";           break;
    case PPA_LAST:        retString = "Last Type (Dummy)"; break;
    default:              retString = "INVALID EPolyProjAlgo";
  }
  return retString;
}

/** \brief  Verifies validity of a PolyProjAlgo enum.
  
    \param  ls  [in]  - enum of the PolyProjAlgo
    \return 1 if the argument is a valid PolyProjAlgo; 0 otherwise.
  */
inline int isValidPolyProjAlgo(EPolyProjAlgo alg){
  return( (alg == PPA_DAIFLETCHER) ||
          (alg == PPA_DYKSTRA)     ||
          (alg == PPA_NEWTON)      ||
          (alg == PPA_RIDDERS)     ||
          (alg == PPA_LAST)
        );
}

inline EPolyProjAlgo & operator++(EPolyProjAlgo &type) {
  return type = static_cast<EPolyProjAlgo>(type+1);
}

inline EPolyProjAlgo operator++(EPolyProjAlgo &type, int) {
  EPolyProjAlgo oldval = type;
  ++type;
  return oldval;
}

inline EPolyProjAlgo & operator--(EPolyProjAlgo &type) {
  return type = static_cast<EPolyProjAlgo>(type-1);
}

inline EPolyProjAlgo operator--(EPolyProjAlgo &type, int) {
  EPolyProjAlgo oldval = type;
  --type;
  return oldval;
}

inline EPolyProjAlgo StringToEPolyProjAlgo(std::string s) {
  s = removeStringFormat(s);
  for ( EPolyProjAlgo alg = PPA_DAIFLETCHER; alg < PPA_LAST; alg++ ) {
    if ( !s.compare(removeStringFormat(EPolyProjAlgoToString(alg))) ) {
      return alg;
    }
  }
  return PPA_DYKSTRA;
}

template<typename Real>
inline Ptr<PolyhedralProjection<Real>> PolyhedralProjectionFactory(const Vector<Real>               &xprim,
                                                                   const Vector<Real>               &xdual,
                                                                   const Ptr<BoundConstraint<Real>> &bnd,
                                                                   const Ptr<Constraint<Real>>      &con,
                                                                   const Vector<Real>               &mul,
                                                                   const Vector<Real>               &res,
                                                                   ParameterList                    &list) {
  EPolyProjAlgo ealg = StringToEPolyProjAlgo(list.sublist("General").sublist("Polyhedral Projection").get("Type","Dykstra"));
  switch(ealg) {
    case PPA_DAIFLETCHER: return makePtr<DaiFletcherProjection<Real>>(xprim,xdual,bnd,con,mul,res,list);      break;
    case PPA_DYKSTRA:     return makePtr<DykstraProjection<Real>>(xprim,xdual,bnd,con,mul,res,list);          break;
    case PPA_NEWTON:      return makePtr<SemismoothNewtonProjection<Real>>(xprim,xdual,bnd,con,mul,res,list); break;
    case PPA_RIDDERS:     return makePtr<RiddersProjection<Real>>(xprim,xdual,bnd,con,mul,res,list);          break;
    default:              return nullPtr;
  }
}
} // namespace ROL

#endif
