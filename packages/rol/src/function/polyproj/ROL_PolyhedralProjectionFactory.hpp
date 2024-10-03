// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_POLYHEDRALPROJECTIONFACTORY_H
#define ROL_POLYHEDRALPROJECTIONFACTORY_H

#include "ROL_DaiFletcherProjection.hpp"
#include "ROL_DykstraProjection.hpp"
#include "ROL_DouglasRachfordProjection.hpp"
#include "ROL_SemismoothNewtonProjection.hpp"
#include "ROL_RiddersProjection.hpp"
#include "ROL_BrentsProjection.hpp"

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
  PPA_DOUGLASRACHFORD,
  PPA_NEWTON,
  PPA_RIDDERS,
  PPA_BRENTS,
  PPA_LAST
};

inline std::string EPolyProjAlgoToString(EPolyProjAlgo alg) {
  std::string retString;
  switch(alg) {
    case PPA_DAIFLETCHER:     retString = "Dai-Fletcher";      break;
    case PPA_DYKSTRA:         retString = "Dysktra";           break;
    case PPA_DOUGLASRACHFORD: retString = "Douglas-Rachford";  break;
    case PPA_NEWTON:          retString = "Semismooth Newton"; break;
    case PPA_RIDDERS:         retString = "Ridders";           break;
    case PPA_BRENTS:          retString = "Brents";            break;
    case PPA_LAST:            retString = "Last Type (Dummy)"; break;
    default:                  retString = "INVALID EPolyProjAlgo";
  }
  return retString;
}

/** \brief  Verifies validity of a PolyProjAlgo enum.
  
    \param  ls  [in]  - enum of the PolyProjAlgo
    \return 1 if the argument is a valid PolyProjAlgo; 0 otherwise.
  */
inline int isValidPolyProjAlgo(EPolyProjAlgo alg){
  return( (alg == PPA_DAIFLETCHER)     ||
          (alg == PPA_DYKSTRA)         ||
          (alg == PPA_DOUGLASRACHFORD) ||
          (alg == PPA_NEWTON)          ||
          (alg == PPA_RIDDERS)         ||
          (alg == PPA_BRENTS)          ||
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
    case PPA_DAIFLETCHER:     return makePtr<DaiFletcherProjection<Real>>(xprim,xdual,bnd,con,mul,res,list);      break;
    case PPA_DYKSTRA:         return makePtr<DykstraProjection<Real>>(xprim,xdual,bnd,con,mul,res,list);          break;
    case PPA_DOUGLASRACHFORD: return makePtr<DouglasRachfordProjection<Real>>(xprim,xdual,bnd,con,mul,res,list);  break;
    case PPA_NEWTON:          return makePtr<SemismoothNewtonProjection<Real>>(xprim,xdual,bnd,con,mul,res,list); break;
    case PPA_RIDDERS:         return makePtr<RiddersProjection<Real>>(xprim,xdual,bnd,con,mul,res,list);          break;
    case PPA_BRENTS:          return makePtr<BrentsProjection<Real>>(xprim,xdual,bnd,con,mul,res,list);           break;
    default:                  return nullPtr;
  }
}
} // namespace ROL

#endif
