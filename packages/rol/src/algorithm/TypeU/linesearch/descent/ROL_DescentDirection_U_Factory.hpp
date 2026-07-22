// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DESCENTDIRECTION_U_FACTORY_H
#define ROL_DESCENTDIRECTION_U_FACTORY_H

#include "ROL_Gradient_U.hpp"
#include "ROL_QuasiNewton_U.hpp"
#include "ROL_NonlinearCG_U.hpp"
#include "ROL_Newton_U.hpp"
#include "ROL_NewtonKrylov_U.hpp"
#include "ROL_LineSearch_U_Types.hpp"

namespace ROL {
template<typename Real>
inline Ptr<DescentDirection_U<Real>> DescentDirectionUFactory(ParameterList &parlist, const Ptr<Secant<Real>> &secant = nullPtr) {
  EDescentU edesc = StringToEDescentU(
    parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").get("Type","Quasi-Newton Method"));
  switch(edesc) {
    case DESCENT_U_STEEPEST:     return makePtr<Gradient_U<Real>>();
    case DESCENT_U_NONLINEARCG:  return makePtr<NonlinearCG_U<Real>>(parlist);
    case DESCENT_U_SECANT:       return makePtr<QuasiNewton_U<Real>>(parlist, secant);
    case DESCENT_U_NEWTON:       return makePtr<Newton_U<Real>>();
    case DESCENT_U_NEWTONKRYLOV:
    {
      Ptr<Krylov<Real>> krylov = nullPtr;
      return makePtr<NewtonKrylov_U<Real>>(parlist,krylov,secant);
    }
    default:                     return nullPtr;
  }
}
} // namespace ROL

#endif
