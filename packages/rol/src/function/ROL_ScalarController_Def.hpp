// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALARCONTROLLER_DEF_H
#define ROL_SCALARCONTROLLER_DEF_H

namespace ROL {

template <class Real, class Key>
ScalarController<Real,Key>::ScalarController(void) : VectorController<Real,Key>() {}

template <class Real, class Key>
bool ScalarController<Real,Key>::get(Real &x, const Key &param) {
  SingletonVector<Real> xv(Real(0));
  bool flag = VectorController<Real,Key>::get(xv,param);
  x = xv.getValue();
  return flag;
}

template <class Real, class Key>
void ScalarController<Real,Key>::set(Real x, const Key &param) {
  SingletonVector<Real> xv(x);
  VectorController<Real,Key>::set(xv,param);
}

} // namespace ROL

#endif
