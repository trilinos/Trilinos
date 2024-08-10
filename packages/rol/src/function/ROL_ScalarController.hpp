// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALARCONTROLLER_H
#define ROL_SCALARCONTROLLER_H

#include "ROL_VectorController.hpp"
#include "ROL_SingletonVector.hpp"

namespace ROL {

template <class Real, class Key=std::vector<Real>>
class ScalarController : public VectorController<Real,Key> {
public:
  ScalarController(void);

  bool get(Real &x, const Key &param);

  void set(Real x, const Key &param);

}; // class ScalarController

} // namespace ROL

#include "ROL_ScalarController_Def.hpp"

#endif
