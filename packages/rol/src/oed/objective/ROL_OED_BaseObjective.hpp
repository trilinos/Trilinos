// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef ROL_OED_BASEOBJECTIVE_HPP
#define ROL_OED_BASEOBJECTIVE_HPP

#include <ostream>
#include "ROL_Objective.hpp"
#include "ROL_BatchManager.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class BaseObjective : public Objective<Real> {
public:
  BaseObjective() {}
  virtual void summarize(std::ostream &stream,
                   const Ptr<BatchManager<Real>> &bman = nullPtr) const {}
  virtual void reset() {}

};

} // End OED Namespace
} // End ROL Namespace

#endif
