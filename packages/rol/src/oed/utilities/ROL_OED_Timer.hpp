// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_TIMER_HPP
#define ROL_OED_TIMER_HPP

#include <vector>
#include <ctime>
#include <iostream>
#include "ROL_Types.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_BatchManager.hpp"

namespace ROL {
namespace OED {

template<typename Real, typename Key>
class Timer {
private:
  mutable std::map<Key,std::clock_t> timer_;
  mutable std::map<Key,Real>         time_;
  mutable std::map<Key,int>          count_;
  std::string                        name_;

public:
  Timer(std::string name = "Default");

  void rename(std::string name);
  void reset();
  void start(Key i);
  void stop(Key i);
  void summarize(std::ostream &stream,
           const Ptr<BatchManager<Real>> &bman = nullPtr) const;
};

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_Timer_Def.hpp"

#endif
