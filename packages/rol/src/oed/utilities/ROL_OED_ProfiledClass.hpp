// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_PROFILEDCLASS_HPP
#define ROL_OED_PROFILEDCLASS_HPP

#include "ROL_OED_Timer.hpp"

namespace ROL {
namespace OED {

template<typename Real, typename Key>
class ProfiledClass {
private:
  const Ptr<Timer<Real,Key>> timer_;

public:
  virtual ~ProfiledClass() {}
  ProfiledClass(std::string name = "Default")
    : timer_(makePtr<Timer<Real,Key>>(name)) {}
  void startTimer(std::string name) const { timer_->start(name); }
  void stopTimer(std::string name) const { timer_->stop(name); }
  void rename(std::string name) const { timer_->rename(name); }
  virtual void reset() { timer_->reset(); }
  virtual void summarize(std::ostream &stream,
                   const Ptr<BatchManager<Real>> &bman = nullPtr) const {
    timer_->summarize(stream,bman);
  }
};

} // End OED Namespace
} // End ROL Namespace

#endif
