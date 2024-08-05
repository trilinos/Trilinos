// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_TRACESAMPLER_HPP
#define ROL_OED_TRACESAMPLER_HPP

#include "ROL_Ptr.hpp"
#include "ROL_SampledVector.hpp"
#include "ROL_Vector.hpp"
#include "ROL_OED_ProfiledClass.hpp"
#include <vector>

namespace ROL {
namespace OED {

template<typename Real>
class TraceSampler : public ProfiledClass<Real,std::string> {
private:
  Ptr<SampledVector<Real>>     g_storage_;

protected:
  void setInStorage(const Vector<Real> &g, const std::vector<Real> &param);
  bool getFromStorage(Vector<Real> &g, const std::vector<Real> &param) const;

  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

public:
  virtual ~TraceSampler() {}
  TraceSampler();
  TraceSampler(const Ptr<Vector<Real>> &theta);

  virtual void generate(Vector<Real> &g) const;
  virtual void get(Vector<Real> &F, const std::vector<Real> &param);

}; // class TraceSampler

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_TraceSampler_Def.hpp"

#endif
