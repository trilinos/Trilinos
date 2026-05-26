// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEC_HPP
#define ROL_OED_HOMOBJECTIVEC_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"

#include "ROL_OED_BaseObjective.hpp"
#include "ROL_OED_MomentOperator.hpp"

namespace ROL::OED::Hom {

template<typename Real>
class ObjectiveC : public BaseObjective<Real> {
private:
  const Ptr<MomentOperator<Real>> M_;
  const Ptr<const Vector<Real>> c_;
  const Ptr<Vector<Real>> u_, ucache_, r_, s_;
  const bool storage_;
  bool isStateComputed_, isStateCached_;

  void solveState(Vector<Real>& u, const Vector<Real>& z);

public:
  ObjectiveC( const Ptr<MomentOperator<Real>>& M,
               const Ptr<const Vector<Real>>& c,
               bool storage = true);

  void update( const Vector<Real>& z, UpdateType type, int iter=-1 ) override;
  Real value( const Vector<Real> &z, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override;
  void summarize(std::ostream &stream,
           const Ptr<BatchManager<Real>> &bman = nullPtr) const {
    if (M_ != nullPtr) {
      stream << std::string(80,'-') << std::endl;
      M_->summarize(stream,bman);
    }
  }
  void reset() {
    M_->reset();
  }
};

} // END ROL::OED::Hom Namespace

#include "ROL_OED_HomObjectiveC_Def.hpp"

#endif
