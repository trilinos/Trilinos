// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HETOBJECTIVEC_HPP
#define ROL_OED_HETOBJECTIVEC_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"

#include "ROL_OED_BaseObjective.hpp"
#include "ROL_OED_MomentOperator.hpp"

namespace ROL::OED::Het {

template<typename Real>
class ObjectiveC : public BaseObjective<Real> {
private:
  const Ptr<MomentOperator<Real>> M0_, M1_;
  const Ptr<const Vector<Real>> c_;
  const Ptr<Vector<Real>> u_, ucache_, Mu_, Mucache_, p_, pcache_, r1_, r2_, s_, q_;
  const bool storage_;
  bool isStateComputed_, isStateCached_, isAdjointComputed_, isAdjointCached_;
  Ptr<Vector<Real>> g_;

  void solveState(Vector<Real>& u, Vector<Real>& Mu, const Vector<Real>& z);
  void solveAdjoint(Vector<Real>& p, const Vector<Real>& Mu, const Vector<Real>& z);

public:
  ObjectiveC( const Ptr<MomentOperator<Real>>& M0,
              const Ptr<MomentOperator<Real>>& M1,
              const Ptr<const Vector<Real>>& c,
              bool storage = true);

  void update( const Vector<Real>& z, UpdateType type, int iter=-1) override;
  Real value( const Vector<Real> &z, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override;
  void summarize(std::ostream &stream,
           const Ptr<BatchManager<Real>> &bman = nullPtr) const {
    if (M0_ != nullPtr) {
      stream << std::string(80,'-') << std::endl;
      M0_->summarize(stream,bman);
    }
    if (M1_ != nullPtr) {
      stream << std::string(80,'-') << std::endl;
      M1_->summarize(stream,bman);
    }
  }
  void reset() {
    M0_->reset();
    M1_->reset();
  }
};

} // End ROL::OED::Het Namespace

#include "ROL_OED_HetObjectiveC_Def.hpp"

#endif
