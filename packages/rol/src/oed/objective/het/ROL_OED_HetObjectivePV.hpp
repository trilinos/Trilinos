// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HETOBJECTIVEPV_HPP
#define ROL_OED_HETOBJECTIVEPV_HPP

#include "ROL_PartitionedVector.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_VectorController.hpp"

#include "ROL_OED_BaseObjective.hpp"
#include "ROL_OED_Factors.hpp"
#include "ROL_OED_TraceSampler.hpp"
#include "ROL_OED_MomentOperator.hpp"

namespace ROL::OED::Het {

template<typename Real>
class ObjectivePV : public BaseObjective<Real> {
private:
  const Ptr<MomentOperator<Real>> M0_, M1_;
  const Ptr<Factors<Real>> F_;
  const Ptr<TraceSampler<Real>> ts_;
  const Ptr<VectorController<Real>> stateStore_, MstateStore_, adjointStore_, rhsStore_;
  const Ptr<Vector<Real>> s_, q_, e_, r1_, r2_;
  const bool storage_;
  const unsigned nobs_;
  Ptr<PartitionedVector<Real>> u_, Mu_, p_, r_;
  Ptr<Vector<Real>> g_;

  void solveState(PartitionedVector<Real>& u, PartitionedVector<Real>& Mu, const Vector<Real>& z, const std::vector<Real>& param);
  void solveAdjoint(PartitionedVector<Real>& p, const PartitionedVector<Real>& Mu, const Vector<Real>& z, const std::vector<Real>& param);

public:
  ObjectivePV( const Ptr<MomentOperator<Real>>& M0,
               const Ptr<MomentOperator<Real>>& M1,
               const Ptr<Factors<Real>>& F,
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
    if (F_ != nullPtr) {
      stream << std::string(80,'-') << std::endl;
      F_->summarize(stream,bman);
    }
    if (ts_ != nullPtr) {
      stream << std::string(80,'-') << std::endl;
      ts_->summarize(stream,bman);
    }
  }
  void reset() {
    M0_->reset();
    M1_->reset();
    F_->reset();
    ts_->reset();
  }
};

} // End ROL::OED::Het Namespace

#include "ROL_OED_HetObjectivePV_Def.hpp"

#endif
