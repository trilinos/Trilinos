// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEPV_HPP
#define ROL_OED_HOMOBJECTIVEPV_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_VectorController.hpp"

#include "ROL_OED_BaseObjective.hpp"
#include "ROL_OED_Factors.hpp"
#include "ROL_OED_TraceSampler.hpp"
#include "ROL_OED_MomentOperator.hpp"

namespace ROL::OED::Hom {

template<typename Real>
class ObjectivePV : public BaseObjective<Real> {
private:
  const Ptr<MomentOperator<Real>> M_;
  const Ptr<Factors<Real>> F_;
  const Ptr<TraceSampler<Real>> ts_;
  const Ptr<VectorController<Real>> stateStore_, rhsStore_;
  const Ptr<Vector<Real>> s_, e_, d_;
  const bool storage_;
  const unsigned nobs_;
  Ptr<PartitionedVector<Real>> u_, r_;
  Ptr<Vector<Real>> g_;

  void solveState(PartitionedVector<Real>& u, PartitionedVector<Real>& r, const Vector<Real>& z, const std::vector<Real>& param);

public:
  ObjectivePV( const Ptr<MomentOperator<Real>>& M,
               const Ptr<Factors<Real>>& F,
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
    M_->reset();
    F_->reset();
    ts_->reset();
  }
};

} // END ROL::OED::Hom Namespace

#include "ROL_OED_HomObjectivePV_Def.hpp"

#endif
