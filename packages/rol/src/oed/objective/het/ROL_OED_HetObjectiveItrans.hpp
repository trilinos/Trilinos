// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HETOBJECTIVEITRANS_HPP
#define ROL_OED_HETOBJECTIVEITRANS_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_VectorController.hpp"

#include "ROL_OED_BaseObjective.hpp"
#include "ROL_OED_Factors.hpp"
#include "ROL_OED_TraceSampler.hpp"
#include "ROL_OED_MomentOperator.hpp"

namespace ROL::OED::Het {

template<typename Real>
class ObjectiveItrans : public BaseObjective<Real> {
private:
  const Ptr<MomentOperator<Real>> M0_, M1_;
  const Ptr<TraceSampler<Real>> ts_;
  const Ptr<VectorController<Real,unsigned>> state0Store_, state1Store_, Mstate0Store_, Mstate1Store_, adjoint0Store_, adjoint1Store_;
  const Ptr<Vector<Real>> u0_, u1_, Mu0_, Mu1_, p0_, p1_, s_, r_, r1_;
  const bool storage_;
  const unsigned dim_;
  std::vector<Ptr<Vector<Real>>> b_;
  Ptr<Vector<Real>> g_;

  void solveState(Vector<Real>& u0, Vector<Real>& u1, Vector<Real>& Mu0, Vector<Real>& Mu1, const Vector<Real>& z, unsigned i);
  void solveAdjoint(Vector<Real>& p0, Vector<Real>& p1, const Vector<Real>& Mu0, const Vector<Real>& Mu1, const Vector<Real>& z, unsigned i);

public:
  ObjectiveItrans( const Ptr<MomentOperator<Real>>& M0,
                   const Ptr<MomentOperator<Real>>& M1,
                   const Ptr<Factors<Real>>& F,
                   const Ptr<SampleGenerator<Real>>& sampler,
                   const Ptr<TraceSampler<Real>>& ts,
                   const std::vector<Real>& wt,
                   bool storage = true);

  void update( const Vector<Real>& z, UpdateType type, int iter=-1 ) override;
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
    if (ts_ != nullPtr) {
      stream << std::string(80,'-') << std::endl;
      ts_->summarize(stream,bman);
    }
  }
  void reset() {
    M0_->reset();
    M1_->reset();
    ts_->reset();
  }
};

} // END ROL::OED::Het Namespace

#include "ROL_OED_HetObjectiveItrans_Def.hpp"

#endif
