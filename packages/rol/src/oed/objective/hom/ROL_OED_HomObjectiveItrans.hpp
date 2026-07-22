// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEITRANS_HPP
#define ROL_OED_HOMOBJECTIVEITRANS_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_VectorController.hpp"

#include "ROL_OED_BaseObjective.hpp"
#include "ROL_OED_Factors.hpp"
#include "ROL_OED_TraceSampler.hpp"
#include "ROL_OED_MomentOperator.hpp"

namespace ROL::OED::Hom {

template<typename Real>
class ObjectiveItrans : public BaseObjective<Real> {
private:
  const Ptr<MomentOperator<Real>> M_;
  const Ptr<TraceSampler<Real>> ts_;
  const Ptr<VectorController<Real,int>> stateStore_, adjointStore_;
  const Ptr<Vector<Real>> u_, p_, s_, r_;
  const bool storage_;
  const int dim_;
  std::vector<Ptr<Vector<Real>>> b_;
  Ptr<Vector<Real>> g_;

  void solveState(Vector<Real>& u, const Vector<Real>& z, int i);
  void solveAdjoint(Vector<Real>& p, const Vector<Real>& z, int i);

public:
  ObjectiveItrans( const Ptr<MomentOperator<Real>>& M,
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
    if (M_ != nullPtr) {
      stream << std::string(80,'-') << std::endl;
      M_->summarize(stream,bman);
    }
    if (ts_ != nullPtr) {
      stream << std::string(80,'-') << std::endl;
      ts_->summarize(stream,bman);
    }
  }
  void reset() {
    M_->reset();
    ts_->reset();
  }
};

} // END ROL::OED::Hom Namespace

#include "ROL_OED_HomObjectiveItrans_Def.hpp"

#endif
