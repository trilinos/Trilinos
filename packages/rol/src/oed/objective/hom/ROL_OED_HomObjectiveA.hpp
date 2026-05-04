// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEA_HPP
#define ROL_OED_HOMOBJECTIVEA_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_VectorController.hpp"

#include "ROL_OED_BaseObjective.hpp"
#include "ROL_OED_TraceSampler.hpp"
#include "ROL_OED_MomentOperator.hpp"

namespace ROL::OED::Hom {

template<typename Real>
class ObjectiveA : public BaseObjective<Real> {
private:
  const Ptr<MomentOperator<Real>> M_;
  const Ptr<TraceSampler<Real>> ts_;
  const std::vector<Real> weight_;
  const Ptr<VectorController<Real,int>> stateStore_;
  const Ptr<Vector<Real>> u_, e_, s_;
  Ptr<Vector<Real>> g_;
  const bool storage_;

  void solveState(Vector<Real>& u, Vector<Real>& e, const Vector<Real>& z, int i);

public:
  ObjectiveA( const Ptr<MomentOperator<Real>>& M,
              const Ptr<const Vector<Real>>& theta,
              const Ptr<TraceSampler<Real>>& ts,
              const std::vector<Real>& weight,
              bool storage = true );

  void update( const Vector<Real>& z, UpdateType type, int iter=-1) override;
  Real value( const Vector<Real>& z, Real& tol ) override;
  void gradient( Vector<Real>& g, const Vector<Real>& z, Real& tol ) override;
  void hessVec( Vector<Real>& hv, const Vector<Real>& v, const Vector<Real>& z, Real& tol ) override;
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

#include "ROL_OED_HomObjectiveA_Def.hpp"

#endif
