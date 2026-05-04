// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEI_HPP
#define ROL_OED_HOMOBJECTIVEI_HPP

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
class ObjectiveI : public BaseObjective<Real> {
private:
  Ptr<BaseObjective<Real>> obj_;

public:
  ObjectiveI( const Ptr<MomentOperator<Real>>& M,
              const Ptr<Factors<Real>>& F,
              const Ptr<SampleGenerator<Real>>& sampler,
              const Ptr<TraceSampler<Real>>& ts,
              const std::vector<Real>& wt,
              bool storage = true);

  ObjectiveI( const Ptr<MomentOperator<Real>>& M,
              const Ptr<Factors<Real>>& F,
              const Ptr<SampleGenerator<Real>>& sampler,
              bool storage = true,
              bool useBasic = false);

  void update( const Vector<Real>& z, UpdateType type, int iter=-1 ) override;
  Real value( const Vector<Real> &z, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override;
  void summarize(std::ostream &stream,
           const Ptr<BatchManager<Real>> &bman = nullPtr) const {
    obj_->summarize(stream,bman);
  }
  void reset() {
    obj_->reset();
  }
};

} // END ROL::OED::Hom Namespace

#include "ROL_OED_HomObjectiveI_Def.hpp"

#endif
