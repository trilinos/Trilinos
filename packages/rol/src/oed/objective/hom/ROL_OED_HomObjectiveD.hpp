// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVED_HPP
#define ROL_OED_HOMOBJECTIVED_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_VectorController.hpp"

#include "ROL_OED_BaseObjective.hpp"
#include "ROL_OED_TraceSampler.hpp"
#include "ROL_OED_MomentOperator.hpp"

/*! \file  OED_D_HomObjective.hpp
    \brief Implements the D-optimal criterion objective for homoscedastic noise.

    Given a covariance operator \f$C(p)\f$ that depends
    on the current design \f$p\f$, this objective function
    implements the logarithm of the determinant of \f$C(p)\f$,
    i.e.,
    \f[
       J(p) = \log(\operatorname{det}(C(p))).
    \f]
*/

namespace ROL::OED::Hom {

template<typename Real>
class ObjectiveD : public BaseObjective<Real> {
private:
  const Ptr<MomentOperator<Real>> M_;
  const Ptr<TraceSampler<Real>> ts_;
  const Ptr<VectorController<Real,int>> stateStore_;
  const Ptr<Vector<Real>> u_, s_, e_, r_;
  const bool storage_;
  const int dim_;
  Ptr<Vector<Real>> g_;
  Real logdet_;
  bool isDetComputed_;

  void solveState(Vector<Real>& u, Vector<Real>& e, const Vector<Real>& z, int i);

public:
  ObjectiveD( const Ptr<MomentOperator<Real>>& M,
              const Ptr<const Vector<Real>>& theta,
              bool storage = true);

  void update(const Vector<Real>& z, UpdateType type, int iter = -1) override;
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

#include "ROL_OED_HomObjectiveD_Def.hpp"

#endif
