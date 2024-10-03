// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_LINEAR_OBJECTIVE_HPP
#define ROL_OED_LINEAR_OBJECTIVE_HPP

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_SampledVector.hpp"
#include "ROL_OED_Factors.hpp"
#include "ROL_OED_TraceSampler.hpp"
#include "ROL_OED_ProfiledClass.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class LinearObjective : public Objective_SimOpt<Real>, public ProfiledClass<Real,std::string> {
private:
  const Ptr<Factors<Real>>           factors_;
  const Ptr<TraceSampler<Real>> traceSampler_;
  const std::string                     type_;
  const Ptr<Vector<Real>>                  g_;

  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

  void computeG(Vector<Real> &g);

public:
  // D Optimality
  LinearObjective();

  // I and R Optimality
  LinearObjective(const Ptr<Factors<Real>> &factors,
                  const std::string &type = "I");

  // C optimality
  LinearObjective(const Ptr<Vector<Real>> &c);

  // A Optimality
  LinearObjective(const Ptr<Vector<Real>>  &theta,
                  const Ptr<TraceSampler<Real>> &traceSampler);

  Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;

  void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;

  void hessVec_11( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;

  void hessVec_12( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;

  void hessVec_21( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;

  void hessVec_22( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;

}; // class LinearObjective

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_LinearObjective_Def.hpp"

#endif
