// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_QUADRATIC_OBJECTIVE_HPP
#define ROL_OED_QUADRATIC_OBJECTIVE_HPP

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_SampledVector.hpp"
#include "ROL_OED_BilinearConstraint.hpp"
#include "ROL_OED_ProfiledClass.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class QuadraticObjective : public Objective_SimOpt<Real>, public ProfiledClass<Real,std::string> {
private:
  const Ptr<BilinearConstraint<Real>> M_;
  Ptr<Vector<Real>>                 g_;
  bool                   isInit_, isG_;

  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

  void initialize(const Vector<Real> &x);
  void apply(const Vector<Real> &u, const Vector<Real> &z, Real &tol);

public:
  QuadraticObjective(const Ptr<BilinearConstraint<Real>> &M);

  Ptr<BilinearConstraint<Real>> getM() const;

  void update( const Vector<Real> &u, const Vector<Real> &z, UpdateType type, int iter = -1 ) override;
  Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void hessVec_11( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void hessVec_12( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void hessVec_21( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void hessVec_22( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) override;
  void setParameter( const std::vector<Real> &param ) override;

}; // class QuadraticObjective

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_QuadraticObjective_Def.hpp"

#endif
