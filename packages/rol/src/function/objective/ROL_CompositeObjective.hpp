// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_COMPOSITE_OBJECTIVE_H
#define ROL_COMPOSITE_OBJECTIVE_H

#include "ROL_StdObjective.hpp"

/** @ingroup func_group
    \class ROL::CompositeObjective
    \brief Provides the interface to evaluate composite objective functions.
*/


namespace ROL {

template<typename Real>
class CompositeObjective : public Objective<Real> {
private:
  const std::vector<Ptr<Objective<Real>>> obj_vec_;
  const Ptr<StdObjective<Real>> std_obj_;

  Ptr<std::vector<Real>> obj_value_, obj_grad_, obj_gv_, obj_hess_;
  Ptr<StdVector<Real>> obj_value_vec_, obj_grad_vec_, obj_gv_vec_, obj_hess_vec_;
  std::vector<Ptr<Vector<Real>>> vec_grad_, vec_hess_;

  bool isInitialized_, isValueComputed_, isGradientComputed_;

public:
  CompositeObjective(const std::vector<Ptr<Objective<Real>>> &obj_vec,
                     const Ptr<StdObjective<Real>> &std_obj);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

// Definitions for parametrized (stochastic) objective functions
public:
  void setParameter(const std::vector<Real> &param) override;

private:
  void initialize(const Vector<Real> &x);
  void computeValue(const Vector<Real> &x, Real &tol);
  void computeGradient(const Vector<Real> &x, Real &tol);
  void computeHessVec(const Vector<Real> &v, const Vector<Real> &x, Real &tol);
};

} // namespace ROL

#include "ROL_CompositeObjective_Def.hpp"

#endif
