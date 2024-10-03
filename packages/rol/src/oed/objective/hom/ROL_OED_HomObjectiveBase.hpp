// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOM_OBJECTIVE_BASE_HPP
#define ROL_OED_HOM_OBJECTIVE_BASE_HPP

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_VectorController.hpp"

#include "ROL_OED_BilinearConstraint.hpp"
#include "ROL_OED_LinearObjective.hpp"

namespace ROL {
namespace OED {
namespace Hom {

//
// DESCRIPTION:
//   Objective:   c(u) where c is a linear form
//   Constraint:  a(u,z,v) = c(v) for all v
//                where a(u,z,l) = a(l,z,u) is a trilinear form
//
//   L(u,z,l) = c(u) + a(u,z,l) - c(l) = c(u-l) + a(u,z,l)
//   dL/du    = c(.) + a(l,z,.)
//   dL/dz    = a(u,.,l)
//   d2L/dudz = a(l,.,.)
//   d2L/du2  = 0
//   d2L/dz2  = 0
//
//   J(z)    = c(S(z)) where S(z) = u solves a(u,z,v) = c(v) for all v
//   J'(z)   = -a(S(z),.,S(z))
//   J''(z)h = -2 a(S(z),.,W(z,h)) where W(z,h) = w solves
//              a(w,z,v) = -a(S(z),h,v) for all v
//

template<typename Real, typename Key>
class ObjectiveBase : public Objective<Real> {
private:
  // Storage for state variables
  const Ptr<VectorController<Real,Key>> stateStore_;

  // OED Objective and Constraints
  Ptr<BilinearConstraint<Real>> con_;
  Ptr<LinearObjective<Real>>  obj_;

  // Vector storage
  Ptr<Vector<Real>> state_, state_sens_, res_;

  // Update information
  UpdateType updateType_;
  int  updateIter_;

  // Storage for state variables
  bool storage_;             
  bool doUpdate_;

protected:
  void setConstraint(const Ptr<BilinearConstraint<Real>> &con);
  void setObjective(const Ptr<LinearObjective<Real>> &obj);
  void setStorage(bool storage);
  void setUpdate(bool doUpdate);
  void initialize(const Ptr<Vector<Real>> &state);
  const Ptr<BilinearConstraint<Real>> getConstraint() const;
  const Ptr<LinearObjective<Real>> getObjective() const;
  const Ptr<Vector<Real>> getState() const;
  const Ptr<Vector<Real>> getStateSens() const;
  void solve_state_equation(const Key &param,
                            const Vector<Real> &z,
                                  Real &tol);
  void solve_state_sensitivity(const Vector<Real> &v,
                               const Vector<Real> &z,
                                     Real &tol);

public:
  ObjectiveBase();

  virtual void update( const Vector<Real> &z,
                       UpdateType type,
                       int iter = -1 );

  virtual void setParameter( const std::vector<Real> &param );
};

} // End Hom Namespace
} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_HomObjectiveBase_Def.hpp"

#endif
