// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_NEW_CONSTRAINT_MANAGER_H
#define ROL_NEW_CONSTRAINT_MANAGER_H

#include "ROL_Constraint_Partitioned.hpp"
#include "ROL_BoundConstraint_Partitioned.hpp"
#include <unordered_map>

/** @ingroup func_group
    \class ROL::NewConstraintManager
    \brief Provides a wrapper for multiple constraints.

    ---
*/

namespace ROL {

template<typename Real>
struct ConstraintData {
  const Ptr<Constraint<Real>>      constraint;
  const Ptr<Vector<Real>>          multiplier;
  const Ptr<Vector<Real>>          residual;
  const Ptr<BoundConstraint<Real>> bounds;

  ConstraintData(const Ptr<Constraint<Real>>      &con,
                 const Ptr<Vector<Real>>          &mul,
                 const Ptr<Vector<Real>>          &res=nullPtr,
                 const Ptr<BoundConstraint<Real>> &bnd=nullPtr)
    : constraint(con), multiplier(mul),
      residual(res==nullPtr ? mul->dual().clone() : res), bounds(bnd) {}
};

template<typename Real>
class NewConstraintManager {
private:
  Ptr<Constraint<Real>>      con_;
  Ptr<Vector<Real>>          mul_;
  Ptr<Vector<Real>>          res_;
  Ptr<Constraint<Real>>      linear_con_;
  Ptr<Vector<Real>>          linear_mul_;
  Ptr<Vector<Real>>          linear_res_;
  Ptr<Vector<Real>>          xprim_;
  Ptr<Vector<Real>>          xdual_;
  Ptr<BoundConstraint<Real>> bnd_;

  std::vector<Ptr<Constraint<Real>>>      cvec_;  // General constraints
  std::vector<Ptr<Vector<Real>>>          lvec_;  // General multiplier vector
  std::vector<Ptr<Vector<Real>>>          rvec_;  // General residual vector
  std::vector<Ptr<Constraint<Real>>>      lcvec_; // Linear constraints
  std::vector<Ptr<Vector<Real>>>          llvec_; // Linear multiplier vector
  std::vector<Ptr<Vector<Real>>>          lrvec_; // Linear residual vector
  std::vector<Ptr<Vector<Real>>>          psvec_; // Primal slack vector
  std::vector<Ptr<Vector<Real>>>          dsvec_; // Dual slack vector
  std::vector<Ptr<BoundConstraint<Real>>> sbnd_;  // Slack bound constraint

  std::vector<bool> isInequality_, isLinearInequality_;

  bool isNull_;
  bool hasInequality_;

  void initializeSlackVariable(const Ptr<Constraint<Real>>      &con,
                               const Ptr<BoundConstraint<Real>> &cbnd,
                               const Ptr<Vector<Real>>          &s,
                               const Ptr<Vector<Real>>          &x) const;

  void initialize(const std::unordered_map<std::string,ConstraintData<Real>> &input_con,
                  const Ptr<Vector<Real>>                                    &xprim,
                  const Ptr<Vector<Real>>                                    &xdual,
                  const Ptr<BoundConstraint<Real>>                           &bnd);

  void initialize(const std::unordered_map<std::string,ConstraintData<Real>> &input_con,
                  const std::unordered_map<std::string,ConstraintData<Real>> &input_lcon,
                  const Ptr<Vector<Real>>                                    &xprim,
                  const Ptr<Vector<Real>>                                    &xdual,
                  const Ptr<BoundConstraint<Real>>                           &bnd);

public:
  virtual ~NewConstraintManager(void) {}

  NewConstraintManager(const std::unordered_map<std::string,ConstraintData<Real>> &con,
                       const Ptr<Vector<Real>>                                    &xprim,
                       const Ptr<Vector<Real>>                                    &xdual,
                       const Ptr<BoundConstraint<Real>>                           &bnd = nullPtr);

  NewConstraintManager(const std::unordered_map<std::string,ConstraintData<Real>> &con,
                       const std::unordered_map<std::string,ConstraintData<Real>> &linear_con,
                       const Ptr<Vector<Real>>                                    &xprim,
                       const Ptr<Vector<Real>>                                    &xdual,
                       const Ptr<BoundConstraint<Real>>                           &bnd = nullPtr);

  const Ptr<Constraint<Real>>      getConstraint(void) const;
  const Ptr<Vector<Real>>          getMultiplier(void) const;
  const Ptr<Vector<Real>>          getResidual(void) const;
  const Ptr<Constraint<Real>>      getLinearConstraint(void) const;
  const Ptr<Vector<Real>>          getLinearMultiplier(void) const;
  const Ptr<Vector<Real>>          getLinearResidual(void) const;
  const Ptr<Vector<Real>>          getOptVector(void) const;
  const Ptr<Vector<Real>>          getDualOptVector(void) const;
  const Ptr<BoundConstraint<Real>> getBoundConstraint(void) const;

  bool isNull(void) const;
  bool hasInequality(void) const;

  void resetSlackVariables(void);

}; // class NewConstraintManager

} // namespace ROL

#include "ROL_NewConstraintManager_Def.hpp"

#endif
