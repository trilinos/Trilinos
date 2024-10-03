// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_INTEGERCONSTRAINT_H
#define ROL_PEBBL_INTEGERCONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_StdVector.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::IntegerConstraint
    \brief Defines the pebbl integer constraint operator interface.

    ROL's pebbl constraint interface is designed to set individual components
    of a vector to a fixed value.  The range space is a ROL::StdVector.

    ---
*/


namespace ROL {
namespace PEBBL {

template <class Real>
class IntegerConstraint : public Constraint<Real> {
private:
  std::map<int,Real> map_;

  Ptr<std::vector<Real>> getData(Vector<Real> &x) const {
    return dynamic_cast<StdVector<Real>&>(x).getVector();
  }

  Ptr<const std::vector<Real>> getConstData(const Vector<Real> &x) const {
    return dynamic_cast<const StdVector<Real>&>(x).getVector();
  }

public:
  IntegerConstraint(void) {}
  IntegerConstraint(const IntegerConstraint &con) : map_(con.map_) {}

  void value(Vector<Real> &c,
       const Vector<Real> &x,
             Real &tol) {
    int cnt = 0;
    for (auto it=map_.begin(); it!=map_.end(); ++it) {
      (*getData(c))[cnt] = x.dot(*x.basis(it->first))-it->second;
      cnt++;
    }
  }

  void applyJacobian(Vector<Real> &jv,
               const Vector<Real> &v,
               const Vector<Real> &x,
                     Real &tol) {
    int cnt = 0;
    for (auto it=map_.begin(); it!=map_.end(); ++it) {
      (*getData(jv))[cnt] = v.dot(*x.basis(it->first));
      cnt++;
    }
  }

  void applyAdjointJacobian(Vector<Real> &ajv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                            Real &tol) {
    int cnt = 0;
    for (auto it=map_.begin(); it!=map_.end(); ++it) {
      ajv.axpy((*getConstData(v))[cnt],*x.basis(it->first));
      cnt++;
    }
  }

  void applyAdjointHessian(Vector<Real> &ahuv,
                     const Vector<Real> &u,
                     const Vector<Real> &v,
                     const Vector<Real> &x,
                           Real &tol) {
    ahuv.zero();
  }

  Ptr<Vector<Real> > makeConstraintVector(void) const {
    return makePtr<StdVector<Real>>(makePtr<std::vector<Real>>(map_.size(),0));
  }

  bool isEmpty(void) const {
    return map_.empty();
  }

  void reset(void) {
    map_.clear();
  }

  void add(const std::map<int,Real> &in) {
    map_.insert(in.begin(),in.end());
  }

  void add(const std::pair<int,Real> &in) {
    map_.insert(in);
  }

}; // class IntegerConstraint

} // namespace PEBBL
} // namespace ROL

#endif
