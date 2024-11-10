// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_INTEGERTRANSFORMATION_H
#define ROL_PEBBL_INTEGERTRANSFORMATION_H

#include "ROL_Constraint.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_PEBBL_MixedVector.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::IntegerTransformation
    \brief Defines the pebbl integer transformation operator interface.

    ROL's pebbl constraint interface is designed to set individual components
    of a vector to a fixed value.  The range space is the same as the domain.

    ---
*/


namespace ROL {
namespace PEBBL {

template <class Real>
class IntegerTransformation : public Constraint<Real> {
private:
  Ptr<Vector<Real>> getOptVector( Vector<Real> &xs ) const {
    try {
      return dynamic_cast<PartitionedVector<Real>&>(xs).get(0);
    }
    catch (std::exception &e) {
      return makePtrFromRef(xs);
    }
  }

  Ptr<Vector<Real>> getIntegerVector(Vector<Real> &xs) const {
    try {
      return dynamic_cast<MixedVector<Real>&>(*getOptVector(xs)).getIntegerVariables();
    }
    catch (std::exception &e) {
      return getOptVector(xs);
    }
  }

protected:
  std::map<int,Real> map_;

public:
  IntegerTransformation(void) {}
  IntegerTransformation(const IntegerTransformation &T) : map_(T.map_) {}

  virtual void fixValues(Vector<Real> &x, bool zero = false) const = 0;

  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
    c.set(x);
    Ptr<Vector<Real>> cp = getIntegerVector(c);
    fixValues(*cp,false);
  }

  void applyJacobian(Vector<Real> &jv,
               const Vector<Real> &v,
               const Vector<Real> &x,
                     Real &tol) {
    jv.set(v);
    Ptr<Vector<Real>> jvp = getIntegerVector(jv);
    fixValues(*jvp,true);
  }

  void applyAdjointJacobian(Vector<Real> &ajv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                            Real &tol) {
    ajv.set(v);
    Ptr<Vector<Real>> ajvp = getIntegerVector(ajv);
    fixValues(*ajvp,true);
  }

  void applyAdjointHessian(Vector<Real> &ahuv,
                     const Vector<Real> &u,
                     const Vector<Real> &v,
                     const Vector<Real> &x,
                           Real &tol) {
    ahuv.zero();
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

}; // class IntegerTransformation

} // namespace PEBBL
} // namespace ROL

#endif
