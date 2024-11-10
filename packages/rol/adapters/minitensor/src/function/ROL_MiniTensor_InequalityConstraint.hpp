// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(ROL_MiniTensor_InequalityConstraint_hpp)
#define ROL_MiniTensor_InequalityConstraint_hpp

#include "MiniTensor_Solvers.h"
#include "ROL_Constraint.hpp"
#include "ROL_MiniTensor_Vector.hpp"

namespace ROL {

using Index = minitensor::Index;

///
/// Function base class that defines the interface to Mini Solvers.
///
template<typename MSIC, typename S, Index M, Index N>
class MiniTensor_InequalityConstraint : public Constraint<S>
{
public:

  MiniTensor_InequalityConstraint(MSIC & msec);

  MiniTensor_InequalityConstraint() = delete;

  virtual
  ~MiniTensor_InequalityConstraint();

  // ROL interface
  virtual
  void
  value(Vector<S> & c, Vector<S> const & x, S & tol) final;

  virtual
  void
  applyJacobian(Vector<S> & jv, Vector<S> const & v,
      Vector<S> const & x, S & tol) final;

  virtual
  void
  applyAdjointJacobian(Vector<S> & ajv, Vector<S> const & v,
      Vector<S> const & x, S & tol) final;

private:
  MSIC
  minisolver_ic_;
};
} // namespace ROL

#include "ROL_MiniTensor_InequalityConstraint_Def.hpp"

#endif // ROL_MiniTensor_InequalityConstraint_hpp
