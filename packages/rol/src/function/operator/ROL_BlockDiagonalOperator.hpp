// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BLOCKDIAGONALOPERATOR_H
#define ROL_BLOCKDIAGONALOPERATOR_H

#include "ROL_LinearOperator.hpp"
#include "ROL_Ptr.hpp"
#include <vector>

/** @ingroup func_group
    \class ROL::BlockDiagonalOperator
    \brief Provides the interface to apply a block diagonal operator to a 
           partitioned vector.

    ---
*/

namespace ROL {

template<class Real>
class BlockDiagonalOperator : public LinearOperator<Real> {
private:
  const std::vector<Ptr<LinearOperator<Real>>> diag_;
  const unsigned numOp_;

public:
  BlockDiagonalOperator( const std::vector<Ptr<LinearOperator<Real>>> &diag );
  const Ptr<LinearOperator<Real>> get(unsigned i) const;
  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const override;
  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const override;
  void applyAdjoint( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const override;
  void applyAdjointInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const override;

}; // class BlockDiagonalOperator

} // namespace ROL

#include "ROL_BlockDiagonalOperator_Def.hpp"

#endif // ROL_BLOCKDIAGONALOPERATOR_H
