// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(ROL_MiniTensor_BoundConstraint_hpp)
#define ROL_MiniTensor_BoundConstraint_hpp

#include "MiniTensor_Solvers.h"
#include "ROL_Bounds.hpp"
#include "ROL_MiniTensor_Vector.hpp"

namespace ROL {

using Index = minitensor::Index;

///
/// Function base class that defines the interface to Mini Solvers.
///
template<typename T, Index N>
class MiniTensor_BoundConstraint: public Bounds<T> {

public:

  MiniTensor_BoundConstraint() = delete;

  MiniTensor_BoundConstraint(
      MiniTensorVector<T, N> & lo,
      MiniTensorVector<T, N> & hi
  );
};

} // namespace ROL

#include "ROL_MiniTensor_BoundConstraint_Def.hpp"

#endif // ROL_MiniTensor_BoundConstraint_hpp
