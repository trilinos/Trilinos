// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_ABS_HPP
#define RTOPPACK_TOP_ABS_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Transformation operator that takes absolute values of elements:
 * <tt>z0[i] = abs(v0[i]), i=0...n-1</tt>.
 */
RTOP_TOP_1_1( TOpAbs )
{
  z0 = ScalarTraits<Scalar>::magnitude(v0);
}


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ABS_HPP
