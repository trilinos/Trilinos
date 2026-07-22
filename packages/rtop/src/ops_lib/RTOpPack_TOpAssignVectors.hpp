// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_ASSIGN_VECTORS_HPP
#define RTOPPACK_TOP_ASSIGN_VECTORS_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {


/** \brief VectorBase assignment transformation operator: <tt>z0[i] = v0[i],
 * i=0...n-1</tt>.
 */
RTOP_TOP_1_1( TOpAssignVectors )
{
  z0 = v0;
}


} // namespace RTOpPack


#endif // RTOPPACK_TOP_ASSIGN_VECTORS_HPP
