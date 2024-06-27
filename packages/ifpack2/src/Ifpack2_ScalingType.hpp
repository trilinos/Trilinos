// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _IFPACK2_SCALINGTYPE_HPP_
#define _IFPACK2_SCALINGTYPE_HPP_
/*! \file Ifpack2_ScalingType.hpp
    \brief Ifpack2::ScalingType enumerable type
 */

namespace Ifpack2 {

//! Ifpack2 scaling type selector.
/*! Selects the type of scaling used (if any) for Ifpack2 preconditioners.
*/
enum ScalingType {
  None,
  LeftDiagonal,
  RightDiagonal, 
  SymmetricDiagonal,
  RowSum,
  ColSum, 
  RowAndColSum
};

}//namespace Ifpack2

#endif /* _IFPACK2_SCALINGTYPE_HPP_ */
