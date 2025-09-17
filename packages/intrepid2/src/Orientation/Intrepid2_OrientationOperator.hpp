// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_OrientationOperator.hpp
    \brief  Header file for the Intrepid2::OrientationOperator class.
    \author Created by Nate Roberts
*/

#ifndef __INTREPID2_ORIENTATION_OPERATOR_HPP__
#define __INTREPID2_ORIENTATION_OPERATOR_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"

namespace Intrepid2 {

  /**
    \brief OrientationOperator
   
   Introduced in 2025, OrientationOperator seeks to make the application of orientations to vectors of local basis coefficients more efficient.
   Each orientation operator corresponds to a subcell (face or edge) of a basis, and applies the specified orientation map to the part of the basis
   corresponding to the subcell.  The assumption is that the input and output vectors supplied prior to application of any orientation operators
   are identical, allowing the identity map to be implemented as a no-op.
  */
  template<class DeviceType>
  class OrientationOperator {
  public:
    // only stores deviations from the identity
    Kokkos::View<ordinal_type*,DeviceType> rowIndices; // index in basis
    Kokkos::View<ordinal_type*,DeviceType> offsetsForRowOrdinal; // ordinal into _rowIndices; offset gives index for _packedColumnIndices and _packedWeights
    Kokkos::View<ordinal_type*,DeviceType> packedColumnIndices; // ordinal into _rowIndices
    Kokkos::View<double*,DeviceType> packedWeights; // ordinal into _rowIndices
  public:
    OrientationOperator(Kokkos::View<ordinal_type*,DeviceType> rowIndices_,
                        Kokkos::View<ordinal_type*,DeviceType> offsetsForRowOrdinal_,
                        Kokkos::View<ordinal_type*,DeviceType> packedColumnIndices_,
                        Kokkos::View<double*,DeviceType>       packedWeights_);
    
    OrientationOperator() {}
  };
}

// include templated function definitions
#include "Intrepid2_OrientationOperatorDef.hpp"

#endif
