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
    using UnmanagedOrdinalView = Kokkos::View<ordinal_type*, typename DeviceType::memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using UnmanagedDoubleView  = Kokkos::View<      double*, typename DeviceType::memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    
    // only stores deviations from the identity
    UnmanagedOrdinalView rowIndices;           // index in basis (the field ordinal)
    UnmanagedOrdinalView offsetsForRowOrdinal; // argument is same as _rowIndices; offset gives index for _packedColumnIndices and _packedWeights
    UnmanagedOrdinalView packedColumnIndices;  // ordinal is the index from offsets
    UnmanagedDoubleView packedWeights;         // ordinal is the index from offsets
    
    // offsetsForRowOrdinal may be empty; if it is, isPermutation is set to true; otherwise it is false
    // if isWeightedPermutation is true, the offsets are all unit-spaced (0,1,…), and packedColumnIndices and packedWeights have the same length as rowIndices
    bool isWeightedPermutation;
  public:
    //! general constructor.  Sets isWeightedPermutation to false.
    KOKKOS_INLINE_FUNCTION
    OrientationOperator(UnmanagedOrdinalView rowIndices_,
                        UnmanagedOrdinalView offsetsForRowOrdinal_,
                        UnmanagedOrdinalView packedColumnIndices_,
                        UnmanagedDoubleView       packedWeights_);
    
    //! weighted-permutation constructor.  Sets isWeightedPermutation to true.
    KOKKOS_INLINE_FUNCTION
    OrientationOperator(UnmanagedOrdinalView rowIndices_,
                        UnmanagedOrdinalView packedColumnIndices_,
                        UnmanagedDoubleView       packedWeights_);
    
    KOKKOS_INLINE_FUNCTION
    OrientationOperator();
  };
}

// include templated function definitions
#include "Intrepid2_OrientationOperatorDef.hpp"

#endif
