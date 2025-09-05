// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_OrientationOperatorDef.hpp
    \brief  Definition file for the Intrepid2::OrientationOperator class.
    \author Created by Nate Roberts
*/
#ifndef __INTREPID2_ORIENTATION_OPERATOR_DEF_HPP__
#define __INTREPID2_ORIENTATION_OPERATOR_DEF_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

template<class DeviceType>
OrientationOperator<DeviceType>::OrientationOperator(typename OrientationOperator<DeviceType>::UnmanagedOrdinalView rowIndices_,
                                                     typename OrientationOperator<DeviceType>::UnmanagedOrdinalView offsetsForRowOrdinal_,
                                                     typename OrientationOperator<DeviceType>::UnmanagedOrdinalView packedColumnIndices_,
                                                     typename OrientationOperator<DeviceType>::UnmanagedDoubleView  packedWeights_)
:
rowIndices(rowIndices_),
offsetsForRowOrdinal(offsetsForRowOrdinal_),
packedColumnIndices(packedColumnIndices_),
packedWeights(packedWeights_)
{}

template<class DeviceType>
OrientationOperator<DeviceType>::OrientationOperator()
{}

} // namespace Intrepid2

#endif
