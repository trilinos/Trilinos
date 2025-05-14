// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_Data.cpp
//
//  Created by Roberts, Nathan V on 5/30/23.
//

#include "Intrepid2_Data.hpp"

using DefaultDeviceType = Kokkos::DefaultExecutionSpace::device_type;

// a HIP ldd bug results in it not being able to find Data::isDiagonal() -- it claims it is hidden and undefined even though readelf -s --wide libintrepid2.a | grep isDiagonal() | c++filt shows that it is WEAK DEFAULT (i.e., visible, and defined) -- when we do ETI on Intrepid2::Data.  So for now we disable this ETI on HIP.
#if !defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HIP)
template class Intrepid2::Data<double,DefaultDeviceType>;
#endif
