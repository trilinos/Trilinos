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

//template class Intrepid2::Data<double,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::Data<double,DefaultDeviceType>;
