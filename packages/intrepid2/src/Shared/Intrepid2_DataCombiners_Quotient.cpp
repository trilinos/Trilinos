// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_DataCombiners.cpp
//
//  Created by Roberts, Nathan V on 6/1/23.
//

#include "Intrepid2_DataCombiners.hpp"

#include "Intrepid2_DataFunctors.hpp"

using DefaultDeviceType = Kokkos::DefaultExecutionSpace::device_type;

namespace Intrepid2
{
  template class DataCombiner<double,DefaultDeviceType,ScalarQuotientFunctor<double> >;
}
