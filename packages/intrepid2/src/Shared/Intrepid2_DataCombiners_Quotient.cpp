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
