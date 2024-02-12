//
//  Intrepid2_Data.cpp
//
//  Created by Roberts, Nathan V on 5/30/23.
//

#include "Intrepid2_Data.hpp"

using DefaultDeviceType = Kokkos::DefaultExecutionSpace::device_type;

//template class Intrepid2::Data<double,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::Data<double,DefaultDeviceType>;
