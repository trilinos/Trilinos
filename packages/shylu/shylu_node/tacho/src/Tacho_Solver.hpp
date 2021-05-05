#ifndef __TACHO_SOLVER_HPP__
#define __TACHO_SOLVER_HPP__

/// \file Tacho_Solver.hpp
/// \brief solver interface
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho.hpp"
#include "Tacho_Driver.hpp"

namespace Tacho {

  template<typename ValueType, typename DeviceType> 
  using Solver = Driver<ValueType,typename UseThisDevice<typename DeviceType::execution_space>::type >;

}

#endif
