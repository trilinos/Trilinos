#ifndef STK_SEARCH_UTIL_INTREPID2_BASIS_FACTORY_H
#define STK_SEARCH_UTIL_INTREPID2_BASIS_FACTORY_H

#include "stk_topology/topology_decl.hpp"
#include <memory>
#include "Intrepid2_Basis.hpp"


namespace stk::search {

using BasisExecSpace = Kokkos::Serial;
using I2DeviceType = typename BasisExecSpace::device_type;
using BasisPtr = std::shared_ptr<Intrepid2::Basis<I2DeviceType,double,double>>;
using BasisViewType = Kokkos::DynRankView<double, typename BasisExecSpace::memory_space>;

BasisPtr lagrangeBasisFactory(stk::topology stkTopo, bool useCompositeTet10 = false);

}

#endif
