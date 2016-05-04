#ifndef PANZER_VIEW_FACTORY_HPP
#define PANZER_VIEW_FACTORY_HPP

#include "KokkosExp_ViewFactory.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer {

  //! Wrapper to simplify Panzer use of Sacado ViewFactory 
  template<typename InputArray,typename ... DimensionPack>
  Kokkos::DynRankView<typename InputArray::value_type,PHX::Device>
  createDynRankView(const InputArray& a, const std::string& name, const DimensionPack... dims)
  {
    using return_type = Kokkos::DynRankView<typename InputArray::value_type,PHX::Device>;
    return Kokkos::ViewFactory<InputArray>::template create_view<return_type>(a,name,dims...);
  }

}

#endif
