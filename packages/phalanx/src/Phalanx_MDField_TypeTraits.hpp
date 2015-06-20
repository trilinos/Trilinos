#ifndef PHALANX_MDFIELD_TYPE_TRAITS_HPP
#define PHALANX_MDFIELD_TYPE_TRAITS_HPP

#include "Kokkos_View_Fad.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Sacado.hpp"

namespace PHX {

  template<typename ViewType>
  struct MDFieldTypeTraits {
    typedef typename ViewType::value_type ScalarT;
    typedef ScalarT& return_type;
  };

  template<typename S, typename L, typename D, typename M>
  struct MDFieldTypeTraits< Kokkos::View<S,L,D,M,Kokkos::Impl::ViewSpecializeSacadoFad> > {
    typedef Kokkos::View<S,L,D,M,Kokkos::Impl::ViewSpecializeSacadoFad> ViewType;
    // dimension and stride
    typedef typename ViewType::reference_type return_type;
  };

}

#endif
