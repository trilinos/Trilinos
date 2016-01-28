#ifndef PHALANX_MDFIELD_TYPE_TRAITS_HPP
#define PHALANX_MDFIELD_TYPE_TRAITS_HPP

#include "Kokkos_View_Fad.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Sacado.hpp"

namespace PHX {

  template<typename ViewType>
  struct MDFieldTypeTraits {
    typedef typename ViewType::value_type ScalarT;
    typedef typename ViewType::reference_type return_type;
  };

}

#endif
