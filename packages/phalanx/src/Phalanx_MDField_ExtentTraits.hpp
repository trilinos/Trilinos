#ifndef PHALANX_MDFIELD_EXTENT_TRAITS_HPP
#define PHALANX_MDFIELD_EXTENT_TRAITS_HPP

#include "Teuchos_TypeNameTraits.hpp"
#include <type_traits>
#include <string>

namespace PHX {

  /// Identifies that a user defined struct is a dimension template parameter for MDFields. Users must specialize for true_types.
  template<typename Extent> struct is_extent : std::false_type {};

  /// Print statement for a user defined struct dimension template parameter for MDFields. Users can specialize this.
  template<typename Extent> std::string print()
  {return Teuchos::demangleName(typeid(Extent).name());}

}



#endif
