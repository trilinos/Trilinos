#ifndef PHALANX_MDFIELD_EXTENT_TRAITS_HPP
#define PHALANX_MDFIELD_EXTENT_TRAITS_HPP

#include <type_traits>

namespace PHX {
  template<typename I> struct is_extent : std::false_type {};
}

#endif
