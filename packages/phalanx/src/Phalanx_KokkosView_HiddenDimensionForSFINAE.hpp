#ifndef PHALANX_KOKKOS_VIEW_HIDDEN_DIMENSION_FOR_SFINAE_HPP
#define PHALANX_KOKKOS_VIEW_HIDDEN_DIMENSION_FOR_SFINAE_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_any.hpp"
#include "Sacado.hpp"
#include <vector>

namespace PHX {
  /// Returns true if the Scalar type has a hidden dimension (e.g. an
  /// AD array for derivatives) where the size must be supplied during
  /// runtime construction of the object. This is typically used for
  /// AD types where the derivative dimension is only known at
  /// runtime. We prune out SFad since its derivative dimension is
  /// known at compiletime.
  template<typename ScalarT> struct requires_dynamic_hidden_dimension;

  // Default implementation for query of dynamic hidden dimension
  template<typename ScalarT> struct requires_dynamic_hidden_dimension : std::false_type {};

  // Specializations for DFad and SLFad. Note that SFad size is chosen
  // at compiletime and can use the default case above for
  // construction.
  template<typename ScalarT> struct requires_dynamic_hidden_dimension<Sacado::Fad::DFad<ScalarT>> : std::true_type {};
  template<typename ScalarT> struct requires_dynamic_hidden_dimension<Sacado::ELRCacheFad::DFad<ScalarT>> : std::true_type {};
  template<typename ScalarT,int N> struct requires_dynamic_hidden_dimension<Sacado::Fad::SLFad<ScalarT,N>> : std::true_type {};
}

#endif
