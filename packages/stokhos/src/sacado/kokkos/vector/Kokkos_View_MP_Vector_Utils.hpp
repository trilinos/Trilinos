// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_MP_VECTOR_UTILS_HPP
#define KOKKOS_VIEW_MP_VECTOR_UTILS_HPP

#include "Kokkos_View_Utils.hpp"
#include "Sacado_MP_Vector.hpp"

namespace Kokkos {

// Type name for a local, unmanaged view with possibly a different static size
template <typename ViewType,
          unsigned LocalSize,
          bool isStatic = Sacado::IsStaticallySized<typename ViewType::value_type>::value>
struct LocalMPVectorView {};

template <typename ViewType, unsigned LocalSize>
struct LocalMPVectorView<ViewType, LocalSize, false> {
  typedef ViewType type;
};

template <typename D, typename ... P, unsigned LocalSize>
struct LocalMPVectorView< View<D,P...>, LocalSize, true > {
  typedef typename Kokkos::Impl::ViewMapping< void, typename Kokkos::ViewTraits<D,P...>, Sacado::MP::VectorPartition<LocalSize> >::type type;
};

namespace Impl {

template< class OldStorageType , class Device >
struct RebindStokhosStorageDevice< Sacado::MP::Vector< OldStorageType > , Device >
{
  typedef typename
    OldStorageType::template apply<
      typename OldStorageType::ordinal_type ,
      typename OldStorageType::value_type ,
      Device >
    NewStorageApply ;

  typedef typename NewStorageApply::type NewStorageType ;
  typedef typename Sacado::MP::Vector< OldStorageType >::template apply< NewStorageType > NewVectorApply ;

  typedef typename NewVectorApply::type type ;
};

template< class OldStorageType , class Device >
struct RebindStokhosStorageDevice< const Sacado::MP::Vector< OldStorageType > , Device >
{
  typedef typename
    OldStorageType::template apply<
      typename OldStorageType::ordinal_type ,
      typename OldStorageType::value_type ,
      Device >
    NewStorageApply ;

  typedef typename NewStorageApply::type NewStorageType ;
  typedef typename Sacado::MP::Vector< OldStorageType >::template apply< NewStorageType > NewVectorApply ;

  typedef const typename NewVectorApply::type type ;
};

} // namespace Impl

// Whether a given type is a view with scalar type Sacado::MP::Vector
template <typename view_type>
struct is_view_mp_vector { static const bool value = false; };

} // namespace Kokkos

#endif // KOKKOS_VIEW_MP_VECTOR_UTILS_HPP
