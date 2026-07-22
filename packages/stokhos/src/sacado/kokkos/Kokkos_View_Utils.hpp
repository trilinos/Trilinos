// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_VIEW_UTILS_HPP
#define KOKKOS_VIEW_UTILS_HPP

#include <stdexcept>

// We are hooking into Kokkos Core internals here
// Need to define this macro since we include non-public headers
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif
#include "Kokkos_View.hpp"
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif

namespace Kokkos {

namespace Impl {

KOKKOS_INLINE_FUNCTION
void raise_error(const char *msg)
{
  KOKKOS_IF_ON_HOST(throw std::runtime_error(msg);)

  KOKKOS_IF_ON_DEVICE(Kokkos::abort(msg);)
}

template< class T , class Device > struct RebindStokhosStorageDevice ;

template< class T , class Device >
struct RebindStokhosStorageDevice< T * , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type * type ;
};

template< class T , class Device >
struct RebindStokhosStorageDevice< T [] , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type * type ;
};

template< class T , unsigned N , class Device >
struct RebindStokhosStorageDevice< T[N] , Device >
{
  typedef typename RebindStokhosStorageDevice< T , Device >::type type[N] ;
};

// Get Sacado size from a list of dimensions
template <unsigned Rank> struct GetSacadoSize {};
template <> struct GetSacadoSize<0> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 = 0 ,
                      const size_t n2 = 0 ,
                      const size_t n3 = 0 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n0;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[0];
  }
};
template <> struct GetSacadoSize<1> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 = 0 ,
                      const size_t n3 = 0 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n1;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[1];
  }
};
template <> struct GetSacadoSize<2> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 = 0 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n2;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[2];
  }
};
template <> struct GetSacadoSize<3> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 = 0 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n3;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[3];
  }
};
template <> struct GetSacadoSize<4> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 = 0 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n4;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[4];
  }
};
template <> struct GetSacadoSize<5> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 ,
                      const size_t n6 = 0 ,
                      const size_t n7 = 0 ) {
    return n5;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[5];
  }
};
template <> struct GetSacadoSize<6> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 ,
                      const size_t n6 ,
                      const size_t n7 = 0 ) {
    return n6;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[6];
  }
};
template <> struct GetSacadoSize<7> {
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const size_t n0 ,
                      const size_t n1 ,
                      const size_t n2 ,
                      const size_t n3 ,
                      const size_t n4 ,
                      const size_t n5 ,
                      const size_t n6 ,
                      const size_t n7 ) {
    return n7;
  }

  template <typename Layout>
  KOKKOS_INLINE_FUNCTION
  static size_t eval( const Layout& layout ) {
    return layout.dimension[7];
  }
};

} // namespace Impl

// Typename of flat array where sacado dimension is folded into neighbor
template <typename view_type, typename Enabled = void>
struct FlatArrayType {
  typedef view_type type;
};

// Typename of the intrinsic scalar type in a view
template <typename view_type, typename Enabled = void>
struct IntrinsicScalarType {
  typedef typename view_type::array_type::non_const_value_type type;
};

template <typename ViewType>
ViewType
make_view(const std::string& label,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(label, N0, N1, N2, N3, N4, N5, N6, N7);
}

template <typename ViewType>
ViewType
make_view(const ViewAllocateWithoutInitializing& init,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(init, N0, N1, N2, N3, N4, N5, N6, N7);
}

template <typename ViewType>
ViewType
make_view(typename ViewType::pointer_type ptr,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(ptr, N0, N1, N2, N3, N4, N5, N6, N7);
}

template <typename ViewType>
ViewType
make_view(const std::string& label,
          const Impl::WithoutInitializing_t& init,
          size_t N0 = 0, size_t N1 = 0, size_t N2 = 0, size_t N3 = 0,
          size_t N4 = 0, size_t N5 = 0, size_t N6 = 0, size_t N7 = 0)
{
  return ViewType(view_alloc(label,init),
                  N0, N1, N2, N3, N4, N5, N6, N7);
}

namespace Impl {

// Specialization for deep_copy( view, view::value_type ) for Cuda

template <class OutputView, typename Enabled = void>
struct StokhosViewFill;

} // namespace Impl

} // namespace Kokkos

#endif // KOKKOS_VIEW_UTILS_HPP
