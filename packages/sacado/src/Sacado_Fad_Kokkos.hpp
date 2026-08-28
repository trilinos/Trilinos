// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_KOKKOS_HPP
#define SACADO_FAD_KOKKOS_HPP

#include "Sacado_ConfigDefs.h"

#if defined(HAVE_SACADO_KOKKOS)

#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Kokkos_Core_fwd.hpp"
#include "Sacado_Fad_Kokkos_View_Support.hpp"

#else

#include "Sacado_Fad_Kokkos_LayoutContiguous.hpp"

// Some default traits definitions that need to be defined even if the view
// specialization is disabled
namespace Sacado {

template <typename ViewType, typename Enabled = void>
struct ThreadLocalScalarType {
  typedef typename ViewType::non_const_value_type type;
};

template <typename ViewType>
struct ViewScalarStride {
  static const unsigned stride =
    Sacado::Impl::LayoutScalarStride< typename ViewType::array_layout>::stride;
  static const bool is_unit_stride =
    Sacado::Impl::LayoutScalarStride< typename ViewType::array_layout>::is_unit_stride;
};

template <typename ViewType> struct is_view_fad_contiguous {
  static const bool value = false;
};
template <typename ViewType> struct is_dynrankview_fad_contiguous {
  static const bool value = false;
};

} // namespace Sacado

// Inject some things into Kokkos for backwards compatibility
#ifdef SACADO_ENABLE_DEPRECATED_CODE
namespace Kokkos {

template<typename ViewType, typename Enabled = void>
using ThreadLocalScalarType SACADO_DEPRECATED_WITH_COMMENT(
    "Use Sacado::ThreadLocalScalarType instead of Kokkos::ThreadLocalScalarType") = Sacado::ThreadLocalScalarType<ViewType, Enabled>;
template<typename ViewType>
using ViewScalarStride SACADO_DEPRECATED_WITH_COMMENT(
    "Use Sacado::ViewScalarStride instead of Kokkos::ViewScalarStride") = Sacado::ViewScalarStride<ViewType>;
template<typename ViewType>
using is_view_fad_contiguous SACADO_DEPRECATED_WITH_COMMENT(
    "Use Sacado::is_view_fad_contiguous instead of Kokkos::is_view_fad_contiguous") = Sacado::is_view_fad_contiguous<ViewType>;
template<typename ViewType>
using is_dynrankview_fad_contiguous SACADO_DEPRECATED_WITH_COMMENT(
    "Use Sacado::is_dynrankview_fad_contiguous instead of Kokkos::is_dynrankview_fad_contiguous") = Sacado::is_dynrankview_fad_contiguous<ViewType>;

} // namespace Kokkos
#endif

namespace Sacado {

  namespace Fad {

    /* Define a partition of a View of Sacado Fad type */
    template <unsigned Size = 0>
    struct Partition {
      static const unsigned PartitionSize = Size;
      unsigned offset ;
      unsigned stride ;

      template< typename iType0 , typename iType1 >
      KOKKOS_INLINE_FUNCTION
      Partition( const iType0 & i0 , const iType1 & i1 ) :
        offset(i0), stride(i1) {
      }
    };

    template <typename T>
    struct is_fad_partition {
      static const bool value = false;
    };

    template <unsigned Stride>
    struct is_fad_partition< Partition<Stride> > {
      static const bool value = true;
    };

  }

  // Type of local scalar type when partitioning a view
  template <typename T, unsigned Stride = 0>
  struct LocalScalarType {
    typedef T type;
  };
  template <typename T, unsigned Stride>
  struct LocalScalarType<const T, Stride> {
    typedef typename LocalScalarType<T,Stride>::type lst;
    typedef const lst type;
  };

  // For DFad, the size is not part of the type, so the default implementation
  // is sufficient

  // Type of local scalar type when partitioning a view
  //
  // For SLFad, divde the array size by the given stride
  namespace Fad {
  namespace Exp {
    template <typename T, int N> class StaticStorage;
    template <typename S> class GeneralFad;
  }
  }
  template <typename T, int N, unsigned Stride>
  struct LocalScalarType< Fad::Exp::GeneralFad< Fad::Exp::StaticStorage<T,N> >,
                          Stride > {
    static const int Ns = (N+Stride-1) / Stride;
    typedef Fad::Exp::GeneralFad< Fad::Exp::StaticStorage<T,Ns> > type;
  };
#ifndef SACADO_NEW_FAD_DESIGN_IS_DEFAULT
  namespace Fad {
    template <typename T, int N> class SLFad;
  }
  template <typename T, int N, unsigned Stride>
  struct LocalScalarType< Fad::SLFad<T,N>, Stride > {
    static const int Ns = (N+Stride-1) / Stride;
    typedef Fad::SLFad<T,Ns> type;
  };
#endif

  // Type of local scalar type when partitioning a view
  //
  // For SFad, divde the array size by the given stride.  If it divides evenly,
  // use SFad, otherwise use SLFad
  namespace Fad {
  namespace Exp {
    template <typename T, int N> class StaticFixedStorage;
    template <typename T, int N> class StaticStorage;
    template <typename S> class GeneralFad;
  }
  }
  template <typename T, int N, unsigned Stride>
  struct LocalScalarType< Fad::Exp::GeneralFad< Fad::Exp::StaticFixedStorage<T,N> >,
                          Stride > {
    static const int Ns = (N+Stride-1) / Stride;
    typedef typename std::conditional<
      Ns == N/Stride ,
      Fad::Exp::GeneralFad< Fad::Exp::StaticFixedStorage<T,Ns> > ,
      Fad::Exp::GeneralFad< Fad::Exp::StaticStorage<T,Ns> >
    >::type type;
  };

#ifndef SACADO_NEW_FAD_DESIGN_IS_DEFAULT
  namespace Fad {
    template <typename T, int N> class SFad;
  }
  template <typename T, int N, unsigned Stride>
  struct LocalScalarType< Fad::SFad<T,N>, Stride > {
    static const int Ns = (N+Stride-1) / Stride;
    typedef typename std::conditional< Ns == N/Stride , Fad::SFad<T,Ns> , Fad::SLFad<T,Ns> >::type type;
  };
#endif

  template <unsigned Stride, typename T>
  KOKKOS_INLINE_FUNCTION
  const T& partition_scalar(const T& x) { return x; }

} // namespace Sacados

#endif // defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#endif // defined(HAVE_SACADO_KOKKOS)

#endif /* #ifndef SACADO_FAD_KOKKOS_HPP */
