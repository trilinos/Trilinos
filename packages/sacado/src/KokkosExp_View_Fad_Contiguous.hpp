// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef KOKKOS_EXPERIMENTAL_VIEW_SACADO_FAD_CONTIGUOUS_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_SACADO_FAD_CONTIGUOUS_HPP

#include "Sacado_ConfigDefs.h"
#if defined(HAVE_SACADO_KOKKOS)

#include "Kokkos_LayoutContiguous.hpp"

// Some default traits definitions that need to be defined even if the view
// specialization is disabled
namespace Kokkos {

template <typename ViewType, typename Enabled = void>
struct ThreadLocalScalarType {
  typedef typename ViewType::non_const_value_type type;
};

template <typename ViewType>
struct ViewScalarStride {
  static const unsigned stride =
    Impl::LayoutScalarStride< typename ViewType::array_layout>::stride;
  static const bool is_unit_stride =
    Impl::LayoutScalarStride< typename ViewType::array_layout>::is_unit_stride;
};

} // namespace Kokkos

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

} // namespace Sacado

// Make sure the user really wants these View specializations
#if defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Sacado_Traits.hpp"
#include "Kokkos_Core.hpp"
#if KOKKOS_VERSION >= 40499
#include "View/Kokkos_ViewMapping.hpp"
#else
#include "impl/Kokkos_ViewMapping.hpp"
#endif

//----------------------------------------------------------------------------

namespace Sacado {

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  namespace Fad {
  namespace Exp {
    template <typename T, typename U> class DynamicStorage;
    template <typename T, int N> class StaticFixedStorage;
    template <typename T, int N> class StaticStorage;
    template <typename S> class GeneralFad;
  }
  }
#ifndef SACADO_VIEW_CUDA_HIERARCHICAL_DFAD
  template <unsigned Stride, typename T, typename U>
  KOKKOS_INLINE_FUNCTION
  typename LocalScalarType< Fad::Exp::GeneralFad< Fad::Exp::DynamicStorage<T,U> >, Stride >::type
  partition_scalar(const Fad::Exp::GeneralFad< Fad::Exp::DynamicStorage<T,U> >& x) {
    typedef typename LocalScalarType< Fad::Exp::GeneralFad< Fad::Exp::DynamicStorage<T,U> >, Stride >::type ret_type;
    const int size = (x.size()+blockDim.x-threadIdx.x-1) / blockDim.x;
    const int offset = threadIdx.x;
    ret_type xp(size, x.val());

    // Note:  we can't use x.dx(offset+i*Stride) if
    // SACADO_VIEW_CUDA_HIERARCHICAL_DFAD_STRIDED is defined because it already
    // uses blockDim.x in its index calculation.  This approach should work
    // regardless
    const T* dx = x.dx();
    for (int i=0; i<size; ++i)
      xp.fastAccessDx(i) = dx[offset+i*Stride];

    return xp;
  }
#endif
  template <unsigned Stride, typename T, int N>
  KOKKOS_INLINE_FUNCTION
  typename LocalScalarType< Fad::Exp::GeneralFad< Fad::Exp::StaticStorage<T,N> >, Stride >::type
  partition_scalar(const Fad::Exp::GeneralFad< Fad::Exp::StaticStorage<T,N> >& x) {
    typedef typename LocalScalarType< Fad::Exp::GeneralFad< Fad::Exp::StaticStorage<T,N> >, Stride >::type ret_type;
    const int size = (x.size()+blockDim.x-threadIdx.x-1) / blockDim.x;
    const int offset = threadIdx.x;
    ret_type xp(size, x.val());
    for (int i=0; i<size; ++i)
      xp.fastAccessDx(i) = x.fastAccessDx(offset+i*Stride);
    return xp;
  }
  template <unsigned Stride, typename T, int N>
  KOKKOS_INLINE_FUNCTION
  typename LocalScalarType< Fad::Exp::GeneralFad< Fad::Exp::StaticFixedStorage<T,N> >, Stride >::type
  partition_scalar(const Fad::Exp::GeneralFad< Fad::Exp::StaticFixedStorage<T,N> >& x) {
    typedef typename LocalScalarType< Fad::Exp::GeneralFad< Fad::Exp::StaticFixedStorage<T,N> >, Stride >::type ret_type;
    const int size = (x.size()+blockDim.x-threadIdx.x-1) / blockDim.x;
    const int offset = threadIdx.x;
    ret_type xp(size, x.val());
    for (int i=0; i<size; ++i)
      xp.fastAccessDx(i) = x.fastAccessDx(offset+i*Stride);
    return xp;
  }

#ifndef SACADO_NEW_FAD_DESIGN_IS_DEFAULT
  namespace Fad {
    template <typename T> class DFad;
    template <typename T, int N> class SLFad;
    template <typename T, int N> class SFad;
  }
#ifndef SACADO_VIEW_CUDA_HIERARCHICAL_DFAD
  template <unsigned Stride, typename T>
  KOKKOS_INLINE_FUNCTION
  typename LocalScalarType< Fad::DFad<T>, Stride >::type
  partition_scalar(const Fad::DFad<T>& x) {
    typedef typename LocalScalarType< Fad::DFad<T>, Stride >::type ret_type;
    const int size = (x.size()+blockDim.x-threadIdx.x-1) / blockDim.x;
    const int offset = threadIdx.x;
    ret_type xp(size, x.val());

    // Note:  we can't use x.dx(offset+i*Stride) if
    // SACADO_VIEW_CUDA_HIERARCHICAL_DFAD_STRIDED is defined because it already
    // uses blockDim.x in its index calculation.  This approach should work
    // regardless
    const T* dx = x.dx();
    for (int i=0; i<size; ++i)
      xp.fastAccessDx(i) = dx[offset+i*Stride];

    return xp;
  }
#endif
  template <unsigned Stride, typename T, int N>
  KOKKOS_INLINE_FUNCTION
  typename LocalScalarType< Fad::SLFad<T,N>, Stride >::type
  partition_scalar(const Fad::SLFad<T,N>& x) {
    typedef typename LocalScalarType< Fad::SLFad<T,N>, Stride >::type ret_type;
    const int size = (x.size()+blockDim.x-threadIdx.x-1) / blockDim.x;
    const int offset = threadIdx.x;
    ret_type xp(size, x.val());
    for (int i=0; i<size; ++i)
      xp.fastAccessDx(i) = x.fastAccessDx(offset+i*Stride);
    return xp;
  }
  template <unsigned Stride, typename T, int N>
  KOKKOS_INLINE_FUNCTION
  typename LocalScalarType< Fad::SFad<T,N>, Stride >::type
  partition_scalar(const Fad::SFad<T,N>& x) {
    typedef typename LocalScalarType< Fad::SFad<T,N>, Stride >::type ret_type;
    const int size = (x.size()+blockDim.x-threadIdx.x-1) / blockDim.x;
    const int offset = threadIdx.x;
    ret_type xp(size, x.val());
    for (int i=0; i<size; ++i)
      xp.fastAccessDx(i) = x.fastAccessDx(offset+i*Stride);
    return xp;
  }
#endif // SACADO_NEW_FAD_DESIGN_IS_DEFAULT

#endif // defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)

} // namespace Sacado

//----------------------------------------------------------------------------

namespace Kokkos {

template< unsigned Stride, typename D, typename ... P  >
KOKKOS_INLINE_FUNCTION
typename Kokkos::Impl::ViewMapping< void, typename Kokkos::ViewTraits<D,P...>, Sacado::Fad::Partition<Stride> >::type
partition( const Kokkos::View<D,P...> & src ,
           const unsigned offset ,
           const unsigned stride )
{
  typedef Kokkos::ViewTraits<D,P...> traits;
  typedef typename Kokkos::Impl::ViewMapping< void, traits, Sacado::Fad::Partition<Stride> >::type DstViewType;
  const Sacado::Fad::Partition<Stride> part( offset , stride );
  return DstViewType(src, part);
}

template <typename ViewType>
struct ThreadLocalScalarType<
  ViewType,
  typename std::enable_if< is_view_fad_contiguous<ViewType>::value >::type > {
  typedef typename ViewType::traits TraitsType;
  typedef Impl::ViewMapping<TraitsType, typename TraitsType::specialize> MappingType;
  typedef typename MappingType::thread_local_scalar_type type;
};

namespace Impl {

#if defined (KOKKOS_ENABLE_CUDA) && defined(SACADO_VIEW_CUDA_HIERARCHICAL)
template< class OutputView >
struct SacadoViewFill<
  OutputView,
  typename std::enable_if<
    ( Kokkos::is_view_fad_contiguous<OutputView>::value &&
      std::is_same<typename OutputView::execution_space, Kokkos::Cuda>::value &&
      !Kokkos::ViewScalarStride<OutputView>::is_unit_stride )
    >::type
  >
{
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename OutputView::execution_space execution_space ;
  typedef Kokkos::TeamPolicy< execution_space> team_policy;
  typedef typename team_policy::member_type team_impl_handle;
  typedef typename Kokkos::ThreadLocalScalarType<OutputView>::type local_scalar_type;
  static const unsigned stride = Kokkos::ViewScalarStride<OutputView>::stride;

  const OutputView output ;
  const_value_type input ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_t i0 ) const
  {
    local_scalar_type input_stride = Sacado::partition_scalar<stride>(input);

    const size_t n1 = output.extent(1);
    const size_t n2 = output.extent(2);
    const size_t n3 = output.extent(3);
    const size_t n4 = output.extent(4);
    const size_t n5 = output.extent(5);
    const size_t n6 = output.extent(6);

    for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_t i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_t i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_t i6 = 0 ; i6 < n6 ; ++i6 ) {
      output.access(i0,i1,i2,i3,i4,i5,i6) = input_stride ;
    }}}}}}
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_impl_handle& team ) const
  {
    const size_t i0 = team.league_rank()*team.team_size() + team.team_rank();
    if (i0 < output.extent(0))
      (*this)(i0);
  }

  SacadoViewFill( const OutputView & arg_out , const_value_type & arg_in )
    : output( arg_out ), input( arg_in )
    {
      const size_t team_size = 256 / stride;
      team_policy policy( (output.extent(0)+team_size-1)/team_size ,
                          team_size , stride );
      Kokkos::parallel_for( policy, *this );
    }
};
#endif

#if defined (KOKKOS_ENABLE_HIP) && defined(SACADO_VIEW_CUDA_HIERARCHICAL)
template< class OutputView >
struct SacadoViewFill<
  OutputView,
  typename std::enable_if<
    ( Kokkos::is_view_fad_contiguous<OutputView>::value &&
      std::is_same<typename OutputView::execution_space, Kokkos::HIP>::value &&
      !Kokkos::ViewScalarStride<OutputView>::is_unit_stride )
    >::type
  >
{
  typedef typename OutputView::const_value_type  const_value_type ;
  typedef typename OutputView::execution_space execution_space ;
  typedef Kokkos::TeamPolicy< execution_space> team_policy;
  typedef typename team_policy::member_type team_impl_handle;
  typedef typename Kokkos::ThreadLocalScalarType<OutputView>::type local_scalar_type;
  static const unsigned stride = Kokkos::ViewScalarStride<OutputView>::stride;

  const OutputView output ;
  const_value_type input ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_t i0 ) const
  {
    local_scalar_type input_stride = Sacado::partition_scalar<stride>(input);

    const size_t n1 = output.extent(1);
    const size_t n2 = output.extent(2);
    const size_t n3 = output.extent(3);
    const size_t n4 = output.extent(4);
    const size_t n5 = output.extent(5);
    const size_t n6 = output.extent(6);
    const size_t n7 = output.extent(7);

    for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
    for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
    for ( size_t i3 = 0 ; i3 < n3 ; ++i3 ) {
    for ( size_t i4 = 0 ; i4 < n4 ; ++i4 ) {
    for ( size_t i5 = 0 ; i5 < n5 ; ++i5 ) {
    for ( size_t i6 = 0 ; i6 < n6 ; ++i6 ) {
      output.access(i0,i1,i2,i3,i4,i5,i6) = input_stride ;
    }}}}}}
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_impl_handle& team ) const
  {
    const size_t i0 = team.league_rank()*team.team_size() + team.team_rank();
    if (i0 < output.extent(0))
      (*this)(i0);
  }

  SacadoViewFill( const OutputView & arg_out , const_value_type & arg_in )
    : output( arg_out ), input( arg_in )
    {
      const size_t team_size = 256 / stride;
      team_policy policy( (output.extent(0)+team_size-1)/team_size ,
                          team_size , stride );
      Kokkos::parallel_for( policy, *this );
    }
};
#endif

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class ... Args >
struct is_ViewSpecializeSacadoFadContiguous { enum { value = false }; };

template< class D , class ... P , class ... Args >
struct is_ViewSpecializeSacadoFadContiguous< Kokkos::View<D,P...> , Args... > {
  enum { value =
    std::is_same< typename Kokkos::ViewTraits<D,P...>::specialize
                , ViewSpecializeSacadoFadContiguous >::value
    &&
    ( ( sizeof...(Args) == 0 ) ||
      is_ViewSpecializeSacadoFadContiguous< Args... >::value ) };
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Compute a partitioned fad size given a stride.  Return 0 if the stride
// does not evenly divide the size
KOKKOS_INLINE_FUNCTION
constexpr unsigned computeFadPartitionSize(unsigned size, unsigned stride)
{
  return
    ((size+stride-1)/stride) == (size/stride) ? ((size+stride-1)/stride) : 0;
}

// Create new Layout with Fad dimension set to the last
// Note:  we use enable_if<> here to handle both LayoutLeft and
// LayoutContiguous<LayoutLeft>
template <unsigned rank, unsigned static_dim, typename Layout>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< !std::is_same<Layout, LayoutLeft>::value &&
                         !std::is_same<Layout, LayoutStride>::value,
                       Layout>::type
create_fad_array_layout(const Layout& layout)
{
  size_t dims[8];
  for (int i=0; i<8; ++i)
    dims[i] = layout.dimension[i];
  if (static_dim > 0)
    dims[rank] = static_dim+1;
  return Layout( dims[0], dims[1], dims[2], dims[3],
                 dims[4], dims[5], dims[6], dims[7] );
}

// Create new Layout with Fad dimension set to the last
// Note:  we use enable_if<> here to handle both LayoutStride and
// LayoutContiguous<LayoutStride>
template <unsigned rank, unsigned static_dim, typename Layout>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< std::is_same<Layout, LayoutStride>::value, Layout>::type
create_fad_array_layout(const Layout& layout)
{
  size_t dims[8], strides[8];
  for (int i=0; i<8; ++i) {
    dims[i] = layout.dimension[i];
    strides[i] = layout.stride[i];
  }
  if (static_dim > 0) {
    dims[rank] = static_dim+1;
    strides[rank] = 1;
  }
  return Layout( dims[0], strides[0],
                 dims[1], strides[1],
                 dims[2], strides[2],
                 dims[3], strides[3],
                 dims[4], strides[4],
                 dims[5], strides[5],
                 dims[6], strides[6],
                 dims[7], strides[7] );
}

// Create new LayoutLeft with Fad dimension shuffled to the first
// Note:  we use enable_if<> here to handle both LayoutLeft and
// LayoutContiguous<LayoutLeft>
  template <unsigned rank, unsigned static_dim, typename Layout>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< std::is_same<Layout, LayoutLeft>::value, Layout >::type
create_fad_array_layout(const Layout& layout)
{
  size_t dims[8];
  for (int i=0; i<8; ++i)
    dims[i] = layout.dimension[i];
  size_t fad_dim = static_dim == 0 ? dims[rank] : static_dim+1;
  for (int i=rank; i>=1; --i)
    dims[i] = dims[i-1];
  dims[0] = fad_dim;
  return Layout( dims[0], dims[1], dims[2], dims[3],
                 dims[4], dims[5], dims[6], dims[7] );
}

template <unsigned Rank, typename Dimension, typename Layout>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< !std::is_same<Layout, LayoutLeft>::value, size_t>::type
getFadDimension(const ViewOffset<Dimension,Layout,void>& offset)
{
  return
    ( Rank == 0 ? offset.dimension_0() :
    ( Rank == 1 ? offset.dimension_1() :
    ( Rank == 2 ? offset.dimension_2() :
    ( Rank == 3 ? offset.dimension_3() :
    ( Rank == 4 ? offset.dimension_4() :
    ( Rank == 5 ? offset.dimension_5() :
    ( Rank == 6 ? offset.dimension_6() :
      offset.dimension_7() )))))));
}

template <unsigned Rank, typename Dimension, typename Layout>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< std::is_same<Layout, LayoutLeft>::value, size_t >::type
getFadDimension(const ViewOffset<Dimension,Layout,void>& offset)
{
  return offset.dimension_0();
}

template< class Traits >
class ViewMapping< Traits , /* View internal mapping */
  typename std::enable_if<
    ( std::is_same< typename Traits::specialize
                  , ViewSpecializeSacadoFadContiguous >::value
      &&
      ( std::is_same< typename Traits::array_layout
                    , Kokkos::LayoutLeft >::value
        ||
        std::is_same< typename Traits::array_layout
                    , Kokkos::LayoutRight >::value
        ||
        std::is_same< typename Traits::array_layout
                    , Kokkos::LayoutStride >::value
      )
    )
    , typename Traits::specialize
    >::type >
{
private:

  template< class , class ... > friend class ViewMapping ;
  template< class , class ... > friend class Kokkos::View ;

  typedef typename Traits::value_type  fad_type ;
  typedef typename Sacado::ValueType< fad_type >::type fad_value_type ;
  typedef typename
    std::add_const< fad_value_type >::type  const_fad_value_type ;

public:
  enum { is_assignable_data_type = true };

  enum { FadStaticDimension = Sacado::StaticSize< fad_type >::value };
  enum { PartitionedFadStride = Traits::array_layout::scalar_stride };

  // The partitioned static size -- this will be 0 if ParitionedFadStride
  // does not evenly divide FadStaticDimension
  enum { PartitionedFadStaticDimension =
           computeFadPartitionSize(FadStaticDimension,PartitionedFadStride) };

#ifdef KOKKOS_ENABLE_CUDA
  typedef typename Sacado::LocalScalarType< fad_type, unsigned(PartitionedFadStride) >::type strided_scalar_type;
  typedef typename std::conditional< std::is_same<typename Traits::execution_space, Kokkos::Cuda>::value, strided_scalar_type, fad_type >::type thread_local_scalar_type;
#elif defined(KOKKOS_ENABLE_HIP)
  typedef typename Sacado::LocalScalarType< fad_type, unsigned(PartitionedFadStride) >::type strided_scalar_type;
  typedef typename std::conditional< std::is_same<typename Traits::execution_space, Kokkos::HIP>::value, strided_scalar_type, fad_type >::type thread_local_scalar_type;
#else
  typedef fad_type thread_local_scalar_type;
#endif

private:
  typedef Sacado::integral_nonzero< unsigned , FadStaticDimension > sacado_size_type;

  typedef fad_value_type * handle_type ;

  typedef ViewArrayAnalysis< typename Traits::data_type > array_analysis ;

  typedef ViewOffset< typename Traits::dimension
                    , typename Traits::array_layout
                    , void
                    >  offset_type ;

  // Prepend/append the fad dimension for the internal offset mapping.
  static constexpr bool is_layout_left =
    std::is_same< typename Traits::array_layout, Kokkos::LayoutLeft>::value;
  typedef ViewOffset<
    typename std::conditional<
      is_layout_left,
      typename array_analysis::dimension::
        template prepend<0>::type,
      typename array_analysis::dimension::
        template append<0>::type >::type,
      typename Traits::array_layout,
      void >
    array_offset_type ;

  handle_type  m_impl_handle ;
  offset_type  m_impl_offset ;
  array_offset_type  m_array_offset ;
  sacado_size_type m_fad_size ;

  // These are for manual partitioning, and will likely be removed
  unsigned m_original_fad_size ;
  unsigned m_fad_stride ;
  unsigned m_fad_index ;

public:

  //----------------------------------------
  // Domain dimensions

  enum { Rank = Traits::dimension::rank };

  // Rank corresponding to the sacado dimension
  enum { Sacado_Rank = std::is_same< typename Traits::array_layout, Kokkos::LayoutLeft >::value ? 0 : Rank+1 };

  // Using the internal offset mapping so limit to public rank:
  template< typename iType >
  KOKKOS_INLINE_FUNCTION constexpr size_t extent( const iType & r ) const
    { return m_impl_offset.m_dim.extent(r); }

  KOKKOS_INLINE_FUNCTION constexpr
  typename Traits::array_layout layout() const
    { return m_impl_offset.layout(); }

  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const
    { return m_impl_offset.dimension_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const
    { return m_impl_offset.dimension_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const
    { return m_impl_offset.dimension_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const
    { return m_impl_offset.dimension_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const
    { return m_impl_offset.dimension_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const
    { return m_impl_offset.dimension_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const
    { return m_impl_offset.dimension_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const
    { return m_impl_offset.dimension_7(); }

  // Is a regular layout with uniform striding for each index.
  // Since we allow for striding within the data type, we can't guarantee
  // regular striding
  using is_regular = std::false_type ;

  // FIXME:  Adjust these for m_stride
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const
    { return m_impl_offset.stride_0(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const
    { return m_impl_offset.stride_1(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const
    { return m_impl_offset.stride_2(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const
    { return m_impl_offset.stride_3(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const
    { return m_impl_offset.stride_4(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const
    { return m_impl_offset.stride_5(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const
    { return m_impl_offset.stride_6(); }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const
    { return m_impl_offset.stride_7(); }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION void stride( iType * const s ) const
    { m_impl_offset.stride(s); }

  // Size of sacado scalar dimension
  KOKKOS_FORCEINLINE_FUNCTION constexpr unsigned dimension_scalar() const
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
    { return PartitionedFadStaticDimension ? PartitionedFadStaticDimension+1 : (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x + 1; }
#else
    { return m_fad_size.value+1; }
#endif

  // trode of sacado scalar dimension
  KOKKOS_FORCEINLINE_FUNCTION constexpr unsigned stride_scalar() const
    { return m_fad_stride; }

  //----------------------------------------
  // Range of mapping

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
  // Return type of reference operators
  // this only works if you are using a team-parallel operation on Cuda or HIP!
  // typedef typename
  //   Sacado::ViewFadType< thread_local_scalar_type , PartitionedFadStaticDimension , (unsigned(ParitionedFadStride) > 1 ? PartitionedFadStride : 0) >::type  reference_type ;
  typedef typename
    Sacado::ViewFadType< thread_local_scalar_type , PartitionedFadStaticDimension , 0 >::type  reference_type ;
#else
  // Return type of reference operators
  typedef typename
    Sacado::ViewFadType< fad_type , FadStaticDimension , 0 >::type  reference_type ;
#endif

  /** \brief Pointer to underlying memory type */
  typedef fad_value_type * pointer_type ;

  /** \brief  Span of the mapped range : [ data() .. data() + span() ) */
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const
    { return m_array_offset.span(); }

  /** \brief  Is the mapped range span contiguous */
  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const
    { return m_array_offset.span_is_contiguous() && (m_fad_stride == 1); }

  /** \brief Raw data access */
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
    { return m_impl_handle + threadIdx.x; }
#else
    { return m_impl_handle + m_fad_index; }
#endif

  //----------------------------------------

  KOKKOS_FORCEINLINE_FUNCTION
  reference_type
  reference() const
    {
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( m_impl_handle + index
                           , m_impl_handle + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<  Kokkos::Impl::are_integral<I0>::value &&
                            is_layout_left, reference_type>::type
  reference( const I0 & i0 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(0,i0);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0>::value &&
                           !is_layout_left, reference_type>::type
  reference( const I0 & i0 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(i0,0);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1>::value &&
                           is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(0,i0,i1);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1>::value &&
                           !is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(i0,i1,0);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }


  template< typename I0 , typename I1 , typename I2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2>::value &&
                           is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(0,i0,i1,i2);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 , typename I2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2>::value &&
                           !is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(i0,i1,i2,0);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 , typename I2 , typename I3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2,I3>::value &&
                           is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(0,i0,i1,i2,i3);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 , typename I2 , typename I3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2,I3>::value &&
                           !is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(i0,i1,i2,i3,0);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2,I3,I4>::value &&
                           is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
           , const I4 & i4 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(0,i0,i1,i2,i3,i4);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2,I3,I4>::value &&
                           !is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
           , const I4 & i4 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(i0,i1,i2,i3,i4,0);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5>::value &&
                           is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
           , const I4 & i4 , const I5 & i5 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(0,i0,i1,i2,i3,i4,i5);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5>::value &&
                           !is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
           , const I4 & i4 , const I5 & i5 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(i0,i1,i2,i3,i4,i5,0);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,I6>::value &&
                           is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
           , const I4 & i4 , const I5 & i5 , const I6 & i6 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(0,i0,i1,i2,i3,i4,i5,i6);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  template< typename I0 , typename I1 , typename I2 , typename I3
          , typename I4 , typename I5 , typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if< Kokkos::Impl::are_integral<I0,I1,I2,I3,I4,I5,I6>::value &&
                           !is_layout_left, reference_type>::type
  reference( const I0 & i0 , const I1 & i1 , const I2 & i2 , const I3 & i3
           , const I4 & i4 , const I5 & i5 , const I6 & i6 ) const
    { pointer_type beg = m_impl_handle + m_array_offset(i0,i1,i2,i3,i4,i5,i6,0);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) && ( defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) )
      const unsigned index = threadIdx.x;
      const unsigned strd = blockDim.x;
      const unsigned size = (m_fad_size.value+blockDim.x-threadIdx.x-1) / blockDim.x;
#else
      const unsigned index = m_fad_index;
      const unsigned strd = m_fad_stride;
      const unsigned size = m_fad_size.value;
#endif
      return reference_type( beg + index
                           , beg + m_original_fad_size
                           , size
                           , strd ); }

  //----------------------------------------

  /** \brief  Span, in bytes, of the required memory */
  KOKKOS_INLINE_FUNCTION
  static constexpr size_t memory_span( typename Traits::array_layout const & layout )
    {
      // Do not introduce padding...
      typedef std::integral_constant< unsigned , 0 >  padding ;
      return array_offset_type(
        padding() ,
        create_fad_array_layout<unsigned(Rank), unsigned(FadStaticDimension)>( layout ) ).span() * sizeof(fad_value_type);
    }

  //----------------------------------------

  KOKKOS_DEFAULTED_FUNCTION ~ViewMapping() = default ;
  KOKKOS_INLINE_FUNCTION ViewMapping() : m_impl_handle(0) , m_impl_offset() , m_array_offset() , m_fad_size(0) , m_original_fad_size(0) , m_fad_stride(1) , m_fad_index(0)  {}

  KOKKOS_DEFAULTED_FUNCTION ViewMapping( const ViewMapping & ) = default ;
  KOKKOS_DEFAULTED_FUNCTION ViewMapping & operator = ( const ViewMapping & ) = default ;

  KOKKOS_DEFAULTED_FUNCTION ViewMapping( ViewMapping && ) = default ;
  KOKKOS_DEFAULTED_FUNCTION ViewMapping & operator = ( ViewMapping && ) = default ;

  template< class ... P >
  KOKKOS_INLINE_FUNCTION
  ViewMapping
    ( ViewCtorProp< P ... > const & prop
    , typename Traits::array_layout const & local_layout
    )
    : m_impl_handle( ( (ViewCtorProp<void,pointer_type> const &) prop ).value )
    , m_impl_offset( std::integral_constant< unsigned , 0 >()
              , local_layout )
    , m_array_offset(
        std::integral_constant< unsigned , 0 >()
        , create_fad_array_layout<unsigned(Rank), unsigned(FadStaticDimension)>( local_layout ) )
    , m_fad_size( getFadDimension<unsigned(Rank)>( m_array_offset ) - 1 )
    , m_original_fad_size( m_fad_size.value )
    , m_fad_stride( 1 )
    , m_fad_index( 0 )
    {
      const unsigned fad_dim =
        getFadDimension<unsigned(Rank)>( m_array_offset );
      if (unsigned(FadStaticDimension) == 0 && fad_dim == 0)
        Kokkos::abort("invalid fad dimension (0) supplied!");
    }

  //----------------------------------------
  /*  Allocate and construct mapped array.
   *  Allocate via shared allocation record and
   *  return that record for allocation tracking.
   */
  template< class ... P >
  SharedAllocationRecord<> *
  allocate_shared( ViewCtorProp< P... > const & prop
                 , typename Traits::array_layout const & local_layout
                 , bool execution_space_specified)
  {
    typedef ViewCtorProp< P... > ctor_prop ;

    typedef typename ctor_prop::execution_space  execution_space ;
    typedef typename Traits::memory_space         memory_space ;
    typedef ViewValueFunctor< execution_space , fad_value_type > functor_type ;
    typedef SharedAllocationRecord< memory_space , functor_type > record_type ;

    // Disallow padding
    typedef std::integral_constant< unsigned , 0 > padding ;

    // Check if ViewCtorProp has CommonViewAllocProp - if so, retrieve the fad_size and append to layout
    enum { test_traits_check = Kokkos::Impl::check_has_common_view_alloc_prop< P... >::value };

    typename Traits::array_layout internal_layout =
      (test_traits_check == true)
      ? Kokkos::Impl::appendFadToLayoutViewAllocHelper< Traits, P... >::returnNewLayoutPlusFad(prop, local_layout)
      : local_layout;

    m_impl_offset = offset_type( padding(), internal_layout );

    m_array_offset =
      array_offset_type( padding() ,
                         create_fad_array_layout<unsigned(Rank), unsigned(FadStaticDimension)>( internal_layout ) );
    const unsigned fad_dim =
      getFadDimension<unsigned(Rank)>( m_array_offset );
    if (unsigned(FadStaticDimension) == 0 && fad_dim == 0)
      Kokkos::abort("invalid fad dimension (0) supplied!");
    m_fad_size = fad_dim - 1 ;
    m_original_fad_size = m_fad_size.value ;
    m_fad_stride = 1;
    m_fad_index = 0;

    const size_t alloc_size = m_array_offset.span() * sizeof(fad_value_type);

    // Create shared memory tracking record with allocate memory from the memory space
    record_type * const record =
      record_type::allocate( ( (ViewCtorProp<void,memory_space> const &) prop ).value
                           , ( (ViewCtorProp<void,std::string>  const &) prop ).value
                           , alloc_size );

    //  Only set the the pointer and initialize if the allocation is non-zero.
    //  May be zero if one of the dimensions is zero.
    if ( alloc_size ) {

      m_impl_handle = handle_type( reinterpret_cast< pointer_type >( record->data() ) );

      if ( ctor_prop::initialize ) {
        // Assume destruction is only required when construction is requested.
        // The ViewValueFunctor has both value construction and destruction operators.
        if (execution_space_specified)
          record->m_destroy = functor_type( ( (ViewCtorProp<void,execution_space> const &) prop).value
                                          , (fad_value_type *) m_impl_handle
                                          , m_array_offset.span()
                                          , record->get_label()
                                          );
        else
          record->m_destroy = functor_type((fad_value_type *) m_impl_handle
                                          , m_array_offset.span()
                                          , record->get_label()
                                          );

        // Construct values
        record->m_destroy.construct_shared_allocation();
      }
    }

    return record ;
  }

};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/**\brief  Assign compatible Sacado FAD view mappings.
 *
 *  View<FAD>      = View<FAD>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    Kokkos::Impl::MemorySpaceAccess
     < typename DstTraits::memory_space
     , typename SrcTraits::memory_space >::assignable
    &&
    // Destination view has FAD
    std::is_same< typename DstTraits::specialize
                , ViewSpecializeSacadoFadContiguous >::value
    &&
    // Source view has FAD
    std::is_same< typename SrcTraits::specialize
                , ViewSpecializeSacadoFadContiguous >::value
  )
  , typename DstTraits::specialize
  >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcFadType ;

  template< class DstType >
  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcFadType & src
             , const TrackType & )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of FAD requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::array_layout
                    , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
                    , Kokkos::LayoutStride >::value ,
        "View assignment must have compatible layout" );

      static_assert(
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::const_value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      static_assert(
        ViewDimensionAssignable
          < typename DstType::offset_type::dimension_type
          , typename SrcFadType::offset_type::dimension_type >::value ,
        "View assignment must have compatible dimensions" );

       static_assert(
        ViewDimensionAssignable
          < typename DstType::array_offset_type::dimension_type
          , typename SrcFadType::array_offset_type::dimension_type >::value ,
        "View assignment must have compatible dimensions" );

      typedef typename DstType::offset_type  dst_offset_type ;
      typedef typename DstType::array_offset_type  dst_array_offset_type ;

      dst.m_impl_handle  = src.m_impl_handle ;
      dst.m_impl_offset  = dst_offset_type( src.m_impl_offset );
      dst.m_array_offset = dst_array_offset_type( src.m_array_offset );
      dst.m_fad_size = src.m_fad_size.value ;
      dst.m_original_fad_size = src.m_original_fad_size ;
      dst.m_fad_stride = src.m_fad_stride ;
      dst.m_fad_index = src.m_fad_index ;
    }
};

/**\brief  Assign compatible Sacado FAD view mappings.
 *
 *  View<FAD,LayoutStride>      = View<FAD,LayoutContiguous>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    std::is_same< typename DstTraits::memory_space
                , typename SrcTraits::memory_space >::value
    &&
    // Destination view has FAD
    std::is_same< typename DstTraits::specialize
                , ViewSpecializeSacadoFad >::value
    &&
    // Source view has FAD contiguous
    std::is_same< typename SrcTraits::specialize
                , ViewSpecializeSacadoFadContiguous >::value
    &&
    // Destination view is LayoutStride
    std::is_same< typename DstTraits::array_layout
                , Kokkos::LayoutStride >::value
  )
  , typename DstTraits::specialize
  >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcFadType ;

  template< class DstType >
  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcFadType & src
             , const TrackType & )
    {
      static_assert(
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutLeft >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutRight >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutStride >::value ,
        "View of FAD requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::value_type >::value ||
        std::is_same< typename DstTraits::value_type
                    , typename SrcTraits::const_value_type >::value ,
        "View assignment must have same value type or const = non-const" );

      static_assert(
        DstTraits::dimension::rank == SrcTraits::dimension::rank,
        "View assignment must have same rank" );

      typedef typename DstType::array_offset_type  dst_offset_type ;

      dst.m_impl_handle  = src.m_impl_handle ;
      dst.m_fad_size = src.m_fad_size.value ;
      dst.m_fad_stride = src.m_fad_stride ;
      dst.m_impl_offset = src.m_impl_offset;

      size_t N[8], S[8];
      N[0] = src.m_array_offset.dimension_0();
      N[1] = src.m_array_offset.dimension_1();
      N[2] = src.m_array_offset.dimension_2();
      N[3] = src.m_array_offset.dimension_3();
      N[4] = src.m_array_offset.dimension_4();
      N[5] = src.m_array_offset.dimension_5();
      N[6] = src.m_array_offset.dimension_6();
      N[7] = src.m_array_offset.dimension_7();
      S[0] = src.m_array_offset.stride_0();
      S[1] = src.m_array_offset.stride_1();
      S[2] = src.m_array_offset.stride_2();
      S[3] = src.m_array_offset.stride_3();
      S[4] = src.m_array_offset.stride_4();
      S[5] = src.m_array_offset.stride_5();
      S[6] = src.m_array_offset.stride_6();
      S[7] = src.m_array_offset.stride_7();

      // For LayoutLeft, we have to move the Sacado dimension from the first
      // to the last
      if (std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value)
      {
        const size_t N_fad = N[0];
        const size_t S_fad = S[0];
        for (int i=0; i<7; ++i) {
          N[i] = N[i+1];
          S[i] = S[i+1];
        }
        N[DstTraits::dimension::rank] = N_fad;
        S[DstTraits::dimension::rank] = S_fad;
      }
      Kokkos::LayoutStride ls( N[0], S[0],
                               N[1], S[1],
                               N[2], S[2],
                               N[3], S[3],
                               N[4], S[4],
                               N[5], S[5],
                               N[6], S[6],
                               N[7], S[7] );
      dst.m_array_offset  = dst_offset_type(std::integral_constant<unsigned,0>(), ls);
    }
};

/**\brief  Assign compatible Sacado FAD view mappings.
 *
 *  View<ordinary> = View<FAD>
 */
template< class DstTraits , class SrcTraits >
class ViewMapping< DstTraits , SrcTraits ,
  typename std::enable_if<(
    std::is_same< typename DstTraits::memory_space
                , typename SrcTraits::memory_space >::value
    &&
    // Destination view has ordinary
    std::is_same< typename DstTraits::specialize , void >::value
    &&
    // Source view has FAD only
    std::is_same< typename SrcTraits::specialize
                , ViewSpecializeSacadoFadContiguous >::value
  )
  , typename DstTraits::specialize
  >::type >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef Kokkos::Impl::SharedAllocationTracker  TrackType ;
  typedef ViewMapping< DstTraits , typename DstTraits::specialize >  DstType ;
  typedef ViewMapping< SrcTraits , typename SrcTraits::specialize >  SrcFadType ;


  // Helpers to assign, and generate if necessary, ViewOffset to the dst map
  // These are necessary to use Kokkos' deep_copy with nested fads
  template < class DstType, class SrcFadType, class Enable = void >
    struct AssignOffset;

  template < class DstType, class SrcFadType >
    struct AssignOffset< DstType, SrcFadType, typename std::enable_if< ((int)DstType::offset_type::dimension_type::rank != (int)SrcFadType::array_offset_type::dimension_type::rank) >::type >
    {
      // ViewOffset's Dimensions Ranks do not match
      KOKKOS_INLINE_FUNCTION
      static void assign( DstType & dst, const SrcFadType & src )
      {
        typedef typename SrcTraits::value_type TraitsValueType;

        if ( Sacado::IsFad<TraitsValueType>::value
            && Sacado::IsStaticallySized< typename Sacado::ValueType< TraitsValueType >::type >::value
           )
        {

          typedef typename DstType::offset_type::array_layout DstLayoutType;
          //typedef typename ViewArrayLayoutSelector<typename DstType::offset_type::array_layout>::type DstLayoutType;
          typedef typename SrcFadType::array_offset_type::dimension_type SrcViewDimension;

          // This is the static dimension of the inner fad, missing from ViewDimension
          const size_t InnerStaticDim = Sacado::StaticSize< typename Sacado::ValueType< TraitsValueType >::type >::value;

          static constexpr bool is_layout_left =
            std::is_same< DstLayoutType, Kokkos::LayoutLeft>::value;

          typedef typename std::conditional< is_layout_left,
                                             typename SrcViewDimension:: template prepend< InnerStaticDim+1 >::type,
                                             typename SrcViewDimension:: template append < InnerStaticDim+1 >::type
                    >::type SrcViewDimensionAppended;

          typedef std::integral_constant< unsigned , 0 >  padding ;

          typedef ViewOffset< SrcViewDimensionAppended, DstLayoutType > TmpOffsetType;

          auto src_layout = src.m_array_offset.layout();

          if ( is_layout_left ) {
            auto prepend_layout = Kokkos::Impl::prependFadToLayout< DstLayoutType >::returnNewLayoutPlusFad(src_layout, InnerStaticDim+1);
            TmpOffsetType offset_tmp( padding(), prepend_layout );
            dst.m_impl_offset = offset_tmp;
          }
          else {
            TmpOffsetType offset_tmp( padding(), src_layout );
            dst.m_impl_offset = offset_tmp;
          }
        } else {
          Kokkos::abort("Sacado error: Applying AssignOffset for case with nested Fads, but without nested Fads - something went wrong");
        }
      }
    };

  template < class DstType, class SrcFadType >
    struct AssignOffset< DstType, SrcFadType, typename std::enable_if< ((int)DstType::offset_type::dimension_type::rank == (int)SrcFadType::array_offset_type::dimension_type::rank) >::type >
    {
      KOKKOS_INLINE_FUNCTION
      static void assign( DstType & dst, const SrcFadType & src )
      {
        typedef typename DstType::offset_type  dst_offset_type ;
        dst.m_impl_offset  = dst_offset_type( src.m_array_offset );
      }
    };

  template< class DstType >
  KOKKOS_INLINE_FUNCTION static
  void assign( DstType & dst
             , const SrcFadType & src
             , const TrackType &
             )
    {
      static_assert(
        (
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename DstTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        &&
        (
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutLeft >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutRight >::value ||
          std::is_same< typename SrcTraits::array_layout
                      , Kokkos::LayoutStride >::value
        )
        , "View of FAD requires LayoutLeft, LayoutRight, or LayoutStride" );

      static_assert(
        std::is_same< typename DstTraits::array_layout
                    , typename SrcTraits::array_layout >::value ||
        std::is_same< typename DstTraits::array_layout
                    , Kokkos::LayoutStride >::value ,
        "View assignment must have compatible layout" );

      if ( src.m_fad_index != 0 || src.m_fad_stride != 1 ) {
        Kokkos::abort("\n\n ******  Kokkos::View< Sacado::Fad ... > Cannot assign to array with partitioned view ******\n\n");
      }

      AssignOffset< DstType, SrcFadType >::assign( dst, src );
      dst.m_impl_handle  = reinterpret_cast< typename DstType::handle_type >(src.m_impl_handle) ;
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Rules for subview arguments and layouts matching

template<class LayoutDest, class LayoutSrc, int RankDest, int RankSrc, int CurrentArg, class ... SubViewArgs>
struct SubviewLegalArgsCompileTime<Kokkos::LayoutContiguous<LayoutDest>,LayoutSrc,RankDest,RankSrc,CurrentArg,SubViewArgs...> {
  enum { value = SubviewLegalArgsCompileTime<LayoutDest,LayoutSrc,RankDest,RankSrc,CurrentArg,SubViewArgs...>::value };
};

template<class LayoutDest, class LayoutSrc, int RankDest, int RankSrc, int CurrentArg, class ... SubViewArgs>
struct SubviewLegalArgsCompileTime<LayoutDest,Kokkos::LayoutContiguous<LayoutSrc>,RankDest,RankSrc,CurrentArg,SubViewArgs...> {
  enum { value = SubviewLegalArgsCompileTime<LayoutDest,LayoutSrc,RankDest,RankSrc,CurrentArg,SubViewArgs...>::value };
};

template<class LayoutDest, class LayoutSrc, int RankDest, int RankSrc, int CurrentArg, class ... SubViewArgs>
struct SubviewLegalArgsCompileTime<Kokkos::LayoutContiguous<LayoutDest>,Kokkos::LayoutContiguous<LayoutSrc>,RankDest,RankSrc,CurrentArg,SubViewArgs...> {
  enum { value = SubviewLegalArgsCompileTime<LayoutDest,LayoutSrc,RankDest,RankSrc,CurrentArg,SubViewArgs...>::value };
};

// Subview mapping

template< class SrcTraits , class Arg0 , class ... Args >
struct ViewMapping
  < typename std::enable_if<(
      // Source view has FAD only
      std::is_same< typename SrcTraits::specialize
                  , ViewSpecializeSacadoFadContiguous >::value
      &&
      (
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutLeft >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutRight >::value ||
        std::is_same< typename SrcTraits::array_layout
                    , Kokkos::LayoutStride >::value
      )
      && !Sacado::Fad::is_fad_partition<Arg0>::value
    )
    >::type
  , SrcTraits
  , Arg0, Args ... >
{
private:

  static_assert( SrcTraits::rank == sizeof...(Args)+1 , "" );

  enum
    { RZ = false
    , R0 = bool(is_integral_extent<0,Arg0,Args...>::value)
    , R1 = bool(is_integral_extent<1,Arg0,Args...>::value)
    , R2 = bool(is_integral_extent<2,Arg0,Args...>::value)
    , R3 = bool(is_integral_extent<3,Arg0,Args...>::value)
    , R4 = bool(is_integral_extent<4,Arg0,Args...>::value)
    , R5 = bool(is_integral_extent<5,Arg0,Args...>::value)
    , R6 = bool(is_integral_extent<6,Arg0,Args...>::value)
    };

  // Public rank
  enum { rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3)
              + unsigned(R4) + unsigned(R5) + unsigned(R6) };

  // Whether right-most non-FAD rank is a range.
  enum { R0_rev = ( 0 == SrcTraits::rank ? RZ : (
                    1 == SrcTraits::rank ? R0 : (
                    2 == SrcTraits::rank ? R1 : (
                    3 == SrcTraits::rank ? R2 : (
                    4 == SrcTraits::rank ? R3 : (
                    5 == SrcTraits::rank ? R4 : (
                    6 == SrcTraits::rank ? R5 : R6 ))))))) };

  // Subview's layout
  typedef typename std::conditional<
      ( /* Same array layout IF */
        ( rank == 0 ) /* output rank zero */
        ||
        // OutputRank 1, InputLayout Left, Interval 0
        // because single stride one
        ( rank <= 1 && R0 && std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutLeft >::value )
        ||
        // OutputRank 1, InputLayout Right, Interval [InputRank-1]
        // because single stride one
        ( rank <= 1 && R0_rev && std::is_same< typename SrcTraits::array_layout , Kokkos::LayoutRight >::value )
        ), typename SrcTraits::array_layout , Kokkos::LayoutContiguous<Kokkos::LayoutStride,SrcTraits::array_layout::scalar_stride>
      >::type array_layout ;

  typedef typename SrcTraits::value_type  fad_type ;

  typedef typename std::conditional< rank == 0 , fad_type ,
          typename std::conditional< rank == 1 , fad_type * ,
          typename std::conditional< rank == 2 , fad_type ** ,
          typename std::conditional< rank == 3 , fad_type *** ,
          typename std::conditional< rank == 4 , fad_type **** ,
          typename std::conditional< rank == 5 , fad_type ***** ,
          typename std::conditional< rank == 6 , fad_type ****** ,
                                                 fad_type *******
          >::type >::type >::type >::type >::type >::type >::type
    data_type ;

public:

  typedef Kokkos::ViewTraits
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > traits_type ;

  typedef Kokkos::View
    < data_type
    , array_layout
    , typename SrcTraits::device_type
    , typename SrcTraits::memory_traits > type ;


  KOKKOS_INLINE_FUNCTION
  static void assign( ViewMapping< traits_type , typename traits_type::specialize > & dst
                    , ViewMapping< SrcTraits , typename SrcTraits::specialize > const & src
                    , Arg0 arg0 , Args ... args )
    {
      typedef ViewMapping< traits_type , typename traits_type::specialize > DstType ;
      typedef typename DstType::offset_type  dst_offset_type ;
      typedef typename DstType::array_offset_type  dst_array_offset_type ;
      typedef typename DstType::handle_type  dst_handle_type ;

      size_t offset;
      if (std::is_same< typename SrcTraits::array_layout, LayoutLeft >::value) {
        const SubviewExtents< SrcTraits::rank + 1 , rank + 1 >
          array_extents( src.m_array_offset.m_dim ,  Kokkos::ALL() , arg0 , args... );
        offset = src.m_array_offset( array_extents.domain_offset(0)
                                   , array_extents.domain_offset(1)
                                   , array_extents.domain_offset(2)
                                   , array_extents.domain_offset(3)
                                   , array_extents.domain_offset(4)
                                   , array_extents.domain_offset(5)
                                   , array_extents.domain_offset(6)
                                   , array_extents.domain_offset(7) );
        dst_array_offset_type dst_array_offset( src.m_array_offset ,
                                                array_extents );
        // For LayoutStride, we always use LayoutRight indexing (because we
        // don't know whether the original array was Left or Right), so we
        // need to swap the Fad dimension to the last and shift all of the
        // other dimensions left by 1
        if constexpr(std::is_same<typename traits_type::array_layout, LayoutStride>::value)
        {
          Kokkos::LayoutStride ls(
            dst_array_offset.m_dim.N0, dst_array_offset.m_stride.S0,
            dst_array_offset.m_dim.N1, dst_array_offset.m_stride.S1,
            dst_array_offset.m_dim.N2, dst_array_offset.m_stride.S2,
            dst_array_offset.m_dim.N3, dst_array_offset.m_stride.S3,
            dst_array_offset.m_dim.N4, dst_array_offset.m_stride.S4,
            dst_array_offset.m_dim.N5, dst_array_offset.m_stride.S5,
            dst_array_offset.m_dim.N6, dst_array_offset.m_stride.S6,
            dst_array_offset.m_dim.N7, dst_array_offset.m_stride.S7);
          auto t1 = ls.dimension[0];
          for (unsigned i=0; i<rank; ++i)
            ls.dimension[i] = ls.dimension[i+1];
          ls.dimension[rank] = t1;
          auto t2 = ls.stride[0];
          for (unsigned i=0; i<rank; ++i)
            ls.stride[i] = ls.stride[i+1];
          ls.stride[rank] = t2;
          dst.m_array_offset = dst_array_offset_type(std::integral_constant<unsigned, 0>(), ls);
        }
        else
          dst.m_array_offset = dst_array_offset;
      }
      else {
        const SubviewExtents< SrcTraits::rank + 1 , rank + 1 >
          array_extents( src.m_array_offset.m_dim , arg0 , args... , Kokkos::ALL() );
        offset = src.m_array_offset( array_extents.domain_offset(0)
                                   , array_extents.domain_offset(1)
                                   , array_extents.domain_offset(2)
                                   , array_extents.domain_offset(3)
                                   , array_extents.domain_offset(4)
                                   , array_extents.domain_offset(5)
                                   , array_extents.domain_offset(6)
                                   , array_extents.domain_offset(7) );
        dst.m_array_offset = dst_array_offset_type( src.m_array_offset ,
                                                    array_extents );
      }

      const SubviewExtents< SrcTraits::rank , rank >
        extents( src.m_impl_offset.m_dim , arg0 , args... );

      dst.m_impl_offset = dst_offset_type( src.m_impl_offset , extents );
      dst.m_impl_handle = dst_handle_type( src.m_impl_handle + offset );
      dst.m_fad_size = src.m_fad_size;
      dst.m_original_fad_size = src.m_original_fad_size;
      dst.m_fad_stride = src.m_fad_stride;
      dst.m_fad_index = src.m_fad_index;
    }

};

} // namespace Impl
} // namespace Kokkos

//---------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Partition mapping

template< class DataType, class ...P, unsigned Stride >
class ViewMapping<
  void,
  ViewTraits<DataType,P...> ,
  Sacado::Fad::Partition<Stride> 
  >
{
public:

  enum { is_assignable = true };
  enum { is_assignable_data_type = true };

  typedef ViewTraits<DataType,P...> src_traits;
  typedef ViewMapping< src_traits , typename src_traits::specialize >  src_type ;

  typedef typename src_type::offset_type::dimension_type src_dimension;
  typedef typename src_traits::value_type fad_type;
  typedef typename Sacado::LocalScalarType<fad_type,Stride>::type strided_fad_type;
  typedef typename
    ViewDataType< strided_fad_type , src_dimension >::type strided_data_type;
  typedef ViewTraits<strided_data_type,P...> dst_traits;
  typedef View<strided_data_type,P...> type;
  typedef ViewMapping< dst_traits , typename dst_traits::specialize >  dst_type ;

  KOKKOS_INLINE_FUNCTION static
  void assign( dst_type & dst
             , const src_type & src
             , const Sacado::Fad::Partition<Stride> & part )
    {
      if ( Stride != part.stride && Stride != 0 ) {
        Kokkos::abort("\n\n ******  Kokkos::View< Sacado::Fad ... > Invalid size in partitioned view assignment ******\n\n");
      }
      if ( src.m_fad_stride != 1 ) {
        Kokkos::abort("\n\n ******  Kokkos::View< Sacado::Fad ... > Can't partition already partitioned view ******\n\n");
      }

      dst.m_impl_handle = src.m_impl_handle ;
      dst.m_impl_offset  = src.m_impl_offset ;
      dst.m_array_offset  = src.m_array_offset ;

      // Assuming the alignment was choosen correctly for the partitioning,
      // each partition should get the same size.  This allows the use of SFad.
      dst.m_fad_size =
        (src.m_fad_size.value + part.stride-part.offset-1) / part.stride ;

      dst.m_original_fad_size = src.m_original_fad_size ;
      dst.m_fad_stride = part.stride ;
      dst.m_fad_index = part.offset ;
    }
};

} // namespace Impl
} // namespace Kokkos

#endif // defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#endif // defined(HAVE_SACADO_KOKKOS)

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_SACADO_FAD_HPP */
