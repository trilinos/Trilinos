// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_KOKKOS_VIEW_SUPPORT_INCLUDES
#error "This file can only be included by Sacado_Fad_Kokkos_View_Support.hpp"
#endif

// =====================================================================
// This file includes helpers to deal with local temporaries of Sacado
// Fad types.
// -- partition_scalar
// -- LocalScalarType
// -- ThreadLocalScalarType

namespace Sacado {

template <unsigned Stride, typename T>
KOKKOS_INLINE_FUNCTION const T &partition_scalar(const T &x) {
  return x;
}

// Type of local scalar type when partitioning a view
template <typename T, unsigned Stride> struct LocalScalarType {
  typedef T type;
};
template <typename T, unsigned Stride> struct LocalScalarType<const T, Stride> {
  typedef typename LocalScalarType<T, Stride>::type lst;
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
} // namespace Exp
} // namespace Fad
template <typename T, int N, unsigned Stride>
struct LocalScalarType<Fad::Exp::GeneralFad<Fad::Exp::StaticStorage<T, N>>,
                       Stride> {
  static const int Ns = (N + Stride - 1) / Stride;
  typedef Fad::Exp::GeneralFad<Fad::Exp::StaticStorage<T, Ns>> type;
};
#ifndef SACADO_NEW_FAD_DESIGN_IS_DEFAULT
namespace Fad {
template <typename T, int N> class SLFad;
}
template <typename T, int N, unsigned Stride>
struct LocalScalarType<Fad::SLFad<T, N>, Stride> {
  static const int Ns = (N + Stride - 1) / Stride;
  typedef Fad::SLFad<T, Ns> type;
};
#endif

// Type of local scalar type when partitioning a view
//
// For SFad, divde the array size by the given stride.  If it divides evenly,
// use SFad, otherwise use SLFad
namespace Fad {
namespace Exp {
template <typename T, typename U> class DynamicStorage;
template <typename T, int N> class StaticFixedStorage;
template <typename T, int N> class StaticStorage;
template <typename S> class GeneralFad;
} // namespace Exp
} // namespace Fad
template <typename T, int N, unsigned Stride>
struct LocalScalarType<Fad::Exp::GeneralFad<Fad::Exp::StaticFixedStorage<T, N>>,
                       Stride> {
  static const int Ns = (N + Stride - 1) / Stride;
  typedef typename std::conditional<
      Ns == N / Stride,
      Fad::Exp::GeneralFad<Fad::Exp::StaticFixedStorage<T, Ns>>,
      Fad::Exp::GeneralFad<Fad::Exp::StaticStorage<T, Ns>>>::type type;
};

#ifndef SACADO_NEW_FAD_DESIGN_IS_DEFAULT
namespace Fad {
template <typename T, int N> class SFad;
}
template <typename T, int N, unsigned Stride>
struct LocalScalarType<Fad::SFad<T, N>, Stride> {
  static const int Ns = (N + Stride - 1) / Stride;
  typedef typename std::conditional<Ns == N / Stride, Fad::SFad<T, Ns>,
                                    Fad::SLFad<T, Ns>>::type type;
};
#endif

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)

#ifndef SACADO_VIEW_CUDA_HIERARCHICAL_DFAD
template <unsigned Stride, typename T, typename U>
KOKKOS_INLINE_FUNCTION typename LocalScalarType<
    Fad::Exp::GeneralFad<Fad::Exp::DynamicStorage<T, U>>, Stride>::type
partition_scalar(
    const Fad::Exp::GeneralFad<Fad::Exp::DynamicStorage<T, U>> &x) {
  typedef typename LocalScalarType<
      Fad::Exp::GeneralFad<Fad::Exp::DynamicStorage<T, U>>, Stride>::type
      ret_type;
  const int size = (x.size() + blockDim.x - threadIdx.x - 1) / blockDim.x;
  const int offset = threadIdx.x;
  ret_type xp(size, x.val());

  // Note:  we can't use x.dx(offset+i*Stride) if
  // SACADO_VIEW_CUDA_HIERARCHICAL_DFAD_STRIDED is defined because it already
  // uses blockDim.x in its index calculation.  This approach should work
  // regardless
  const T *dx = x.dx();
  for (int i = 0; i < size; ++i)
    xp.fastAccessDx(i) = dx[offset + i * Stride];

  return xp;
}
#endif
template <unsigned Stride, typename T, int N>
KOKKOS_INLINE_FUNCTION typename LocalScalarType<
    Fad::Exp::GeneralFad<Fad::Exp::StaticStorage<T, N>>, Stride>::type
partition_scalar(const Fad::Exp::GeneralFad<Fad::Exp::StaticStorage<T, N>> &x) {
  typedef typename LocalScalarType<
      Fad::Exp::GeneralFad<Fad::Exp::StaticStorage<T, N>>, Stride>::type
      ret_type;
  const int size = (x.size() + blockDim.x - threadIdx.x - 1) / blockDim.x;
  const int offset = threadIdx.x;
  ret_type xp(size, x.val());
  for (int i = 0; i < size; ++i)
    xp.fastAccessDx(i) = x.fastAccessDx(offset + i * Stride);
  return xp;
}
template <unsigned Stride, typename T, int N>
KOKKOS_INLINE_FUNCTION typename LocalScalarType<
    Fad::Exp::GeneralFad<Fad::Exp::StaticFixedStorage<T, N>>, Stride>::type
partition_scalar(
    const Fad::Exp::GeneralFad<Fad::Exp::StaticFixedStorage<T, N>> &x) {
  typedef typename LocalScalarType<
      Fad::Exp::GeneralFad<Fad::Exp::StaticFixedStorage<T, N>>, Stride>::type
      ret_type;
  const int size = (x.size() + blockDim.x - threadIdx.x - 1) / blockDim.x;
  const int offset = threadIdx.x;
  ret_type xp(size, x.val());
  for (int i = 0; i < size; ++i)
    xp.fastAccessDx(i) = x.fastAccessDx(offset + i * Stride);
  return xp;
}
#endif

template <typename ViewType, typename Enabled = void>
struct ThreadLocalScalarType {
  typedef typename ViewType::non_const_value_type type;
};

template <typename ViewType>
struct ThreadLocalScalarType<
    ViewType,
    typename std::enable_if<is_view_fad_contiguous<ViewType>::value>::type> {
  typedef typename ViewType::traits TraitsType;
  // typedef Impl::ViewMapping<TraitsType, typename TraitsType::specialize>
  // MappingType; typedef typename MappingType::thread_local_scalar_type type;

  using fad_type = typename ViewType::value_type;
  enum { FadStaticDimension = Sacado::StaticSize<fad_type>::value };
  enum { PartitionedFadStride = TraitsType::array_layout::scalar_stride };

  // The partitioned static size -- this will be 0 if ParitionedFadStride
  // does not evenly divide FadStaticDimension
  enum {
    PartitionedFadStaticDimension =
        Impl::computeFadPartitionSize(FadStaticDimension, PartitionedFadStride)
  };
#ifdef KOKKOS_ENABLE_CUDA
  typedef typename Sacado::LocalScalarType<
      fad_type, unsigned(PartitionedFadStride)>::type strided_scalar_type;
  typedef typename std::conditional<
      std::is_same<typename TraitsType::execution_space, Kokkos::Cuda>::value,
      strided_scalar_type, fad_type>::type thread_local_scalar_type;
#elif defined(KOKKOS_ENABLE_HIP)
  typedef typename Sacado::LocalScalarType<
      fad_type, unsigned(PartitionedFadStride)>::type strided_scalar_type;
  typedef typename std::conditional<
      std::is_same<typename TraitsType::execution_space, Kokkos::HIP>::value,
      strided_scalar_type, fad_type>::type thread_local_scalar_type;
#else
  typedef fad_type thread_local_scalar_type;
#endif
  typedef thread_local_scalar_type type;
};
} // namespace Sacado
