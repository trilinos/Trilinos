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

//
// How wide a vector must a partitioned Fad be sized for?
//
// Not the width that was requested.  Kokkos does not guarantee a TeamPolicy
// gets the vector width it asked for on SYCL -- it narrows the request in two
// places, and the second cannot be seen from the host:
//
//   Kokkos_SYCL_TeamPolicy.hpp:100
//       determine_vector_length() clamps to the device's largest sub-group
//       size and rounds down to a power of two.
//   Kokkos_SYCL_ParallelFor_Team.hpp:86
//       final_vector_size = min(requested, max_sub_group_size), where that
//       maximum belongs to the BUILT KERNEL, not the device, and falls as the
//       kernel's register pressure rises.  No host-side API reports it.
//
// The loop bounds in FadAccessor::access() and partition_scalar() use the
// width actually delivered, but the local Fad's capacity is fixed at compile
// time from the layout stride.  When Kokkos narrows the width, each thread
// becomes responsible for MORE derivative components than its local Fad can
// hold.  For a statically sized Fad that is a write past the end of a
// fixed-size member array, in device code, with no diagnostic.
//
// Register pressure makes the exposure worse than it first appears:  the
// kernels most likely to be narrowed are the ones carrying the most derivative
// components, which are exactly the kernels the hierarchical scheme exists to
// make fast.
//
// So on SYCL we size for the NARROWEST width Kokkos could deliver rather than
// the one requested, and keep the component count a runtime value within that
// capacity.  That is correct whatever width arrives.  It costs registers:  a
// stride of 32 sized against a floor of 8 reserves 4x the derivative storage
// per thread.  This is a deliberate trade of speed for correctness, made
// because the failure it prevents is silent.
//
// SACADO_SYCL_MIN_VECTOR_LENGTH is that floor.  The default of 8 is the
// smallest SIMD width Intel GPUs compile to, so it is safe everywhere.  A
// build targeting a device whose sub-group sizes start higher can raise it and
// recover the registers:  SyclIndexProbe reports sub_group_sizes for the
// device, and on Ponte Vecchio it is [16, 32], where 16 is safe and halves the
// cost.
//
#if defined(KOKKOS_ENABLE_SYCL)
// Configured, not defaulted here:  this value selects a type, so every
// translation unit in a build has to agree on it.  Sacado_SYCL_MIN_VECTOR_LENGTH
// is a CMake cache variable and reaches us through Sacado_config.h.
#ifndef SACADO_SYCL_MIN_VECTOR_LENGTH
#error "SACADO_SYCL_MIN_VECTOR_LENGTH is missing from the Sacado_config.h this file is seeing.  Usually that header is stale:  re-run CMake so the build tree regenerates it, and check that no previously installed Sacado_config.h comes first on the include path.  Do not define the macro in a source file -- a per-file value would give different translation units different View types."
#endif
// Stride itself when it is already at or below the floor:  Kokkos only ever
// narrows a request, never widens it.
#define SACADO_IMPL_SIZING_STRIDE(Stride)                                      \
  ((Stride) < SACADO_SYCL_MIN_VECTOR_LENGTH ? (Stride)                         \
                                            : SACADO_SYCL_MIN_VECTOR_LENGTH)
#endif

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
template <typename T, int N> class StaticStorage;
template <typename S> class GeneralFad;
} // namespace Fad
template <typename T, int N, unsigned Stride>
struct LocalScalarType<Fad::GeneralFad<Fad::StaticStorage<T, N>>,
                       Stride> {
#if defined(KOKKOS_ENABLE_SYCL)
  static const unsigned SizingStride = SACADO_IMPL_SIZING_STRIDE(Stride);
#else
  static const unsigned SizingStride = Stride;
#endif
  static const int Ns = (N + SizingStride - 1) / SizingStride;
  typedef Fad::GeneralFad<Fad::StaticStorage<T, Ns>> type;
};
// Type of local scalar type when partitioning a view
//
// For SFad, divde the array size by the given stride.  If it divides evenly,
// use SFad, otherwise use SLFad
namespace Fad {
template <typename T, typename U> class DynamicStorage;
template <typename T, int N> class StaticFixedStorage;
template <typename T, int N> class StaticStorage;
template <typename S> class GeneralFad;
} // namespace Fad
template <typename T, int N, unsigned Stride>
struct LocalScalarType<Fad::GeneralFad<Fad::StaticFixedStorage<T, N>>,
                       Stride> {
#if defined(KOKKOS_ENABLE_SYCL)
  static const unsigned SizingStride = SACADO_IMPL_SIZING_STRIDE(Stride);
  static const int Ns = (N + SizingStride - 1) / SizingStride;
  // Always the runtime-sized storage, never StaticFixedStorage.  How many
  // components a thread actually holds depends on the width Kokkos delivered,
  // so a type whose size is fixed at compile time is wrong for every width but
  // one.  Capacity is compile time, count is runtime.
  typedef Fad::GeneralFad<Fad::StaticStorage<T, Ns>> type;
#else
  static const int Ns = (N + Stride - 1) / Stride;
  typedef typename std::conditional<
      Ns == N / Stride,
      Fad::GeneralFad<Fad::StaticFixedStorage<T, Ns>>,
      Fad::GeneralFad<Fad::StaticStorage<T, Ns>>>::type type;
#endif
};

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) ||               \
    defined(__SYCL_DEVICE_ONLY__)

#ifndef SACADO_VIEW_CUDA_HIERARCHICAL_DFAD
template <unsigned Stride, typename T, typename U>
KOKKOS_INLINE_FUNCTION typename LocalScalarType<
    Fad::GeneralFad<Fad::DynamicStorage<T, U>>, Stride>::type
partition_scalar(
    const Fad::GeneralFad<Fad::DynamicStorage<T, U>> &x) {
  typedef typename LocalScalarType<
      Fad::GeneralFad<Fad::DynamicStorage<T, U>>, Stride>::type
      ret_type;
  if (Stride == 1u)
    return x;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  const int size = (x.size() + blockDim.x - threadIdx.x - 1) / blockDim.x;
  const int offset = threadIdx.x;
#else
  // See the note in Sacado_Fad_Kokkos_View_Support.hpp:  this is the SYCL
  // analogue of the threadIdx.x / blockDim.x pair above, and it reads as lane 0
  // of a width-1 vector in a flat kernel.
  const auto item = sycl::ext::oneapi::this_work_item::get_nd_item<2>();
  const int lane = item.get_local_id(1);
  const int vec = item.get_local_range(1);
  const int size = (x.size() + vec - lane - 1) / vec;
  const int offset = lane;
#endif
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
    Fad::GeneralFad<Fad::StaticStorage<T, N>>, Stride>::type
partition_scalar(const Fad::GeneralFad<Fad::StaticStorage<T, N>> &x) {
  typedef typename LocalScalarType<
      Fad::GeneralFad<Fad::StaticStorage<T, N>>, Stride>::type
      ret_type;
  if (Stride == 1u)
    return x;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  const int size = (x.size() + blockDim.x - threadIdx.x - 1) / blockDim.x;
  const int offset = threadIdx.x;
#else
  // See the note in Sacado_Fad_Kokkos_View_Support.hpp:  this is the SYCL
  // analogue of the threadIdx.x / blockDim.x pair above, and it reads as lane 0
  // of a width-1 vector in a flat kernel.
  const auto item = sycl::ext::oneapi::this_work_item::get_nd_item<2>();
  const int lane = item.get_local_id(1);
  const int vec = item.get_local_range(1);
  const int size = (x.size() + vec - lane - 1) / vec;
  const int offset = lane;
#endif
  ret_type xp(size, x.val());
  for (int i = 0; i < size; ++i)
    xp.fastAccessDx(i) = x.fastAccessDx(offset + i * Stride);
  return xp;
}
template <unsigned Stride, typename T, int N>
KOKKOS_INLINE_FUNCTION typename LocalScalarType<
    Fad::GeneralFad<Fad::StaticFixedStorage<T, N>>, Stride>::type
partition_scalar(
    const Fad::GeneralFad<Fad::StaticFixedStorage<T, N>> &x) {
  typedef typename LocalScalarType<
      Fad::GeneralFad<Fad::StaticFixedStorage<T, N>>, Stride>::type
      ret_type;
  if (Stride == 1u)
    return x;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  const int size = (x.size() + blockDim.x - threadIdx.x - 1) / blockDim.x;
  const int offset = threadIdx.x;
#else
  // See the note in Sacado_Fad_Kokkos_View_Support.hpp:  this is the SYCL
  // analogue of the threadIdx.x / blockDim.x pair above, and it reads as lane 0
  // of a width-1 vector in a flat kernel.
  const auto item = sycl::ext::oneapi::this_work_item::get_nd_item<2>();
  const int lane = item.get_local_id(1);
  const int vec = item.get_local_range(1);
  const int size = (x.size() + vec - lane - 1) / vec;
  const int offset = lane;
#endif
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
#elif defined(KOKKOS_ENABLE_SYCL)
  typedef typename Sacado::LocalScalarType<
      fad_type, unsigned(PartitionedFadStride)>::type strided_scalar_type;
  typedef typename std::conditional<
      std::is_same<typename TraitsType::execution_space, Kokkos::SYCL>::value,
      strided_scalar_type, fad_type>::type thread_local_scalar_type;
#else
  typedef fad_type thread_local_scalar_type;
#endif
  typedef thread_local_scalar_type type;
};
} // namespace Sacado
