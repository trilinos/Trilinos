// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_TEMPVIEWUTILS_HPP
#define TPETRA_DETAILS_TEMPVIEWUTILS_HPP

#include "Kokkos_Core.hpp"
#include "Tpetra_Details_isInterComm.hpp"

namespace Tpetra
{
namespace Details
{
namespace TempView
{

template<typename MemorySpace>
struct AlwaysMPISafe
{
  enum : bool {value = false};
};

template<>
struct AlwaysMPISafe<Kokkos::HostSpace>
{
  enum : bool {value = true};
};

#ifdef KOKKOS_ENABLE_CUDA
template<>
struct AlwaysMPISafe<Kokkos::CudaHostPinnedSpace>
{
  enum : bool {value = true};
};
#endif

#ifdef KOKKOS_ENABLE_HIP
template<>
struct AlwaysMPISafe<Kokkos::HIPHostPinnedSpace>
{
  enum : bool {value = true};
};
#endif

/// Get the contiguous layout that matches as many of the given views as possible. If neither or both arguments are contiguous, favor LayoutLeft.
template<typename View1, typename View2>
struct UnifiedContiguousLayout
{
  using L1 = typename View1::array_layout;
  using L2 = typename View2::array_layout;
  enum : bool {EitherLeft = std::is_same<L1, Kokkos::LayoutLeft>::value || std::is_same<L2, Kokkos::LayoutLeft>::value};
enum : bool {BothStride = std::is_same<L1, Kokkos::LayoutStride>::value && std::is_same<L2, Kokkos::LayoutStride>::value};
  using type = typename std::conditional<EitherLeft || BothStride, Kokkos::LayoutLeft, Kokkos::LayoutRight>::type;
};

/// Get a deep copy of SrcView with the requested contiguous layout and default memory traits. If SrcView already has that layout, just return it.
template<typename SrcView, typename Layout, typename std::enable_if<!std::is_same<typename SrcView::array_layout, Layout>::value>::type* = nullptr>
Kokkos::View<typename SrcView::data_type, Layout, typename SrcView::device_type>
toLayout(const SrcView& src)
{
  static_assert(!std::is_same<Kokkos::LayoutStride, Layout>::value,
      "TempView::toLayout: Layout must be contiguous (not LayoutStride)");
  Layout layout(src.extent(0), src.extent(1));
  Kokkos::View<typename SrcView::non_const_data_type, Layout, typename SrcView::device_type> dst(Kokkos::ViewAllocateWithoutInitializing(src.label()), layout);
  Kokkos::deep_copy(dst, src);
  return dst;
}

template<typename SrcView, typename Layout, typename std::enable_if<std::is_same<typename SrcView::array_layout, Layout>::value>::type* = nullptr>
Kokkos::View<typename SrcView::data_type, Layout, typename SrcView::device_type>
toLayout(const SrcView& src)
{
  if(src.span_is_contiguous())
  {
    return src;
  }
  else
  {
    //Even though the layout is already correct, it's not contiguous.
    Layout layout(src.extent(0), src.extent(1));
    Kokkos::View<typename SrcView::non_const_data_type, Layout, typename SrcView::device_type>
      result(Kokkos::ViewAllocateWithoutInitializing(src.label()), layout);
    Kokkos::deep_copy(result, src);
    return result;
  }
}

/// Get a copy of SrcView that is safe to use with MPI.
/// If already safe, just returns src.
/// If a new copy must be made, it will always be in HostSpace.
template<typename SrcView, bool AssumeGPUAware, typename = typename std::enable_if<AssumeGPUAware || AlwaysMPISafe<typename SrcView::memory_space>::value>::type>
SrcView
toMPISafe(const SrcView& src)
{
  using SrcLayout = typename SrcView::array_layout;
  static_assert(!std::is_same<SrcLayout, Kokkos::LayoutStride>::value, "toMPISafe requires that SrcView is contiguous");
  return toLayout<SrcView, SrcLayout>(src);
}

template<typename SrcView, bool AssumeGPUAware, typename = typename std::enable_if<!(AssumeGPUAware || AlwaysMPISafe<typename SrcView::memory_space>::value)>::type>
decltype(Kokkos::create_mirror_view_and_copy(std::declval<Kokkos::HostSpace>(), std::declval<SrcView>()))
toMPISafe(const SrcView& src)
{
  using SrcLayout = typename SrcView::array_layout;
  static_assert(!std::is_same<SrcLayout, Kokkos::LayoutStride>::value, "toMPISafe requires that SrcView is contiguous");
  auto srcContig = toLayout<SrcView, SrcLayout>(src);
  return Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), srcContig);
}

}}} //namespace Tpetra::Details::TempView

#endif

