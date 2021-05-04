// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
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
struct AlwaysMPISafe<Kokkos::Experimental::HIPHostPinnedSpace>
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
  Layout layout(src.extent(0), src.extent(1), src.extent(2), src.extent(3), src.extent(4), src.extent(5), src.extent(6), src.extent(7)); 
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
    Layout layout(src.extent(0), src.extent(1), src.extent(2), src.extent(3), src.extent(4), src.extent(5), src.extent(6), src.extent(7)); 
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

