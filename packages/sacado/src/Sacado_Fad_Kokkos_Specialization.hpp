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

// This file contains pieces of Kokkos that still need to be
// overloaded/specialized for LayoutContiguous!
// -- Kokkos::Impl::SubviewLegalArgsCompileTime
// -- Kokkos::subview
// -- Kokkos::subdynrankview
// -- Kokkos::resize (when using Hierarchical Fad)

// Deal with subview of LayoutContiguous which isn't yet a real mdspan Layout
namespace Kokkos {
namespace Impl {

// Rules for subview arguments and layouts matching

template <class LayoutDest, unsigned StrideDst, class LayoutSrc, int RankDest,
          int RankSrc, int CurrentArg, class... SubViewArgs>
struct SubviewLegalArgsCompileTime<
    Kokkos::LayoutContiguous<LayoutDest, StrideDst>, LayoutSrc, RankDest,
    RankSrc, CurrentArg, SubViewArgs...> {
  enum {
    value =
        SubviewLegalArgsCompileTime<LayoutDest, LayoutSrc, RankDest, RankSrc,
                                    CurrentArg, SubViewArgs...>::value
  };
};

template <class LayoutDest, class LayoutSrc, unsigned StrideSrc, int RankDest,
          int RankSrc, int CurrentArg, class... SubViewArgs>
struct SubviewLegalArgsCompileTime<
    LayoutDest, Kokkos::LayoutContiguous<LayoutSrc, StrideSrc>, RankDest,
    RankSrc, CurrentArg, SubViewArgs...> {
  enum {
    value =
        SubviewLegalArgsCompileTime<LayoutDest, LayoutSrc, RankDest, RankSrc,
                                    CurrentArg, SubViewArgs...>::value
  };
};

template <class LayoutDest, unsigned StrideDest, class LayoutSrc,
          unsigned StrideSrc, int RankDest, int RankSrc, int CurrentArg,
          class... SubViewArgs>
struct SubviewLegalArgsCompileTime<
    Kokkos::LayoutContiguous<LayoutDest, StrideDest>,
    Kokkos::LayoutContiguous<LayoutSrc, StrideSrc>, RankDest, RankSrc,
    CurrentArg, SubViewArgs...> {
  enum {
    value =
        SubviewLegalArgsCompileTime<LayoutDest, LayoutSrc, RankDest, RankSrc,
                                    CurrentArg, SubViewArgs...>::value
  };
};

template <class DstT, class DstL, unsigned DstS, class... DstArgs, class SrcT,
          class SrcL, unsigned SrcS, class... SrcArgs, class... Args>
struct CommonSubview<
    Kokkos::View<DstT, Kokkos::LayoutContiguous<DstL, DstS>, DstArgs...>,
    Kokkos::View<SrcT, Kokkos::LayoutContiguous<SrcL, SrcS>, SrcArgs...>,
    Args...> {
  using DstType =
      Kokkos::View<DstT, Kokkos::LayoutContiguous<DstL, DstS>, DstArgs...>;
  using SrcType =
      Kokkos::View<SrcT, Kokkos::LayoutContiguous<SrcL, SrcS>, SrcArgs...>;
  using dst_subview_type =
      decltype(subview(std::declval<DstType>(), std::declval<Args>()...));
  using src_subview_type =
      decltype(subview(std::declval<SrcType>(), std::declval<Args>()...));
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType &dst, const SrcType &src, const Args &...args)
      : dst_sub(subview(dst, args...)), src_sub(subview(src, args...)) {}
};

} // namespace Impl
} // namespace Kokkos

namespace {
template <class T, size_t N> struct data_type_construct {
  using type = typename data_type_construct<T *, N - 1>::type;
};
template <class T> struct data_type_construct<T, 0> {
  using type = T;
};
} // namespace
namespace Kokkos {
// This is needed to deal with the return Layout Deduction for LayoutContiguous
// ...
template <class D, class LayoutSrc, unsigned StrideSrc, class... P,
          class... Args>
KOKKOS_INLINE_FUNCTION auto subview(
    const View<D, Kokkos::LayoutContiguous<LayoutSrc, StrideSrc>, P...> &src,
    Args... args) {
  using view_t = View<D, Kokkos::LayoutContiguous<LayoutSrc, StrideSrc>, P...>;
  auto submapping_result = submdspan_mapping(
      src.mapping(), Impl::transform_kokkos_slice_to_mdspan_slice(args)...);
  using sub_data_type = typename data_type_construct<
      typename view_t::value_type,
      decltype(submapping_result.mapping)::extents_type::rank()>::type;
  using layout_t = std::conditional_t<
      std::is_same_v<typename decltype(submapping_result.mapping)::layout_type,
                     layout_stride>,
      LayoutStride, LayoutSrc>;
  return View<sub_data_type, LayoutContiguous<layout_t, StrideSrc>,
              typename view_t::device_type, typename view_t::memory_traits>(
      src.accessor().offset(src.data_handle(), submapping_result.offset),
      submapping_result.mapping, src.accessor());
}

// This is needed to deal with the return Layout Deduction for LayoutContiguous
// ...
template <class D, class LayoutSrc, unsigned StrideSrc, class... P,
          class... Args>
KOKKOS_INLINE_FUNCTION auto
subview(const DynRankView<D, Kokkos::LayoutContiguous<LayoutSrc, StrideSrc>,
                          P...> &src,
        Args... args) {
  static_assert(View<D, P...>::rank == sizeof...(Args),
                "subview requires one argument for each source View rank");

  using sub_mdspan_t = decltype(submdspan(
      src.to_mdspan(), Impl::transform_kokkos_slice_to_mdspan_slice(args)...));
  if constexpr (std::is_same_v<typename sub_mdspan_t::layout_type,
                               layout_stride>) {
    return typename Kokkos::Impl::ViewMapping<
        void /* deduce subview type from source view traits */
        ,
        typename Impl::RemoveAlignedMemoryTrait<
            D, Kokkos::LayoutContiguous<LayoutStride, StrideSrc>, P...>::type,
        Args...>::type(src, args...);
  } else {
    return typename Kokkos::Impl::ViewMapping<
        void /* deduce subview type from source view traits */
        ,
        typename Impl::RemoveAlignedMemoryTrait<
            D, Kokkos::LayoutContiguous<LayoutSrc, StrideSrc>, P...>::type,
        Args...>::type(src, args...);
  }
}

template <class T, class LayoutSrc, unsigned StrideSrc, class... DRVArgs,
          class SubArg0 = int, class SubArg1 = int, class SubArg2 = int,
          class SubArg3 = int, class SubArg4 = int, class SubArg5 = int,
          class SubArg6 = int>
KOKKOS_INLINE_FUNCTION auto
subdynrankview(const DynRankView<T, LayoutContiguous<LayoutSrc, StrideSrc>,
                                 DRVArgs...> &drv,
               SubArg0 arg0 = SubArg0{}, SubArg1 arg1 = SubArg1{},
               SubArg2 arg2 = SubArg2{}, SubArg3 arg3 = SubArg3{},
               SubArg4 arg4 = SubArg4{}, SubArg5 arg5 = SubArg5{},
               SubArg6 arg6 = SubArg6{}) {
  auto sub = subview(drv.DownCast(), arg0, arg1, arg2, arg3, arg4, arg5, arg6);
  using sub_t = decltype(sub);
  size_t new_rank = (drv.rank() > 0 && !std::is_integral_v<SubArg0> ? 1 : 0) +
                    (drv.rank() > 1 && !std::is_integral_v<SubArg1> ? 1 : 0) +
                    (drv.rank() > 2 && !std::is_integral_v<SubArg2> ? 1 : 0) +
                    (drv.rank() > 3 && !std::is_integral_v<SubArg3> ? 1 : 0) +
                    (drv.rank() > 4 && !std::is_integral_v<SubArg4> ? 1 : 0) +
                    (drv.rank() > 5 && !std::is_integral_v<SubArg5> ? 1 : 0) +
                    (drv.rank() > 6 && !std::is_integral_v<SubArg6> ? 1 : 0);

  using return_type = DynRankView<
      typename sub_t::value_type,
      typename sub_t::array_layout, // LayoutContiguous<LayoutStride,
                                    // StrideSrc>,
      typename sub_t::device_type, typename sub_t::memory_traits>;

  auto layout = sub.layout().base_layout();
  for (int i = new_rank; i < 8; i++)
    layout.dimension[i] = 1;
  if constexpr (std::is_same_v<decltype(layout), LayoutStride>)
    for (int i = new_rank; i < 8; i++)
      layout.stride[i] = 1;

  return return_type{
      typename return_type::view_type(
          sub.data_handle(),
          Impl::mapping_from_array_layout<typename return_type::mapping_type>(
              layout),
          sub.accessor()),
      new_rank};
}

template <class T, class LayoutSrc, unsigned StrideSrc, class... DRVArgs,
          class SubArg0 = int, class SubArg1 = int, class SubArg2 = int,
          class SubArg3 = int, class SubArg4 = int, class SubArg5 = int,
          class SubArg6 = int>
KOKKOS_INLINE_FUNCTION auto
subview(const DynRankView<T, LayoutContiguous<LayoutSrc, StrideSrc>, DRVArgs...>
            &drv,
        SubArg0 arg0 = SubArg0{}, SubArg1 arg1 = SubArg1{},
        SubArg2 arg2 = SubArg2{}, SubArg3 arg3 = SubArg3{},
        SubArg4 arg4 = SubArg4{}, SubArg5 arg5 = SubArg5{},
        SubArg6 arg6 = SubArg6{}) {
  return subdynrankview(drv, arg0, arg1, arg2, arg3, arg4, arg5, arg6);
}
} // namespace Kokkos

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL)
namespace Kokkos {
template <class SrcT, class SrcL, unsigned SrcS, class... SrcArgs>
void resize(
    Kokkos::View<SrcT, Kokkos::LayoutContiguous<SrcL, SrcS>, SrcArgs...> &src,
    const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
    const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
    const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
    const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
    const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
    const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
    const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
    const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
  const size_t new_extents[8] = {n0, n1, n2, n3, n4, n5, n6, n7};
  bool size_mismatch = false;
  using view_t =
      Kokkos::View<SrcT, Kokkos::LayoutContiguous<SrcL, SrcS>, SrcArgs...>;
  for (int r = 0; r < (int)src.rank(); r++) {
    if (new_extents[r] != src.extent(r))
      size_mismatch = true;
  }
  if (size_mismatch) {
    using exec_space = typename view_t::execution_space;
    auto dst = view_t(src.label(), n0, n1, n2);
    size_t total_extent = 1;
    for (size_t r = 0; r < view_t::rank(); r++)
      total_extent *= Kokkos::min(src.extent(r), dst.extent(r));
    // It looks like SFAD only works with 64 wide vector in HIP
#ifdef KOKKOS_ENABLE_HIP
    size_t vector_size = 64;
#else
    size_t vector_size = 32;
#endif

    // Just arbitraryly using team_size = 1 for low concurrency backends (i.e.
    // CPUs)
    size_t team_size =
        exec_space().concurrency() > 1000 ? 512 / vector_size : 1;

    size_t num_teams = (total_extent + team_size - 1) / team_size;

    Kokkos::parallel_for(
        "Sacado::view_copy Hierarchical",
        Kokkos::TeamPolicy<exec_space>(num_teams, team_size, vector_size),
        KOKKOS_LAMBDA(
            const typename Kokkos::TeamPolicy<exec_space>::member_type &team) {
          size_t ii = team.league_rank() * team.team_size() + team.team_rank();
          if (ii >= total_extent)
            return;
          // work around capture restriction
          if (dst.data() == src.data())
            return;
          if constexpr (view_t::rank() == 0)
            dst() == src();
          else if constexpr (view_t::rank() == 1) {
            dst(ii) == src(ii);
          } else if constexpr (view_t::rank() == 2) {
            int i1 = ii % src.extent(1);
            int i0 = ii / src.extent(1);
            dst(i0, i1) = src(i0, i1);
          } else if constexpr (view_t::rank() == 3) {
            int i2 = ii % src.extent(2);
            int i1 = (ii / src.extent(2)) % src.extent(1);
            int i0 = ii / (src.extent(2) * src.extent(1));
            dst(i0, i1, i2) = src(i0, i1, i2);
          } else if constexpr (view_t::rank() == 4) {
            int i3 = ii % src.extent(3);
            int i2 = (ii / src.extent(3)) % src.extent(2);
            int i1 = (ii / (src.extent(3) * src.extent(2))) % src.extent(1);
            int i0 = (ii / (src.extent(3) * src.extent(2) * src.extent(1)));
            dst(i0, i1, i2, i3) = src(i0, i1, i2, i3);
          } else if constexpr (view_t::rank() == 5) {
            int i4 = ii % src.extent(4);
            int i3 = (ii / src.extent(4)) % src.extent(3);
            int i2 = (ii / (src.extent(4) * src.extent(3))) % src.extent(2);
            int i1 = (ii / (src.extent(4) * src.extent(3) * src.extent(2))) %
                     src.extent(1);
            int i0 = (ii / (src.extent(4) * src.extent(3) * src.extent(2) *
                            src.extent(1)));
            dst(i0, i1, i2, i3, i4) = src(i0, i1, i2, i3, i4);
          } else if constexpr (view_t::rank() == 5) {
            int i5 = ii % src.extent(5);
            int i4 = (ii / src.extent(5)) % src.extent(4);
            int i3 = (ii / (src.extent(5) * src.extent(4))) % src.extent(3);
            int i2 = (ii / (src.extent(5) * src.extent(4) * src.extent(3))) %
                     src.extent(2);
            int i1 = (ii / (src.extent(5) * src.extent(4) * src.extent(3) *
                            src.extent(2))) %
                     src.extent(1);
            int i0 = (ii / (src.extent(5) * src.extent(4) * src.extent(3) *
                            src.extent(2) * src.extent(1)));
            dst(i0, i1, i2, i3, i4, i5) = src(i0, i1, i2, i3, i4, i5);
          }
        });
    Kokkos::fence();
    src = dst;
  }
}
} // namespace Kokkos
#endif // SACADO_VIEW_CUDA_HIERARCHICAL
