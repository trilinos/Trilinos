// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_KOKKOS_VIEW_FACTORY_HPP
#define SACADO_FAD_KOKKOS_VIEW_FACTORY_HPP

#include <type_traits>

#include "Sacado_Traits.hpp"
#include "Sacado_Fad_Kokkos.hpp"

namespace Sacado {

namespace Impl {

// Class to determine the value_type for a view as a function of one or more
// input views
template <class ... ViewPack>
struct ViewFactoryType {};

template <class View>
struct ViewFactoryType<View> {
  typedef typename View::value_type type;
};

template <class View, class ... ViewPack>
struct ViewFactoryType<View,ViewPack...> {
  typedef typename Sacado::Promote<
    typename View::value_type,
    typename ViewFactoryType<ViewPack...>::type
    >::type type;
};

// Whether a view constructor property refers to existing memory, i.e., it is
// a pointer (that is not a label) or a Kokkos::Impl::ViewCtorProp holding a
// pointer, in which case the view factory creates an unmanaged view.  Unlike
// allocating views, this can be done on device as well.
template <class CtorProp>
struct IsUnmanagedViewCtorProp :
    std::bool_constant<std::is_pointer_v<CtorProp> &&
                       !std::is_convertible_v<CtorProp, const char*>> {};

template <class ... P>
struct IsUnmanagedViewCtorProp<Kokkos::Impl::ViewCtorProp<P...>> :
    std::bool_constant<Kokkos::Impl::ViewCtorProp<P...>::has_pointer> {};

// Pointer held by a constructor property of an unmanaged view
template <class T>
KOKKOS_INLINE_FUNCTION
T* view_ctor_prop_pointer(T* ptr) { return ptr; }

template <class ... P>
KOKKOS_INLINE_FUNCTION
auto view_ctor_prop_pointer(const Kokkos::Impl::ViewCtorProp<P...>& prop) {
  return Kokkos::Impl::get_property<Kokkos::Impl::PointerTag>(prop);
}

// Scalar (e.g., Fad) dimension given by the Kokkos::Impl::AccessorArg_t of a
// constructor property if it has one, otherwise default_value
template <class CtorProp>
KOKKOS_INLINE_FUNCTION
size_t view_ctor_prop_accessor_arg(const CtorProp&,
                                   const size_t default_value) {
  return default_value;
}

template <class ... P>
KOKKOS_INLINE_FUNCTION
size_t view_ctor_prop_accessor_arg(
  [[maybe_unused]] const Kokkos::Impl::ViewCtorProp<P...>& prop,
  const size_t default_value) {
  if constexpr (Kokkos::Impl::ViewCtorProp<P...>::has_accessor_arg)
    return Kokkos::Impl::get_property<Kokkos::Impl::AccessorArgTag>(prop).value;
  else
    return default_value;
}

// Enables a ViewFactory wrapper for the given input view and constructor
// property types, the latter creating an unmanaged view or not
template <class InputViewType, class CtorProp, bool Unmanaged, class T>
using EnableIfViewFactoryArgs = typename std::enable_if<
  (Kokkos::is_view<InputViewType>::value ||
   Kokkos::is_dyn_rank_view<InputViewType>::value) &&
  (IsUnmanagedViewCtorProp<CtorProp>::value == Unmanaged),
  T>::type;

}

// Function to compute the scalar dimension (e.g., Fad dimesion) from one or
// more views.  It relies on the overload for a single view provided by Sacado

// Traits class used to create a view for a given rank and dimension as a
// function of one or more views.  The value_type for the view is determined
// by value_type, and the view is created through the create_view() function.
// The calling code must determine the rank and dimensions of the view to create
// however internal Sacado dimension will be determined automatically.
template <class ... ViewPack>
struct ViewFactory {

  typedef typename Impl::ViewFactoryType<ViewPack...>::type value_type;

  template<class T>
  struct IsCtorProp : std::false_type{};
  template<class ... Args>
  struct IsCtorProp<Kokkos::Impl::ViewCtorProp<Args...>> : std::true_type{};

  // Create a view allocating its memory, prop being a label or a
  // Kokkos::Impl::ViewCtorProp without a pointer (host only)
  template <class ResultView, class CtorProp, class ... Dims>
  static
  typename std::enable_if<!Impl::IsUnmanagedViewCtorProp<CtorProp>::value,
                          ResultView>::type
  create_view(const ViewPack& ... views,
              const CtorProp& prop,
              const Dims ... dims) {

    auto layout = create_layout<ResultView>(views..., dims...);

    if constexpr (is_view_fad<ResultView>::value) {
      const size_t fad_size = fad_dimension(views...);
      // turns out CtorProp can be just a label...
      if constexpr (IsCtorProp<CtorProp>::value) {
        return ResultView(Kokkos::Impl::with_properties_if_unset(prop, Kokkos::Impl::AccessorArg_t{fad_size}), layout);
      } else {
        return ResultView(Kokkos::view_alloc(prop, Kokkos::Impl::AccessorArg_t{fad_size}), layout);
      }
    } else {
      return ResultView(prop, layout);
    }
  }

  // Create an unmanaged view of the memory given by prop, a pointer or a
  // Kokkos::Impl::ViewCtorProp holding one (host and device)
  template <class ResultView, class CtorProp, class ... Dims>
  KOKKOS_INLINE_FUNCTION
  static
  typename std::enable_if<Impl::IsUnmanagedViewCtorProp<CtorProp>::value,
                          ResultView>::type
  create_view(const ViewPack& ... views,
              const CtorProp& prop,
              const Dims ... dims) {

    auto layout = create_layout<ResultView>(views..., dims...);

    if constexpr (is_view_fad<ResultView>::value) {
      // Fad dimension given by prop, if any, takes precedence
      const size_t fad_size =
        Impl::view_ctor_prop_accessor_arg(prop, fad_dimension(views...));
      const auto ptr = Impl::view_ctor_prop_pointer(prop);
#if KOKKOS_VERSION_LESS(5,3,0)
      // Before Kokkos 5.3, the DynRankView constructor wrapping a pointer
      // computes a rank that is one too small when the Fad dimension is passed
      // through an AccessorArg_t (fixed by kokkos/kokkos#9651).  Pass the Fad
      // dimension as the trailing extent of the layout instead, which that
      // constructor handles correctly.
      constexpr bool fad_size_in_layout =
        Kokkos::is_dyn_rank_view<ResultView>::value;
#else
      constexpr bool fad_size_in_layout = false;
#endif
      if constexpr (fad_size_in_layout) {
        layout.dimension[sizeof...(Dims)] = fad_size;
        return ResultView(Kokkos::view_wrap(ptr), layout);
      } else {
        return ResultView(Kokkos::view_wrap(ptr, Kokkos::Impl::AccessorArg_t{fad_size}), layout);
      }
    } else {
      return ResultView(prop, layout);
    }
  }

private:

  // Fad dimension (derivative dimension + 1) of the view to create
  KOKKOS_INLINE_FUNCTION
  static size_t fad_dimension(const ViewPack& ... views) {
    const size_t fad_size = Sacado::dimension_scalar(views...);
    return fad_size == 0 ? 1 : fad_size;
  }

  // Layout of the view to create.  For views of Fad types, the Fad dimension
  // is not part of the layout, as it is passed through an AccessorArg_t.
  template <class ResultView, class ... Dims>
  KOKKOS_INLINE_FUNCTION
  static typename ResultView::array_layout
  create_layout(const ViewPack& ... views,
                const Dims ... dims) {

    using nc_value_type = typename ResultView::non_const_value_type;
    constexpr bool is_scalar = Sacado::IsScalarType<nc_value_type>::value;
    constexpr bool is_dyn_rank = Kokkos::is_dyn_rank_view<ResultView>::value;

    // rank == number of arguments
    constexpr unsigned rank = sizeof...(Dims);

    // Check rank is valid
    static_assert( rank <= 7, "Invalid rank...too many dimension arguments" );

    // Create layout from our dimension arguments
    typename ResultView::array_layout layout(dims...);

    // Set scalar dimension
    layout.dimension[rank] = Sacado::dimension_scalar(views...);

    // Handle the case where all of the input view's are scalar's, but the
    // result isn't (e.g., a Fad), in which case we have to specify a valid
    // scalar dimension
    if (!is_scalar && layout.dimension[rank] == 0)
      layout.dimension[rank] = 1;

    // Reconstruct layout for dynamic rank
    if (is_dyn_rank) {
      constexpr unsigned r = is_scalar ? rank : rank + 1;
      layout = Kokkos::Impl::reconstructLayout(layout, r);
    }

    if constexpr (is_view_fad<ResultView>::value)
      layout.dimension[rank] = is_dyn_rank ? KOKKOS_INVALID_INDEX : KOKKOS_IMPL_CTOR_DEFAULT_ARG;

    return layout;
  }

};

//! Wrapper to simplify use of Sacado ViewFactory
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename ... Dims>
Impl::EnableIfViewFactoryArgs<InputViewType, CtorProp, false, ResultViewType>
createDynRankViewWithType(const InputViewType& a,
                          const CtorProp& prop,
                          const Dims... dims)
{
  using view_factory = Sacado::ViewFactory<InputViewType>;
  return view_factory::template create_view<ResultViewType>(a,prop,dims...);
}

//! Wrapper to simplify use of Sacado ViewFactory for unmanaged views
//! (prop is a pointer or holds one), callable on device as well
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename ... Dims>
KOKKOS_INLINE_FUNCTION
Impl::EnableIfViewFactoryArgs<InputViewType, CtorProp, true, ResultViewType>
createDynRankViewWithType(const InputViewType& a,
                          const CtorProp& prop,
                          const Dims... dims)
{
  using view_factory = Sacado::ViewFactory<InputViewType>;
  return view_factory::template create_view<ResultViewType>(a,prop,dims...);
}

namespace Impl {
  // Helper type trait to determine type of resulting DynRankView from
  // createDynRankView below
  template <typename InputView>
  struct ResultDynRankView {
    // Allow for use of LayoutStride in InputViewType.  We don't want to create
    // a new view with LayoutStride, so replace it with the default layout
    // instead.
    using input_value    = typename InputView::non_const_value_type;
    using input_layout   = typename InputView::array_layout;
    using input_device   = typename InputView::device_type;
    using default_layout = typename input_device::execution_space::array_layout;
    using result_layout  =
      typename std::conditional<
        std::is_same< input_layout, Kokkos::LayoutStride >::value,
        default_layout,
        input_layout >::type;
    using type =
      Kokkos::DynRankView<input_value, result_layout, input_device>;
  };

}

//! Wrapper to simplify use of Sacado ViewFactory
template <typename InputViewType, typename CtorProp, typename ... Dims >
Impl::EnableIfViewFactoryArgs<InputViewType, CtorProp, false,
                              typename Impl::ResultDynRankView<InputViewType>::type>
createDynRankView(const InputViewType& a,
                  const CtorProp& prop,
                  const Dims... dims)
{
  using ResultViewType = typename Impl::ResultDynRankView<InputViewType>::type;
  return Sacado::createDynRankViewWithType<ResultViewType>(a, prop, dims...);
}

//! Wrapper to simplify use of Sacado ViewFactory for unmanaged views
//! (prop is a pointer or holds one), callable on device as well
template <typename InputViewType, typename CtorProp, typename ... Dims >
KOKKOS_INLINE_FUNCTION
Impl::EnableIfViewFactoryArgs<InputViewType, CtorProp, true,
                              typename Impl::ResultDynRankView<InputViewType>::type>
createDynRankView(const InputViewType& a,
                  const CtorProp& prop,
                  const Dims... dims)
{
  using ResultViewType = typename Impl::ResultDynRankView<InputViewType>::type;
  return Sacado::createDynRankViewWithType<ResultViewType>(a, prop, dims...);
}

//! Wrapper to simplify use of Sacado ViewFactory
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename ... Dims>
Impl::EnableIfViewFactoryArgs<InputViewType, CtorProp, false, ResultViewType>
createViewWithType(const InputViewType& a,
                   const CtorProp& prop,
                   const Dims... dims)
{
  using view_factory = Sacado::ViewFactory<InputViewType>;
  return view_factory::template create_view<ResultViewType>(a,prop,dims...);
}

//! Wrapper to simplify use of Sacado ViewFactory for unmanaged views
//! (prop is a pointer or holds one), callable on device as well
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename ... Dims>
KOKKOS_INLINE_FUNCTION
Impl::EnableIfViewFactoryArgs<InputViewType, CtorProp, true, ResultViewType>
createViewWithType(const InputViewType& a,
                   const CtorProp& prop,
                   const Dims... dims)
{
  using view_factory = Sacado::ViewFactory<InputViewType>;
  return view_factory::template create_view<ResultViewType>(a,prop,dims...);
}

}

#ifdef SACADO_ENABLE_DEPRECATED_CODE
namespace Kokkos {

// Deprecated wrapper function for createDynRankViewWithType
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename... Dims>
SACADO_DEPRECATED_WITH_COMMENT("Use Sacado::createDynRankViewWithType() instead of Kokkos::createDynRankViewWithType()")
typename std::enable_if<
  Kokkos::is_view<InputViewType>::value || Kokkos::is_dyn_rank_view<InputViewType>::value,
  ResultViewType>::type
createDynRankViewWithType(const InputViewType& a,
                          const CtorProp& prop,
                          const Dims... dims)
{
  return Sacado::createDynRankViewWithType<ResultViewType>(a, prop, dims...);
}

// Deprecated wrapper function for createDynRankView
template <typename InputViewType, typename CtorProp, typename... Dims>
SACADO_DEPRECATED_WITH_COMMENT("Use Sacado::createDynRankView() instead of Kokkos::createDynRankView()")
typename std::enable_if<
  Kokkos::is_view<InputViewType>::value || Kokkos::is_dyn_rank_view<InputViewType>::value,
  typename Sacado::Impl::ResultDynRankView<InputViewType>::type
  >::type
createDynRankView(const InputViewType& a,
                  const CtorProp& prop,
                  const Dims... dims)
{
  return Sacado::createDynRankView(a, prop, dims...);
}

// Deprecated wrapper function for createViewWithType
template <typename ResultViewType, typename InputViewType, typename CtorProp,
          typename... Dims>
SACADO_DEPRECATED_WITH_COMMENT("Use Sacado::createViewWithType() instead of Kokkos::createViewWithType()")
typename std::enable_if<
  Kokkos::is_view<InputViewType>::value || Kokkos::is_dyn_rank_view<InputViewType>::value,
  ResultViewType>::type
createViewWithType(const InputViewType& a,
                   const CtorProp& prop,
                   const Dims... dims)
{
  return Sacado::createViewWithType<ResultViewType>(a, prop, dims...);
}

} // namespace Kokkos
#endif

#endif /* #ifndef SACADO_FAD_KOKKOS_VIEW_FACTORY_HPP */
