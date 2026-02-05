// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_KOKKOS_VIEW_SUPPORT_HPP
#define SACADO_FAD_KOKKOS_VIEW_SUPPORT_HPP

#include "Sacado_ConfigDefs.h"
#if defined(HAVE_SACADO_KOKKOS)

// Only include forward declarations so any overloads appear before they
// might be used inside Kokkos
#include "Kokkos_View_Fad_Fwd.hpp"

#include "Sacado_Traits.hpp"
#include <Kokkos_DynRankView.hpp>
#include <Sacado_Fad_Ops_Fwd.hpp>

// Rename this file
#include <Kokkos_LayoutContiguous.hpp>

// ====================================================================
// Kokkos customization points for the mdspan based View implementation
// ====================================================================
// - FadAccessor: mdspan accessor returning ViewFadType
// - customize_view_arguments: inject accessor into View for FadTypes
// - allocation_size_from_mapping_and_accessor: compute allocation size
// - accessor_from_mapping_and_accessor_arg: construct accessor
//
// These functions are injected into Kokkos implementation via ADL

namespace Sacado {

template <typename T, unsigned Stride = 0> struct LocalScalarType;

namespace Impl {
KOKKOS_INLINE_FUNCTION
constexpr unsigned computeFadPartitionSize(unsigned size, unsigned stride) {
  return ((size + stride - 1) / stride) == (size / stride)
             ? ((size + stride - 1) / stride)
             : 0;
}
} // namespace Impl

namespace Fad {

namespace Exp {

// PartitionedFadStride > 0 implies LayoutContiguous
template <class ElementType, class MemorySpace, size_t FadStaticStride,
          size_t PartitionedFadStride, bool IsUnmanaged>
class FadAccessor {
  using fad_type = ElementType;
  using fad_value_type = typename Sacado::ValueType<fad_type>::type;
  constexpr static size_t FadStaticDimension =
      Sacado::StaticSize<fad_type>::value;

public:
  using element_type = ElementType;
  using data_handle_type =
      std::conditional_t<IsUnmanaged,
      fad_value_type*,
      Kokkos::Impl::ReferenceCountedDataHandle<fad_value_type, MemorySpace>>;

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) &&                                  \
    (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
  // avoid division by zero later
  constexpr static size_t partitioned_fad_stride =
      PartitionedFadStride > 0 ? PartitionedFadStride : 1;
  // The partitioned static size -- this will be 0 if PartitionedFadStride
  // does not evenly divide FadStaticDimension
  constexpr static size_t PartitionedFadStaticDimension =
      Sacado::Impl::computeFadPartitionSize(FadStaticDimension,
                                            partitioned_fad_stride);

#if defined(KOKKOS_ENABLE_CUDA)
  typedef typename Sacado::LocalScalarType<
      fad_type, unsigned(partitioned_fad_stride)>::type strided_scalar_type;
  typedef typename std::conditional_t<
      std::is_same<typename MemorySpace::execution_space, Kokkos::Cuda>::value,
      strided_scalar_type, fad_type>
      thread_local_scalar_type;
#elif defined(KOKKOS_ENABLE_HIP)
  typedef typename Sacado::LocalScalarType<
      fad_type, unsigned(partitioned_fad_stride)>::type strided_scalar_type;
  typedef typename std::conditional_t<
      std::is_same<typename MemorySpace::execution_space, Kokkos::HIP>::value,
      strided_scalar_type, fad_type>
      thread_local_scalar_type;
#else
  typedef fad_type thread_local_scalar_type;
#endif
  using reference = std::conditional_t<
      (PartitionedFadStride > 0),
      typename Sacado::ViewFadType<thread_local_scalar_type,
                                   PartitionedFadStaticDimension, 0>::type,
      typename Sacado::ViewFadType<fad_type, FadStaticDimension,
                                   FadStaticStride>::type>;
#else
  using reference = typename Sacado::ViewFadType<fad_type, FadStaticDimension,
                                                 FadStaticStride>::type;
#endif

  using offset_policy = FadAccessor<ElementType, MemorySpace, FadStaticStride,
                                    PartitionedFadStride, IsUnmanaged>;
  using memory_space = MemorySpace;

  using scalar_type = fad_value_type;

  typedef Sacado::integral_nonzero<unsigned, FadStaticDimension>
      sacado_size_type;
  typedef Sacado::integral_nonzero<unsigned, FadStaticStride>
      sacado_stride_type;

  sacado_size_type m_fad_size = {};
  sacado_stride_type m_fad_stride = {};

  KOKKOS_FUNCTION
  constexpr auto fad_size() const { return m_fad_size; }

  KOKKOS_DEFAULTED_FUNCTION
  constexpr FadAccessor() noexcept = default;

  KOKKOS_FUNCTION
  constexpr FadAccessor(size_t fad_size, size_t fad_stride) noexcept {
    if constexpr (FadStaticDimension == 0) {
      m_fad_size = fad_size - 1;
    }
    if constexpr (FadStaticStride == 0) {
      m_fad_stride = fad_stride;
    }
  }

  template <class OtherElementType, class OtherSpace,
            size_t OtherFadStaticStride, size_t OtherPartitionedFadStride,
	    bool OtherIsUnmanaged,
            class =
	       std::enable_if_t<std::is_same_v<std::remove_cv_t<ElementType>,
	                        std::remove_cv_t<OtherElementType>>>
		// In ISO C++ we generally don't allow const->non-const conversion
		//std::enable_if_t<std::is_convertible_v<
                //                 OtherElementType (*)[], element_type (*)[]>>
		    >
  KOKKOS_FUNCTION constexpr FadAccessor(
      const FadAccessor<OtherElementType, OtherSpace, OtherFadStaticStride,
                        OtherPartitionedFadStride, OtherIsUnmanaged> &other) {
    m_fad_size = other.m_fad_size;
    m_fad_stride = (OtherPartitionedFadStride == 0)
                       ? static_cast<sacado_stride_type>(other.m_fad_stride)
                       : static_cast<sacado_stride_type>(1);
  }

private:
  KOKKOS_FUNCTION 
  fad_value_type* get_ptr(const data_handle_type &p) const {
    if constexpr (IsUnmanaged) {
      return p;
    } else {
      return p.get();
    }
  }
public:

  KOKKOS_FUNCTION
  constexpr reference access(const data_handle_type &p, size_t i) const {
    if constexpr (PartitionedFadStride == 0)
      return reference(
          get_ptr(p) + i * (m_fad_stride.value <= 1 ? m_fad_size.value + 1 : 1),
          m_fad_size.value, m_fad_stride.value != 0 ? m_fad_stride.value : 1);
    else {
      size_t base_offset = i * (m_fad_size.value + 1);
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) &&                                  \
    (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
      return reference(get_ptr(p) + base_offset + threadIdx.x,
                       get_ptr(p) + base_offset + m_fad_size.value,
                       (m_fad_size.value + blockDim.x - threadIdx.x - 1) /
                           blockDim.x,
                       blockDim.x);
#else
      return reference(get_ptr(p) + base_offset,
                       get_ptr(p) + base_offset + m_fad_size.value,
                       m_fad_size.value, 1);
#endif
    }
  }

  KOKKOS_FUNCTION
  constexpr data_handle_type offset(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type &p,
#else
      // FIXME OpenACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const {
    if constexpr (std::is_pointer_v<data_handle_type>) {
      if constexpr (PartitionedFadStride != 0) {
        return data_handle_type{get_ptr(p) + i * (m_fad_size.value + 1)};
      } else {
        return data_handle_type{
          get_ptr(p) + i * (m_fad_stride.value <= 1 ? m_fad_size.value + 1 : 1)};
      }
    } else {
      if constexpr (PartitionedFadStride != 0)
        return data_handle_type{p, get_ptr(p) + i * (m_fad_size.value + 1)};
      else
        return data_handle_type{
          p,
          get_ptr(p) + i * (m_fad_stride.value <= 1 ? m_fad_size.value + 1 : 1)};
    }
  }
};

template <class MappingType, class ElementType, class MemorySpace,
          size_t FadStaticStride, size_t PartitionedFadStride, bool IsUnmanaged>
KOKKOS_INLINE_FUNCTION size_t allocation_size_from_mapping_and_accessor(
    const MappingType &mapping,
    const FadAccessor<ElementType, MemorySpace, FadStaticStride,
                      PartitionedFadStride, IsUnmanaged> &acc) {
  size_t element_size = acc.m_fad_size.value + 1;
  return mapping.required_span_size() * element_size;
}

template <class StorageType> struct GeneralFad;

template <class T, class LayoutType, class DeviceType, class MemoryTraits>
KOKKOS_INLINE_FUNCTION constexpr auto customize_view_arguments(
    Kokkos::Impl::ViewArguments<GeneralFad<T>, LayoutType, DeviceType,
                                MemoryTraits>) {
  constexpr int static_stride =
      std::is_same_v<Kokkos::LayoutRight, LayoutType> ? 1 : 0;
  constexpr size_t partitioned_fad_stride = []() {
    if constexpr (Kokkos::is_layout_contiguous<LayoutType>::value)
      return LayoutType::scalar_stride;
    else
      return 0;
  }();

  return Kokkos::Impl::ViewCustomArguments<
      size_t, FadAccessor<GeneralFad<T>, typename DeviceType::memory_space,
                          static_stride, partitioned_fad_stride, MemoryTraits::is_unmanaged>>();
}

template <class T, class LayoutType, class DeviceType, class MemoryTraits>
KOKKOS_INLINE_FUNCTION constexpr auto customize_view_arguments(
    Kokkos::Impl::ViewArguments<const GeneralFad<T>, LayoutType, DeviceType,
                                MemoryTraits>) {
  constexpr int static_stride =
      std::is_same_v<Kokkos::LayoutRight, LayoutType> ? 1 : 0;
  constexpr size_t partitioned_fad_stride = []() {
    if constexpr (Kokkos::is_layout_contiguous<LayoutType>::value)
      return LayoutType::scalar_stride;
    else
      return 0;
  }();

  return Kokkos::Impl::ViewCustomArguments<
      size_t,
      FadAccessor<const GeneralFad<T>, typename DeviceType::memory_space,
                  static_stride, partitioned_fad_stride, MemoryTraits::is_unmanaged>>();
}

template <class MappingType, class ElementType, class MemorySpace,
          size_t FadStaticStride, size_t PartitionedFadStride, bool IsUnmanaged>
KOKKOS_INLINE_FUNCTION auto accessor_from_mapping_and_accessor_arg(
    const Kokkos::Impl::AccessorTypeTag<FadAccessor<
        ElementType, MemorySpace, FadStaticStride, PartitionedFadStride, IsUnmanaged>> &,
    const MappingType &mapping,
    const Kokkos::Impl::AccessorArg_t &accessor_arg) {
  using fad_acc_t = FadAccessor<ElementType, MemorySpace, FadStaticStride,
                                PartitionedFadStride, IsUnmanaged>;
  return fad_acc_t(accessor_arg.value, mapping.required_span_size());
}

} // namespace Exp
} // namespace Fad
} // namespace Sacado

// =======================================================================
// Sacado View helpers
// =======================================================================

namespace Sacado {
// Whether a given type is a view with Sacado FAD scalar type
template <typename view_type> struct is_view_fad {
  static const bool value = false;
};

template <typename T, typename... P> struct is_view_fad<Kokkos::View<T, P...>> {
  typedef Kokkos::View<T, P...> view_type;
  static const bool value =
      Sacado::IsFad<typename view_type::value_type>::value;
};

template <class... ViewArgs>
struct is_view_fad<Kokkos::DynRankView<ViewArgs...>> {
  constexpr static bool value =
      is_view_fad<typename Kokkos::DynRankView<ViewArgs...>::view_type>::value;
};

template <typename view_type> struct is_view_fad_contiguous {
  static const bool value =
      is_view_fad<view_type>::value &&
      Kokkos::is_layout_contiguous<typename view_type::array_layout>::value;
};

template <typename ViewType> struct is_dynrankview_fad_contiguous {
  static const bool value = false;
};

template <class... Args>
struct is_dynrankview_fad_contiguous<Kokkos::DynRankView<Args...>> {
  using view_type = Kokkos::DynRankView<Args...>;
  static const bool value =
      is_view_fad<view_type>::value &&
      Kokkos::is_layout_contiguous<typename view_type::array_layout>::value;
};

template <typename ViewType> struct ViewScalarStride {
  static constexpr unsigned stride = 1;
  static constexpr bool is_unit_stride = 1;
};

template <class T, class L, unsigned S, class... Args>
struct ViewScalarStride<
    Kokkos::View<T, Kokkos::LayoutContiguous<L, S>, Args...>> {
  static constexpr unsigned stride = S;
  static constexpr bool is_unit_stride = (stride == 1u);
};

template <class T, class L, unsigned S, class... Args>
struct ViewScalarStride<
    Kokkos::DynRankView<T, Kokkos::LayoutContiguous<L, S>, Args...>> {
  static constexpr unsigned stride = S;
  static constexpr bool is_unit_stride = (stride == 1u);
};

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto
as_scalar_view(const Kokkos::View<DataType, Properties...> &view) {
  using view_t = Kokkos::View<DataType, Properties...>;
  if constexpr (Kokkos::is_view_fad<view_t>::value) {
    using value_type = typename view_t::value_type::value_type;
    return Kokkos::View<value_type *, Properties...>(
        view.data(),
        view.mapping().required_span_size() * view.accessor().fad_size());
  } else {
    return view;
  }
}

KOKKOS_INLINE_FUNCTION size_t dimension_scalar() { return 0; }

template <class View>
KOKKOS_FUNCTION size_t dimension_scalar(const View &view) {
  if constexpr (Kokkos::is_view_fad<View>::value) {
    return static_cast<size_t>(view.accessor().fad_size() + 1);
  } else {
    return 0;
  }
}

template <class... Views>
KOKKOS_FUNCTION size_t dimension_scalar(const Views &...views) {
  return Kokkos::max(dimension_scalar(views)...);
}

template<class T>
KOKKOS_FUNCTION
auto data_address_of(T& val) { return &val; }

template<class T>
KOKKOS_FUNCTION
auto data_address_of(Fad::Exp::GeneralFad<T>& val) { return static_cast<int>(val.size()) > 0 ? &val.fastAccessDx(0) : &val.val(); }


template <class... ViewArgs, class... OtherViews>
KOKKOS_FUNCTION
auto common_view_alloc_prop(const Kokkos::View<ViewArgs...> &view,
                            OtherViews... views) {
  constexpr bool any_fad_view =
      Kokkos::is_view_fad<Kokkos::View<ViewArgs...>>::value ||
      (Kokkos::is_view_fad<OtherViews>::value || ... || false);
  if constexpr (any_fad_view) {
    if constexpr (sizeof...(OtherViews) == 0) {
      return Kokkos::Impl::AccessorArg_t{
          static_cast<size_t>(Sacado::dimension_scalar(view))};
    } else {
      return Kokkos::Impl::AccessorArg_t{static_cast<size_t>(Kokkos::max(
          Sacado::dimension_scalar(view), Sacado::dimension_scalar(views)...))};
    }
  } else {
    using value_type =
        std::common_type_t<typename Kokkos::View<ViewArgs...>::value_type,
                           typename OtherViews::value_type...>;
    return Kokkos::Impl::CommonViewAllocProp<void, value_type>();
  }
}
template <class... ViewArgs, class... OtherViews>
KOKKOS_FUNCTION
auto common_view_alloc_prop(const Kokkos::DynRankView<ViewArgs...> &view,
                            OtherViews... views) {
  constexpr bool any_fad_view =
      Kokkos::is_view_fad<Kokkos::DynRankView<ViewArgs...>>::value ||
      (Kokkos::is_view_fad<OtherViews>::value || ... || false);
  if constexpr (any_fad_view) {
    if constexpr ((sizeof ... (OtherViews)) == 0) {
      return Kokkos::Impl::AccessorArg_t{static_cast<size_t>(Sacado::dimension_scalar(view))};
    } else {
      return Kokkos::Impl::AccessorArg_t{static_cast<size_t>(Kokkos::max(
          Sacado::dimension_scalar(view), Sacado::dimension_scalar(views)...))};
    }
  } else {
    using value_type =
        std::common_type_t<typename Kokkos::View<ViewArgs...>::value_type,
                           typename OtherViews::value_type...>;
    return Kokkos::Impl::CommonViewAllocProp<void, value_type>();
  }
}

} // namespace Sacado

#define SACADO_FAD_KOKKOS_VIEW_SUPPORT_INCLUDES
#include "Sacado_Fad_Kokkos_ThreadLocalScalar.hpp"
#include "Sacado_Fad_Kokkos_Specialization.hpp"

#if defined(HAVE_SACADO_KOKKOS) && defined(HAVE_SACADO_TEUCHOSKOKKOSCOMM) &&   \
    defined(HAVE_SACADO_VIEW_SPEC) && !defined(SACADO_DISABLE_FAD_VIEW_SPEC)

#include "Sacado_Fad_Kokkos_TeuchosComm.hpp"
#endif
#undef SACADO_FAD_KOKKOS_VIEW_SUPPORT_INCLUDES

namespace Kokkos {

using Sacado::common_view_alloc_prop;
using Sacado::dimension_scalar;
using Sacado::is_dynrankview_fad_contiguous;
using Sacado::is_view_fad;
using Sacado::is_view_fad_contiguous;
using Sacado::ThreadLocalScalarType;
using Sacado::ViewScalarStride;
namespace Impl {
using Sacado::Impl::computeFadPartitionSize;
}
} // namespace Kokkos

namespace std {
  template<class T, class T2>
  struct common_type<Sacado::Fad::Exp::GeneralFad<T>, T2> {
    using type = typename Sacado::Promote<Sacado::Fad::Exp::GeneralFad<T>, T2>::type;
  };
  template<class T, class T2>
  struct common_type<const Sacado::Fad::Exp::GeneralFad<T>, T2> {
    using type = typename Sacado::Promote<const Sacado::Fad::Exp::GeneralFad<T>, T2>::type;
  };
  template<class T1, class T>
  struct common_type<T1, Sacado::Fad::Exp::GeneralFad<T>> {
    using type = typename Sacado::Promote<Sacado::Fad::Exp::GeneralFad<T>, T1>::type;
  };
  template<class T1, class T>
  struct common_type<T1, const Sacado::Fad::Exp::GeneralFad<T>> {
    using type = typename Sacado::Promote<const Sacado::Fad::Exp::GeneralFad<T>, T1>::type;
  };
  template<class T1, class T2>
  struct common_type<Sacado::Fad::Exp::GeneralFad<T1>, Sacado::Fad::Exp::GeneralFad<T2>> {
    using type = typename Sacado::Promote<Sacado::Fad::Exp::GeneralFad<T1>, Sacado::Fad::Exp::GeneralFad<T2>>::type;
  };
  template<class T1, class T2>
  struct common_type<const Sacado::Fad::Exp::GeneralFad<T1>, Sacado::Fad::Exp::GeneralFad<T2>> {
    using type = typename Sacado::Promote<const Sacado::Fad::Exp::GeneralFad<T1>, Sacado::Fad::Exp::GeneralFad<T2>>::type;
  };
  template<class T1, class T2>
  struct common_type<Sacado::Fad::Exp::GeneralFad<T1>, const Sacado::Fad::Exp::GeneralFad<T2>> {
    using type = typename Sacado::Promote<Sacado::Fad::Exp::GeneralFad<T1>, const Sacado::Fad::Exp::GeneralFad<T2>>::type;
  };
  template<class T1, class T2>
  struct common_type<const Sacado::Fad::Exp::GeneralFad<T1>, const Sacado::Fad::Exp::GeneralFad<T2>> {
    using type = typename Sacado::Promote<const Sacado::Fad::Exp::GeneralFad<T1>, const Sacado::Fad::Exp::GeneralFad<T2>>::type;
  };
}
#endif

#endif // SACADO_FAD_KOKKOS_VIEW_SUPPORT_HPP
