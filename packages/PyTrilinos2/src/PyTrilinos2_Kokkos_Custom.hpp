#ifndef PYTRILINOS2_KOKKOS_CUSTOM_HPP
#define PYTRILINOS2_KOKKOS_CUSTOM_HPP

#include "Kokkos_Core.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <torch/torch.h>
#include <torch/extension.h>


template <typename T>
pybind11::array_t<typename T::value_type>
view_as_numpy(T kokkos_array_host) {

  pybind11::array_t<typename T::value_type> result;

  T *c = new T(std::move(kokkos_array_host));
  auto capsule =
      pybind11::capsule(c, [](void *x) { delete reinterpret_cast<T *>(x); });

  if constexpr (T::rank == 1) {
    result = pybind11::array_t<typename T::value_type>(
        kokkos_array_host.extent(0), kokkos_array_host.data(), capsule);
  } else if constexpr (T::rank == 2) {
    if constexpr (std::is_same_v<typename T::array_layout, Kokkos::LayoutLeft>) {
      result = pybind11::array_t<typename T::value_type, pybind11::array::f_style | pybind11::array::forcecast>(
                                                                                                                {kokkos_array_host.extent(0), kokkos_array_host.extent(1)},
        kokkos_array_host.data(), capsule);
    } else if constexpr (std::is_same_v<typename T::array_layout, Kokkos::LayoutRight>) {
      result = pybind11::array_t<typename T::value_type, pybind11::array::c_style | pybind11::array::forcecast>(
                                                                                                                {kokkos_array_host.extent(0), kokkos_array_host.extent(1)},
        kokkos_array_host.data(), capsule);
    } else
      result = pybind11::array_t<typename T::value_type>(0);
  } else {
    result = pybind11::array_t<typename T::value_type>(0);
  }
  return result;
}


template <class T>
torch::Tensor view_as_torch(T kokkos_array) {
  torch::Tensor result;
  if constexpr (T::rank == 1)
    result = torch::from_blob(kokkos_array.data(), {kokkos_array.extent(0)}, {kokkos_array.stride(0)}, torch::CppTypeToScalarType<typename T::value_type>::value);
  else if constexpr (T::rank == 2)
    result = torch::from_blob(kokkos_array.data(), {kokkos_array.extent(0), kokkos_array.extent(1)}, {kokkos_array.stride(0), kokkos_array.stride(1)}, torch::CppTypeToScalarType<typename T::value_type>::value);
  return result;
}


template <typename T, typename V = void> constexpr bool is_incomplete = true;
template <typename T> constexpr bool is_incomplete<T, std::enable_if_t<sizeof(T)>> = false;

template <class DataType, class... Properties> class binder_Kokkos_View {
  using view_type = Kokkos::View<DataType, Properties...>;
  using holder_type = std::shared_ptr<view_type>;
  using py_class_type = pybind11::class_<view_type, holder_type>;

public:
  binder_Kokkos_View(pybind11::module &m, std::string const &name,
                     std::string const &prop_name) {
    py_class_type cl(m, ("View_" + name + prop_name).c_str());

    cl.def("extent", &view_type::extent);

    if constexpr (Kokkos::SpaceAccessibility<Kokkos::HostSpace, typename view_type::memory_space>::accessible)
      cl.def("numpy", &view_as_numpy<view_type>);

    if constexpr (!is_incomplete<torch::CppTypeToScalarType<typename view_type::value_type>>)
      cl.def("torch", &view_as_torch<view_type>);

    // cl.def("rank", &view_type::rank);
    // cl.def("rank_dynamic", &view_type::rank_dynamic);
    // cl.def("stride", &view_type::stride);
  }
};

#endif
