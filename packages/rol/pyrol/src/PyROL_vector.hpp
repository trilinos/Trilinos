// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYROL_VECTOR
#define PYROL_VECTOR

#include <ROL_Ptr.hpp>
#include <pybind11/stl_bind.h>

#include <string>
#include <type_traits>
#include <vector>

namespace pyrol {

// Binder calls this for every std::vector type the bindings use (see
// "+binder std::vector" in PyROL_common.cfg). pybind11's bind_vector makes it
// a Python sequence that can be built from any iterable, and Python lists,
// tuples and NumPy arrays convert to it wherever ROL takes a std::vector.
// Vectors of numbers also expose their data to NumPy without a copy; such a
// view is valid only until the vector changes size (see pyrol/README.md). The
// C++-style methods keep code written for the previous bindings working.
template <typename T, typename Allocator>
void vector_binder(pybind11::module &m, std::string const &name, std::string const & /*allocator*/) {
  using Vector = std::vector<T, Allocator>;
  using Size = typename Vector::size_type;
  using Class = pybind11::class_<Vector, ROL::Ptr<Vector>>;

  // Binder's names: vector_double_t, vector_ROL_TimeStamp_double_t, ...
  std::string const suffix = name.size() > 2 && name.compare(name.size() - 2, 2, "_t") == 0 ? "" : "_t";
  std::string const class_name = "vector_" + name + suffix;
  Class cl = [&]() -> Class {
    if constexpr (std::is_arithmetic<T>::value && !std::is_same<T, bool>::value)
      return pybind11::bind_vector<Vector, ROL::Ptr<Vector>>(m, class_name, pybind11::buffer_protocol());
    else
      return pybind11::bind_vector<Vector, ROL::Ptr<Vector>>(m, class_name);
  }();
  pybind11::implicitly_convertible<pybind11::iterable, Vector>();

  if constexpr (std::is_default_constructible<T>::value) {
    cl.def(pybind11::init<Size>());
    cl.def("resize", [](Vector &v, Size n) { v.resize(n); });
  }
  cl.def(pybind11::init<Size, const T &>());
  cl.def("push_back", [](Vector &v, const T &x) { v.push_back(x); });
  cl.def("pop_back", [](Vector &v) { if (v.empty()) throw pybind11::index_error(); v.pop_back(); });
  cl.def("size", [](const Vector &v) { return v.size(); });
  cl.def("empty", [](const Vector &v) { return v.empty(); });
  cl.def("at", [](const Vector &v, Size i) -> T { return v.at(i); });
  cl.def("front", [](const Vector &v) -> T { if (v.empty()) throw pybind11::index_error(); return v.front(); });
  cl.def("back", [](const Vector &v) -> T { if (v.empty()) throw pybind11::index_error(); return v.back(); });
  cl.def("reserve", [](Vector &v, Size n) { v.reserve(n); });
  cl.def("capacity", [](const Vector &v) { return v.capacity(); });
  cl.def("swap", [](Vector &v, Vector &w) { v.swap(w); });
}

} // namespace pyrol

#endif // PYROL_VECTOR
