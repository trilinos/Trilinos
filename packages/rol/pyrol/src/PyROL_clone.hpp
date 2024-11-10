// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_RCP.hpp"
#include "pybind11/pybind11.h"

template <typename A, typename trampoline_A>
Teuchos::RCP<A> customClone(const trampoline_A *ptr_to_this, const std::string &class_name, const std::string &function_name) {
    pybind11::gil_scoped_acquire gil;
    pybind11::function overload = pybind11::get_overload(static_cast<const A *>(ptr_to_this), function_name.c_str());
    if (overload) {
        auto o = overload.operator()<pybind11::return_value_policy::reference>();
        auto self = pybind11::cast(ptr_to_this);
        auto cloned = self.attr(function_name.c_str())();

        auto keep_python_state_alive = Teuchos::rcp<pybind11::object>(new pybind11::object(cloned));
        auto ptr = pybind11::cast<trampoline_A *>(cloned);

        // aliasing shared_ptr: points to `trampoline_A* ptr` but refcounts the Python object
        return Teuchos::RCP<A>(keep_python_state_alive, ptr);		
    }
    const std::string error_message = "Tried to call pure virtual function \"customClone\" called from " + class_name + "::"+ function_name;
    pybind11::pybind11_fail(error_message.c_str());
}