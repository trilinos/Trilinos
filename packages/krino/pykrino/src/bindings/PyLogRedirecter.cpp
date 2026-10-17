// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <PyLogRedirecter.hpp>
#include <pybind11/pybind11.h>

namespace pykrino {

namespace py = pybind11;

PyLogRedirecter::PyLogRedirecter() : myRedirecter() {}

PyLogRedirecter::~PyLogRedirecter() {
  // Destructor will automatically restore original buffer
}

void PyLogRedirecter::clear() {
  myRedirecter.clear();
}

std::string PyLogRedirecter::get_log() const {
  return myRedirecter.get_log();
}

void bind_log_redirecter(py::module_ &m) {
  py::class_<PyLogRedirecter>(m, "LogRedirecter")
    .def(py::init<>(), "Create a LogRedirecter to capture krino log output")
    .def("clear", &PyLogRedirecter::clear, "Clear the captured log buffer")
    .def("get_log", &PyLogRedirecter::get_log, "Get the captured log output as a string");
}

} // namespace pykrino
