// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYROL_TEUCHOS_CUSTOM
#define PYROL_TEUCHOS_CUSTOM

#include <Teuchos_ParameterList.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <ROL_Ptr.hpp>

namespace py = pybind11;

template<typename T>
pybind11::array_t<T> copyTeuchosArrayToNumPy(Teuchos::Array< T > & tArray) {

    pybind11::array_t<T> array(tArray.size());
    auto data = array.template mutable_unchecked<1>();
    for (int i=0; i < tArray.size(); ++i)
      data(i) = tArray[i];
    return array;
}

// Implementation based on:
// https://github.com/trilinos/Trilinos/tree/master/packages/PyTrilinos/src/PyTrilinos_Teuchos_Util.cpp
bool setPythonParameter(Teuchos::ParameterList & plist,
			const std::string      & name,
			py::object             value);

inline void setPythonParameterOrThrow(Teuchos::ParameterList & plist,
                                      const std::string & name,
                                      py::object value)
{
  if (!setPythonParameter(plist, name, value))
    throw py::type_error("Cannot store a value of type '" +
                         py::type::of(value).attr("__name__").cast<std::string>() +
                         "' in ParameterList entry '" + name +
                         "'; supported types are bool, int, float, str and ParameterList");
}


// Implementation based on:
// https://github.com/trilinos/Trilinos/tree/master/packages/PyTrilinos/src/PyTrilinos_Teuchos_Util.cpp
py::object getPythonParameter(const Teuchos::ParameterList & plist,
			      const std::string            & name);

// The sublist called name, as a pointer that also keeps plist alive.
inline Teuchos::RCP<Teuchos::ParameterList> sublistPtr(const Teuchos::RCP<Teuchos::ParameterList> &plist, const std::string &name) {
  return Teuchos::sublist(plist, name);
}

inline std::shared_ptr<Teuchos::ParameterList> sublistPtr(const std::shared_ptr<Teuchos::ParameterList> &plist, const std::string &name) {
  return std::shared_ptr<Teuchos::ParameterList>(plist, &plist->sublist(name));
}

template <typename T>
void def_ParameterList_member_functions(T cl) {
  using Ptr = ROL::Ptr<Teuchos::ParameterList>;
  cl.def("__setitem__", [](Ptr &m, const std::string &name, Teuchos::ParameterList value) { m->set(name,value);  });
  cl.def("set", [](Ptr &m, const std::string &name, Teuchos::ParameterList value) { m->set(name,value);  });
  cl.def("sublist", [](Ptr &m, const std::string &name) { if (!m->isSublist(name)) { throw py::key_error(name); } return pybind11::cast(sublistPtr(m, name)); }, pybind11::return_value_policy::reference);
  cl.def("__setitem__", [](Ptr &m, const std::string &name, pybind11::object value) { setPythonParameterOrThrow(*m,name,value);  });
  cl.def("__getitem__", [](Ptr &m, const std::string &name) {
    // Sublist
    if (m->isSublist(name))
      return py::cast(sublistPtr(m, name));
    return getPythonParameter(*m,name);
  });
  cl.def("set", [](Ptr &m, const std::string &name, pybind11::object value) { setPythonParameterOrThrow(*m,name,value);  });
  cl.def("get", [](Ptr &m, const std::string &name) {
    // Sublist
    if (m->isSublist(name))
      return py::cast(sublistPtr(m, name));
    return getPythonParameter(*m,name);
  });
  cl.def("__contains__", [](Ptr &m, const std::string &name) { return m->isParameter(name); });
  cl.def("__len__", [](Ptr &m) { return m->numParams(); });
}
#endif // PYROL_TEUCHOS_CUSTOM
