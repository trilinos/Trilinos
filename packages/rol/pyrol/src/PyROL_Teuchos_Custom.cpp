// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <PyROL_Teuchos_Custom.hpp>

// Implementation based on:
// https://github.com/trilinos/Trilinos/tree/master/packages/PyTrilinos/src/PyTrilinos_Teuchos_Util.cpp
bool setPythonParameter(Teuchos::ParameterList & plist,
			const std::string      & name,
			py::object             value)
{
  py::handle h = value;

  // Boolean values
  if (PyBool_Check(value.ptr ()))
  {
    plist.set(name, h.cast<bool>());
  }

  // Integer values
  else if (PyLong_Check(value.ptr ()))
  {
    plist.set(name, h.cast<int>());
  }

  // Floating point values
  else if (PyFloat_Check(value.ptr ()))
  {
    plist.set(name, h.cast<double>());
  }

  // Unicode values (stored as UTF-8)
  else if (PyUnicode_Check(value.ptr ()))
  {
    plist.set(name, h.cast<std::string>());
  }

  // Bytes values (stored as strings)
  else if (PyBytes_Check(value.ptr ()))
  {
    plist.set(name, h.cast<std::string>());
  }

  // All other value types, None included, are unsupported
  else
  {
    return false;
  }

  // Successful type conversion
  return true;
}    // setPythonParameter


// Implementation based on:
// https://github.com/trilinos/Trilinos/tree/master/packages/PyTrilinos/src/PyTrilinos_Teuchos_Util.cpp
py::object getPythonParameter(const Teuchos::ParameterList & plist,
			      const std::string            & name)
{
  // Get the parameter entry.  I now deal with the Teuchos::ParameterEntry
  // objects so that I can query the Teuchos::ParameterList without setting
  // the "used" flag to true.
  const Teuchos::ParameterEntry * entry = plist.getEntryPtr(name);
  if (entry == nullptr)
    throw py::key_error(name);
  // Boolean parameter values
  if (entry->isType< bool >())
  {
    bool value = Teuchos::any_cast< bool >(entry->getAny(false));
    return py::cast(value);
  }
  // Integer parameter values
  else if (entry->isType< int >())
  {
    int value = Teuchos::any_cast< int >(entry->getAny(false));
    return py::cast(value);
  }
  // Double parameter values
  else if (entry->isType< double >())
  {
    double value = Teuchos::any_cast< double >(entry->getAny(false));
    return py::cast(value);
  }
  // String parameter values
  else if (entry->isType< std::string >())
  {
    std::string value = Teuchos::any_cast< std::string >(entry->getAny(false));
    return py::cast(value.c_str());
  }
  // Char * parameter values
  else if (entry->isType< char * >())
  {
    char * value = Teuchos::any_cast< char * >(entry->getAny(false));
    return py::cast(value);
  }

  else if (entry->isArray())
  {
    try
    {
      Teuchos::Array< int > tArray =
        Teuchos::any_cast< Teuchos::Array< int > >(entry->getAny(false));
      return copyTeuchosArrayToNumPy(tArray);
    }
    catch(Teuchos::bad_any_cast &e)
    {
      try
      {
        Teuchos::Array< long > tArray =
          Teuchos::any_cast< Teuchos::Array< long > >(entry->getAny(false));
        return copyTeuchosArrayToNumPy(tArray);
      }
      catch(Teuchos::bad_any_cast &e)
      {
        try
        {
          Teuchos::Array< float > tArray =
            Teuchos::any_cast< Teuchos::Array< float > >(entry->getAny(false));
          return copyTeuchosArrayToNumPy(tArray);
        }
        catch(Teuchos::bad_any_cast &e)
        {
          try
          {
            Teuchos::Array< double > tArray =
              Teuchos::any_cast< Teuchos::Array< double > >(entry->getAny(false));
            return copyTeuchosArrayToNumPy(tArray);
          }
          catch(Teuchos::bad_any_cast &e)
          {
            return py::none();
          }
        }
      }
    }
  }

  // All  other types are unsupported
  return py::none();
}    // getPythonParameter
