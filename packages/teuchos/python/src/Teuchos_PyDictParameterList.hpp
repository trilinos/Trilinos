// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef PYDICTPARAMETERLIST_H
#define PYDICTPARAMETERLIST_H

#include "Python.h"

#include "Teuchos_ParameterList.hpp"

class PyDictParameterList : public Teuchos::ParameterList {

public:
  PyDictParameterList();
  PyDictParameterList(PyObject * dict);
  PyDictParameterList(const PyDictParameterList & source);
  ~PyDictParameterList();

  // Python methods that provide PyDict-like behavior
  int        __cmp__(PyObject * dict) const;
  int        __cmp__(const PyDictParameterList) const;
  bool       __contains__(const char * key) const;
  void       __delitem__(const char * key);
  bool       __eq__(PyObject * dict) const;
  bool       __eq__(const PyDictParameterList) const;
  bool       __ge__(PyObject * dict) const;
  bool       __ge__(const PyDictParameterList) const;
  PyObject * __getattribute__(char * name) const;
  PyObject * __getitem__(const char * key) const;
  bool       __gt__(PyObject * dict) const;
  bool       __gt__(const PyDictParameterList) const;
  PyObject * __iter__() const;
  bool       __le__(PyObject * dict) const;
  bool       __le__(const PyDictParameterList) const;
  int        __len__() const;
  bool       __lt__(PyObject * dict) const;
  bool       __lt__(const PyDictParameterList) const;
  bool       __ne__(PyObject * dict) const;
  bool       __ne__(const PyDictParameterList) const;
  char *     __repr__() const;
  int        __setitem__(const char * key, PyObject * value);
  char *     __str__() const;
  void       clear();
  PyObject * copy() const;
  PyObject * get(const char * key, PyObject * value = Py_None) const;
  bool       has_key(const char * key) const;
  PyObject * items() const;
  PyObject * iteritems() const;
  PyObject * iterkeys() const;
  PyObject * itervalues() const;
  PyObject * keys() const;
  PyObject * pop(const char * key, PyObject * default = Py_None);
  PyObject * popitem();
  void       setdefault(const char * key, PyObject * value = Py_None);
  void       update(PyObject * dict);
  void       update(const PyDictParameterList &);
  PyObject * values() const;

private:
  PyObject * __pyDict__;

};

#endif
