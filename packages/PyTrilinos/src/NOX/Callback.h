// @HEADER
// ***********************************************************************
//
//                 PyTrilinos: Rapid Prototyping Package
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

#ifndef Callback_h
#define Callback_h

#include "Python.h"

class Callback
{
public:
  Callback();
  ~Callback();

  PyObject       * setFunction(PyObject * args);
  PyObject       * getFunction();
  const PyObject * getFunction() const;

private:
  // Private and not implemented so never can be called
  Callback(const Callback & a_ref);
  const Callback & operator = (const Callback & a_ref); 

  PyObject * mp_callback;
};

#endif //Callback_h
