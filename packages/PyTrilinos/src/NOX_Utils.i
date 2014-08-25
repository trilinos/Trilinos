// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2010) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// NOX_Utils is accessed from several places within PyTrilinos.NOX,
// and it also has a nested class issue to work around.  So I
// concentrate all of the NOX::Utils wrapper logic here.

#if SWIG_VERSION >= 0x030000
%feature("flatnested");
#else
// Handle the NOX::Utils:Fill and Sci nested classes by defining them
// exclusively for SWIG as though they were not nested.
namespace NOX
{
class Fill {
public:
  Fill(int ntimes, char ch);
  ~Fill();
  int n;
  char c;
};
%nestedworkaround Utils::Fill;

class Sci {
public:
  Sci(double val, int precision=-1);
  ~Sci();
  double d;
  int p;
};
%nestedworkaround Utils::Sci;
}
#endif

%{
#include "NOX_Utils.H"
%}

///////////////////////
// NOX Utils support //
///////////////////////
%rename(_print) NOX::Utils::print;
%include "NOX_Utils.H"

// SWIG thinks that Fill and Sci are un-nested NOX classes, so we
// need to trick the C++ compiler into understanding these so called
// un-nested NOX types.
%{
namespace NOX
{
typedef Utils::Fill Fill;
typedef Utils::Sci  Sci ;
}
%}
