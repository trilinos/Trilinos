// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//           PyTrilinos.Anasazi: Python Interface to Anasazi
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

// This documentation string will be the python help facility help
// string
%define DOCSTRING
"The Anasazi module allows access to The Trilinos package
Anasazi.  Use the python help() facility for local documentation
on classes and methods, or see the on-line documentation for more
in-depth information."
%enddef

// Define the module name, its package and documentation string
%module(package="PyTrilinos", docstring=DOCSTRING) Anasazi

// Code within the percent-bracket delimiters is copied verbatim to
// the C++ wrapper source file.  Anything that is %include-ed later
// needs to be #include-ed here.
%{
// System includes
#include <sstream>

// Anasazi includes
#include "AnasaziVersion.cpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziOutputManager.hpp"
%}

// Ignore directives.  Here we use them to prevent wrapping the Print
// methods supported by Newp_Hello and Newp_Jambo.  Instead, we will
// later define __str__() methods for the classes, which is standard
// python technique for output.  The ignore directive is also good for
// eliminating warnings about methods python cannot wrap, such as
// operator=.
%ignore *::operator=;
%ignore *::print;

// Auto-documentation feature.  This ensures that calling help() on a
// wrapped method returns an argument list (or lists, in the case of
// overloaded methods).  The "1" option includes type information in
// the list.  While not as extensive as say, doxygen documentation,
// this is often enough to greatly increase the wrapper's usability.
%feature("autodoc", "1");

// Rename directives
%rename (OutputManager_) Anasazi::OutputManager;

// C++ STL support.  If the wrapped class uses standard template
// library containers, the following %include-s wraps the containers
// and makes certain conversions seamless, such as between std::string
// and python strings.
%include "std_string.i"

 //%import "Epetra.i"

// Anasazi interface includes.  Create a %include line for every
// header file with a prototype you want wrapped.  In this example,
// Anasazi_Version contains a function, Newp_Hello contains a
// class, and Newp_Jambo contains an optional class.
using namespace std;
%include "AnasaziVersion.cpp"
%include "AnasaziBasicSort.hpp"
%include "AnasaziOutputManager.hpp"

 //%template (OutputManager) Anasazi::OutputManager<double>;

// Extensions.  The %extend directive allows you to extend classes to
// include additional methods.  In this example, we are adding a
// __str__() method, which is a standard python method for classes.
// When the python "print" command encounters a non-string object, it
// calls the str() function, which in turn looks for a __str__()
// method in the class.  Thus if hello is a Newp_Hello object, then
// "print hello" will return a string obtained from the C++ Print()
// method.

// Python code.  This code is added to the end of the python proxy
// file created by swig.  Here we set the namespace attribute
// "__version__" equal to the value returned by the
// Anasazi_Version function.
%pythoncode %{

  __version__ = Anasazi_Version().split()[2]

%}
