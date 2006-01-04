// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//             PyTrilinos.Thyra: Python Interface to Thyra
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
"The Thyra module allows access to The Trilinos package
Thyra.  Use the python help() facility for local documentation
on classes and methods, or see the on-line documentation for more
in-depth information."
%enddef

// Define the module name, its package and documentation string
%module(package="PyTrilinos", docstring=DOCSTRING) Thyra

// Code within the percent-bracket delimiters is copied verbatim to
// the C++ wrapper source file.  Anything that is %include-ed later
// needs to be #include-ed here.
%{
// System includes
#include <sstream>

// Thyra includes
#include "Thyra_Config.h"
#include "Thyra_SerialVectorSpaceStd.hpp"
%}

// Ignore directives
%ignore Thyra::VectorSpaceBase<double>::createMember() const;
%ignore Thyra::VectorSpaceBase<double>::createMembers(int) const;
%ignore Thyra::VectorSpaceBase<double>::createMemberView(const RTOpPack::MutableSubVectorT<double> &raw_v) const;
%ignore Thyra::VectorSpaceBase<double>::createMemberView(const RTOpPack::SubVectorT<double> &raw_v) const;
%ignore Thyra::VectorSpaceBase<double>::createMembersView(const RTOpPack::MutableSubMultiVectorT<double> &raw_mv) const;
%ignore Thyra::VectorSpaceBase<double>::createMembersView(const RTOpPack::SubMultiVectorT<double> &raw_mv) const;

// Auto-documentation feature.
%feature("autodoc", "1");

// C++ STL support.
%include "std_string.i"

// Thyra interface includes.
using namespace std;
#define TEMPLATE_FRIENDS_NOT_SUPPORTED
%import  "Teuchos_Describable.hpp"
%import  "RTOp_config.h"
%import  "Thyra_OperatorVectorTypes.hpp"
%import  "Thyra_ScalarProdBaseDecl.hpp"
%include "Thyra_VectorSpaceBaseDecl.hpp"
%include "Thyra_VectorSpaceDefaultBaseDecl.hpp"
%include "Thyra_ScalarProdVectorSpaceBaseDecl.hpp"
%include "Thyra_SerialVectorSpaceBaseDecl.hpp"
%include "Thyra_SerialVectorSpaceStdDecl.hpp"

using namespace Thyra;
%template (VectorSpaceBaseDouble)           Thyra::VectorSpaceBase<double>;
%template (VectorSpaceDefaultBaseDouble)    Thyra::VectorSpaceDefaultBase<double>;
%template (ScalarProdVectorSpaceBaseDouble) Thyra::ScalarProdVectorSpaceBase<double>;
%template (SerialVectorSpaceBaseDouble)     Thyra::SerialVectorSpaceBase<double>;
%template (SerialVectorSpaceStdDouble)      Thyra::SerialVectorSpaceStd<double>;

%template (createMemberDouble) Thyra::createMember<double>;

// Extensions.

// Python code.
