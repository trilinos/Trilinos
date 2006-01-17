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
// #include "Thyra_Config.h"
#include "Thyra_SerialVectorSpaceStd.hpp"
#include "Thyra_VectorBase.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Thyra_VectorStdOps.hpp"
%}

// Ignore directives
%ignore *::operator=;

// Auto-documentation feature.
%feature("autodoc", "1");

// C++ STL support.
%include "std_string.i"


// Thyra interface imports and includes
using namespace std;
#define TEMPLATE_FRIENDS_NOT_SUPPORTED
%import  "Teuchos_VerbosityLevel.hpp"
%import  "Teuchos_Describable.hpp"
%import  "RTOp_config.h"
%import  "Thyra_OperatorVectorTypes.hpp"
%import  "Thyra_ScalarProdBaseDecl.hpp"
%include "Thyra_VectorSpaceBaseDecl.hpp"
%include "Thyra_LinearOpBaseDecl.hpp"
%include "Thyra_MultiVectorBaseDecl.hpp"
%include "Thyra_VectorBaseDecl.hpp"
%include "Teuchos_RefCountPtrDecl.hpp"
%include "Thyra_VectorSpaceDefaultBaseDecl.hpp"
%include "Thyra_ScalarProdVectorSpaceBaseDecl.hpp"
%include "Thyra_SerialVectorSpaceBaseDecl.hpp"
%include "Thyra_SerialVectorSpaceStdDecl.hpp"
%include "Thyra_VectorStdOpsDecl.hpp"

// Macro for an interface, templated on type
%define INTERFACE(type)

%ignore Teuchos::RefCountPtr<Thyra::VectorBase<type> >::access_node() const;
%ignore Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<type> >::access_node() const;

%ignore Thyra::VectorSpaceBase<type>::createMember() const;
%ignore Thyra::VectorSpaceBase<type>::createMembers(int) const;
%ignore Thyra::VectorSpaceBase<type>::createMemberView(const RTOpPack::MutableSubVectorT<type> &raw_v) const;
%ignore Thyra::VectorSpaceBase<type>::createMemberView(const RTOpPack::SubVectorT<type> &raw_v) const;
%ignore Thyra::VectorSpaceBase<type>::createMembersView(const RTOpPack::MutableSubMultiVectorT<type> &raw_mv) const;
%ignore Thyra::VectorSpaceBase<type>::createMembersView(const RTOpPack::SubMultiVectorT<type> &raw_mv) const;

// This is the start of a typemap for converting RefCountPtr's for a
// VectorBase<double> to a raw pointer
// %typemap(in) Teuchos::RefCountPtr<Thyra::VectorBase<type> >* {
//   
// }

%template (RCPVectorBase_ ## type)      Teuchos::RefCountPtr<Thyra::VectorBase<type> >;
%template (RCPVectorSpaceBase_ ## type) Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<type> >;

%template (VectorSpaceBase_ ## type)           Thyra::VectorSpaceBase<type>;
%template (LinearOpBase_ ## type)              Thyra::LinearOpBase<type,type>;
%template (MultiVectorBase_ ## type)           Thyra::MultiVectorBase<type>;
%template (VectorBase_ ## type)      	       Thyra::VectorBase<type>;
%template (VectorSpaceDefaultBase_ ## type)    Thyra::VectorSpaceDefaultBase<type>;
%template (ScalarProdVectorSpaceBase_ ## type) Thyra::ScalarProdVectorSpaceBase<type>;
%template (SerialVectorSpaceBase_ ## type)     Thyra::SerialVectorSpaceBase<type>;
%template (SerialVectorSpaceStd_ ## type)      Thyra::SerialVectorSpaceStd<type>;
%template (createMember_ ## type)              Thyra::createMember<type>;
%template (V_S_ ## type)                       Thyra::V_S<type>;
%template (sum_ ## type)                       Thyra::sum<type>;

%enddef

// Instantiations of interfaces
INTERFACE(double)

// Extensions.

// Python code.
