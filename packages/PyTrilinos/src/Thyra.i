// -*- c++ -*-

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
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define %thyra_docstring
"
PyTrilinos.Thyra is the python interface to the Trilinos abstract
interface package Thyra:

    http://trilinos.sandia.gov/packages/thyra

The purpose of Thyra is to provide a set of interfaces and supporting
code that defines basic interoperability mechanisms between different
types of numerical software.  The python interface is currently
experimental, and at this time is inoperable.
"
%enddef

%module(package   = "PyTrilinos",
	docstring = %thyra_docstring) Thyra

// Code within the percent-bracket delimiters is copied verbatim to
// the C++ wrapper source file.  Anything that is %include-ed later
// needs to be #include-ed here.
%{
// System includes
#include <sstream>

// Configuration includes
#include "PyTrilinos_config.h"

// Thyra includes
#include "Teuchos_Range1D.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_OperatorVectorTypes.hpp"

// Namespace flattening
using Teuchos::RCP;
%}

// Ignore directives
%ignore *::operator=;

// Auto-documentation feature.
%feature("autodoc", "1");

// C++ STL support.
%include "stl.i"

namespace std
{
class logic_error;
}

// Thyra interface imports and includes
#define TEMPLATE_FRIENDS_NOT_SUPPORTED
%import  "Teuchos_ConfigDefs.hpp"
%import  "Teuchos_Range1D.hpp"
%import  "Teuchos_VerbosityLevel.hpp"
%import  "Teuchos_Describable.hpp"
%import  "Teuchos_TypeNameTraits.hpp"
%import  "Teuchos_RCPDecl.hpp"
%import  "RTOp_MPI_config.h"
%import  "Thyra_OperatorVectorTypes.hpp"
%import  "Thyra_ScalarProdBaseDecl.hpp"
%include "Thyra_VectorSpaceBaseDecl.hpp"
%include "Thyra_LinearOpBaseDecl.hpp"
%include "Thyra_MultiVectorBaseDecl.hpp"
%include "Thyra_VectorBaseDecl.hpp"
%include "Thyra_VectorSpaceDefaultBaseDecl.hpp"
%include "Thyra_ScalarProdVectorSpaceBaseDecl.hpp"
%include "Thyra_SpmdVectorSpaceBaseDecl.hpp"
%include "Thyra_SpmdVectorSpaceDefaultBaseDecl.hpp"
%include "Thyra_DefaultSpmdVectorSpaceDecl.hpp"
%include "Thyra_VectorStdOpsDecl.hpp"

// Macro for an interface, templated on type
%define %thyra_interface(type)

%ignore RCP<Thyra::VectorBase<type> >::access_node() const;
%ignore RCP<const Thyra::VectorSpaceBase<type> >::access_node() const;

%ignore Thyra::VectorSpaceBase<type>::createMember() const;
%ignore Thyra::VectorSpaceBase<type>::createMembers(int) const;
%ignore Thyra::VectorSpaceBase<type>::createMemberView(const RTOpPack::SubVectorView<type> &raw_v) const;
%ignore Thyra::VectorSpaceBase<type>::createMemberView(const RTOpPack::ConstSubVectorView<type> &raw_v) const;
%ignore Thyra::VectorSpaceBase<type>::createMembersView(const RTOpPack::SubMultiVectorView<type> &raw_mv) const;
%ignore Thyra::VectorSpaceBase<type>::createMembersView(const RTOpPack::ConstSubMultiVectorView<type> &raw_mv) const;

%feature("pythonappend") RCP<Thyra::VectorBase<type> >::
  RCP<Thyra::VectorBase<type> >
{print "Hello World!"}

// Teuchos templates
%template (HandleableVectorSpaceBase_ ## type) Teuchos::Handleable<Thyra::VectorSpaceBase<type> >;
%template (RCPVectorBase_ ## type)             RCP<Thyra::VectorBase<type> >;
%template (RCPVectorSpaceBase_ ## type)        RCP<const Thyra::VectorSpaceBase<type> >;

// Thyra templates
%template (VectorSpaceBase_ ## type)           	  Thyra::VectorSpaceBase<type>;
%template (LinearOpBase_ ## type)              	  Thyra::LinearOpBase<type,type>;
%template (MultiVectorBase_ ## type)           	  Thyra::MultiVectorBase<type>;
%template (VectorBase_ ## type)      	       	  Thyra::VectorBase<type>;
%template (VectorSpaceDefaultBase_ ## type)    	  Thyra::VectorSpaceDefaultBase<type>;
%template (ScalarProdVectorSpaceBase_ ## type) 	  Thyra::ScalarProdVectorSpaceBase<type>;
%template (SerialVectorSpaceBase_ ## type)     	  Thyra::SpmdVectorSpaceBase<type>;
%template (SerialVectorSpaceDefaultBase_ ## type) Thyra::SpmdVectorSpaceDefaultBase<type>;
%template (SerialVectorSpaceStd_ ## type)      	  Thyra::DefaultSpmdVectorSpace<type>;
%template (createMember_ ## type)              	  Thyra::createMember<type>;
%template (V_S_ ## type)                       	  Thyra::V_S<type>;
%template (sum_ ## type)                       	  Thyra::sum<type>;

%enddef

// Instantiations of interfaces
%thyra_interface(double)

// Extensions.



// Python code.
// %pythoncode
// {
// def createMember(arg):
//     rcpvb = createMember_double(arg)
//     d     = rcpvb.__deref__()
//     try:
//         rcpvb.this.append(d)
//     except SystemError:
//         rcpvb.this = [rcpvb.this, d]
//     return rcpvb
// }
