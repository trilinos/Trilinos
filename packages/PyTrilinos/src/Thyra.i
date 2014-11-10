// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
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
