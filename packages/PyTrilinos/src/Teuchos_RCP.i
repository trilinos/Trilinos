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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifdef HAVE_TEUCHOS

// The purpose of Teuchos_RCP.i is to provide a SWIG interface file
// that provides all of the infrastructure necessary for proper
// implementation of the Teuchos::RCP<> reference-counted pointer
// class.
%{
#include "Teuchos_RCP.hpp"
%}

// The Problem: Teuchos::RCP<> provides a mechanism for reference
// counting C++ objects, and preventing destruction until the
// reference count goes to zero.  However, all python objects are
// referenced counted already, so it is decidedly "unpythonic" to ask
// a python programmer to encapsulate a python object in a
// reference-counting wrapper.  We therefore want RCP handling to
// occur entirely "under the covers" and not be seen by the python
// programmer.  This means that when the C++ documentation calls for a
// Teuchos::RCP-wrapped object, the corresponding python class expects
// an unadorned python object of the analagous type, and the
// conversions in the compiled code are performed correctly.

// SWIG stores object instances by default using raw pointers.  To
// properly provide Teuchos::RCP<> wrapped objects to SWIG-wrapped
// functions or methods that require them, the object should be stored
// instead as a Teuchos::RCP<>*.  SWIG provides this alternate
// internal storage capability for boost::shared_ptr<> and
// std::shared_ptr<>.  We can modify this for use with Teuchos::RCP<>
// by implementing the following two #define-s:

#define SWIG_SHARED_PTR_NAMESPACE Teuchos
#define shared_ptr RCP

// and %include-ing the boost_shared_ptr.i file:

%include <boost_shared_ptr.i>

// This implements all of the boost::shared_ptr<> logic for
// Teuchos::RCP<>.  That is, the %shared_ptr() macro, when used, will
// apply to Teuchos::RCP<>.

// Most of the logic for implementing this infrastructure is in
// %typemap definitions to properly convert between unadorned python
// objects and Teuchos::RCP<> wrapped C++ objects.  Some of these
// typemaps need to be overridden in order to account for some
// differences in the Teuchos::RCP<> interface compared to the
// boost::shared_ptr<> interface.  This is done by introducing a
// %teuchos_rcp() macro that calls the %shared_ptr() macro and then
// implements all of the overrides.

%define %teuchos_rcp_typemaps_overrides(CONST, TYPE...)
// Output a plain pointer
%typemap(out) CONST TYPE *
{
  Teuchos::RCP< CONST TYPE > *smartresult = $1 ? new Teuchos::RCP< CONST TYPE >($1, bool($owner)): 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(Teuchos::RCP< TYPE > *),
				 $owner | SWIG_POINTER_OWN));
}
// Output a plain reference
%typemap(out) CONST TYPE &
{
  Teuchos::RCP< CONST TYPE > *smartresult = new Teuchos::RCP< CONST TYPE >($1, bool($owner));
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(Teuchos::RCP< TYPE > *),
				 SWIG_POINTER_OWN));
}
// Output a Teuchos::RCP< >
%typemap(out) Teuchos::RCP< CONST TYPE >
{
  Teuchos::RCP< CONST TYPE > *smartresult = $1.is_null() ? 0 : new Teuchos::RCP< CONST TYPE >($1);
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(Teuchos::RCP< TYPE > *),
				 SWIG_POINTER_OWN));
}
// Output a reference to a Teuchos::RCP< >
%typemap(out) Teuchos::RCP< CONST TYPE > &
{
  Teuchos::RCP< CONST TYPE > *smartresult = $1->is_null() ? 0 : new Teuchos::RCP< CONST TYPE >(*$1);
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor(Teuchos::RCP< TYPE > *),
				 SWIG_POINTER_OWN));
}
// Convert virtual method input argument to a plain pointer
%typemap(directorin) CONST TYPE *
{
  Teuchos::RCP< CONST TYPE > *temp$argnum = new Teuchos::RCP< CONST TYPE >($1_name, false);
  $input = SWIG_NewPointerObj(%as_voidptr(temp$argnum), $descriptor(Teuchos::RCP< TYPE > *), 0);
}
// Convert virtual method input argument to a plain reference
%typemap(directorin) CONST TYPE &
{
  Teuchos::RCP< CONST TYPE > *temp$argnum = new Teuchos::RCP< CONST TYPE >(&$1_name, false);
  $input = SWIG_NewPointerObj(%as_voidptr(temp$argnum), $descriptor(Teuchos::RCP< TYPE > *), 0);
}
// Convert virtual method python output to a plain pointer
%typemap(directorout) CONST TYPE * (void *argp = 0, int res = 0, Teuchos::RCP< CONST TYPE > temp)
{
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($1, &argp, $descriptor(Teuchos::RCP< TYPE > *),
			      %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res))
  {
    %dirout_fail(res, "$type"); 
  }
  if (newmem & SWIG_CAST_NEW_MEMORY)
  {
    temp = *%reinterpret_cast(argp, Teuchos::RCP< CONST TYPE > *);
    delete %reinterpret_cast(argp, Teuchos::RCP< CONST TYPE > *);
    $result = %const_cast(temp.get(), $1_ltype);
  }
  else
  {
    $result = argp ? %const_cast(%reinterpret_cast(argp, Teuchos::RCP< CONST TYPE > *)->get(),
				 $1_ltype) : 0;
  }
}
// Convert virtual method python output to a plain reference
%typemap(directorout) CONST TYPE &  (void *argp = 0, int res = 0, Teuchos::RCP< CONST TYPE > temp)
{
  int newmem = 0;
  res = SWIG_ConvertPtrAndOwn($1, &argp, $descriptor(Teuchos::RCP< TYPE > *),
			      %convertptr_flags, &newmem);
  if (!SWIG_IsOK(res))
  {
    %dirout_fail(res, "$type"); 
  }
  if (!argp)
  {
    %dirout_nullref("$type");
  }
  if (newmem & SWIG_CAST_NEW_MEMORY)
  {
    temp = *%reinterpret_cast(argp, Teuchos::RCP< CONST TYPE > *);
    delete %reinterpret_cast(argp, Teuchos::RCP< CONST TYPE > *);
    $result = %const_cast(temp.get(), $1_ltype);
  }
  else
  {
    $result = %const_cast(%reinterpret_cast(argp, Teuchos::RCP< CONST TYPE > *)->get(), $1_ltype);
  }
}
%enddef

// Define the %teuchos_rcp() macro to implement all of the typemaps
// that %shared_ptr() does, and then override the appropriate
// typemaps for both const and non-const versions.
%define %teuchos_rcp(CLASS...)
  %shared_ptr(CLASS)
  %teuchos_rcp_typemaps_overrides(SWIGEMPTYHACK, CLASS)
  %teuchos_rcp_typemaps_overrides(const, CLASS)
%enddef

#else

// If HAVE_TEUCHOS is not defined, then define %teuchos_rcp() to be
// an empty macro
%define %teuchos_rcp(CLASS...)
%enddef

#endif
