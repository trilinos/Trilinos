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

// Set SHARED_PTR_DISOWN to $disown if required, for example
// #define SHARED_PTR_DISOWN $disown
#if !defined(SHARED_PTR_DISOWN)
#define SHARED_PTR_DISOWN 0
#endif

%{
#include "Teuchos_DefaultComm.hpp"
%}

/* %fragment("SWIG_null_deleter_python", "header", fragment="SWIG_null_deleter") { */
/* %#define SWIG_NO_NULL_DELETER_SWIG_BUILTIN_INIT */
/* } */

// Language specific macro implementing all the customisations for handling the smart pointer
%define %teuchos_rcp_dap_typemaps(CONST, CONVERTER, CLASS...)

// Mark CLASS as a smart pointer
%feature("smartptr", noblock=1) CLASS {Teuchos::RCP< CLASS >}

// %naturalvar is as documented for member variables
%naturalvar CLASS;
%naturalvar Teuchos::RCP< CONST CLASS >;

// destructor wrapper customisation
%feature("unref") CLASS "(void)arg1; delete smartarg1;"

// Typemap customizations...

// plain value
%typemap(in) CONST CLASS
{
  int newmem = 0;
  Teuchos::RCP< CLASS > * smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  $1 = %const_cast(*smartarg->get(), $1_ltype);
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(out) CONST CLASS
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    new Teuchos::RCP< CONST CLASS >(new $1_ltype(($1_ltype &)$1));
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                 $descriptor(Teuchos::RCP< CLASS > *),
                                 SWIG_POINTER_OWN));
}

%typemap(varin) CONST CLASS
{
  int newmem = 0;
  Teuchos::RCP< CLASS > * smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  $1 = %const_cast(*smartarg->get(), $1_ltype);
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(varout) CONST CLASS
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    new Teuchos::RCP< CONST CLASS >(new $1_ltype(($1_ltype &)$1));
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                    $descriptor(Teuchos::RCP< CLASS > *),
                                    SWIG_POINTER_OWN));
}

// plain pointer
// Note: $disown not implemented by default as it will lead to a
// memory leak of the RCP instance
%typemap(in) CONST CLASS * (Teuchos::RCP< CLASS > * smartarg = 0)
{
  int newmem = 0;
  smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  $1 = %const_cast(smartarg->get(), $1_ltype);
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(out, fragment="SWIG_null_deleter_python") CONST CLASS *
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    $1 ? new Teuchos::RCP< CONST CLASS >($1 SWIG_NO_NULL_DELETER_$owner) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                 $descriptor(Teuchos::RCP< CLASS > *),
                                 $owner | SWIG_POINTER_OWN));
}

%typemap(varin) CONST CLASS *
{
  int newmem = 0;
  Teuchos::RCP< CLASS > * smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  $1 = %const_cast(smartarg->get(), $1_ltype);
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(varout, fragment="SWIG_null_deleter_python") CONST CLASS *
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    $1 ? new Teuchos::RCP< CONST CLASS >($1 SWIG_NO_NULL_DELETER_0) : 0;
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                    $descriptor(Teuchos::RCP< CLASS > *),
                                    SWIG_POINTER_OWN));
}

// plain reference
%typemap(in) CONST CLASS & (Teuchos::RCP< CLASS > * smartarg = 0)
{
  int newmem = 0;
  smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  $1 = %const_cast(smartarg->get(), $1_ltype);
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(out, fragment="SWIG_null_deleter_python") CONST CLASS &
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    new Teuchos::RCP< CONST CLASS >($1 SWIG_NO_NULL_DELETER_$owner);
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                 $descriptor(Teuchos::RCP< CLASS > *),
                                 SWIG_POINTER_OWN));
}

%typemap(varin) CONST CLASS &
{
  int newmem = 0;
  Teuchos::RCP< CLASS > * smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  $1 = %const_cast(smartarg->get(), $1_ltype);
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(varout, fragment="SWIG_null_deleter_python") CONST CLASS &
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    new Teuchos::RCP< CONST CLASS >(&$1 SWIG_NO_NULL_DELETER_0);
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                    $descriptor(Teuchos::RCP< CLASS > *),
                                    SWIG_POINTER_OWN));
}

// plain pointer by reference
// Note: $disown not implemented by default as it will lead to a memory leak of the RCP instance
%typemap(in) CLASS *CONST& ($*1_ltype temp = 0,
                           Teuchos::RCP< CONST CLASS > tempshared)
{
  int newmem = 0;
  Teuchos::RCP< CLASS > * smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  tempshared = *smartarg;
  temp = %const_cast(tempshared.get(), $*1_ltype);
  $1 = &temp;
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(out, fragment="SWIG_null_deleter_python") CLASS *CONST&
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    new Teuchos::RCP< CONST CLASS >(*$1 SWIG_NO_NULL_DELETER_$owner);
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                 $descriptor(Teuchos::RCP< CLASS > *),
                                 SWIG_POINTER_OWN));
}

%typemap(varin) CLASS *CONST& %{
#error "varin typemap not implemented"
%}

%typemap(varout) CLASS *CONST& %{
#error "varout typemap not implemented"
%}

// RCP by value
%typemap(in) Teuchos::RCP< CONST CLASS >
{
  int newmem = 0;
  Teuchos::RCP< CLASS > * smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  $1 = Teuchos::rcp_const_cast< CONST CLASS >(*smartarg);
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(out) Teuchos::RCP< CONST CLASS > {
  Teuchos::RCP< CONST CLASS > *smartresult =
    $1 ? new Teuchos::RCP< CONST CLASS >($1) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                 $descriptor(Teuchos::RCP< CLASS > *),
                                 SWIG_POINTER_OWN));
}

%typemap(varin) Teuchos::RCP< CONST CLASS >
{
  int newmem = 0;
  Teuchos::RCP< CLASS > smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  $1 = Teuchos::rcp_const_cast< CONST CLASS >(*smartarg);
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(varout) Teuchos::RCP< CONST CLASS >
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    $1 ? new Teuchos::RCP< CONST CLASS >($1) : 0;
  %set_varoutput(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                    $descriptor(Teuchos::RCP< CLASS > *),
                                    SWIG_POINTER_OWN));
}

// RCP by reference
%typemap(in) Teuchos::RCP< CONST CLASS > & (Teuchos::RCP< CONST CLASS > tempshared)
{
  int newmem = 0;
  Teuchos::RCP< CLASS > * smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  tempshared = Teuchos::rcp_const_cast< CONST CLASS >(*smartarg);
  $1 = &tempshared;
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(out) Teuchos::RCP< CONST CLASS > &
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    *$1 ? new Teuchos::RCP< CONST CLASS >(*$1) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                 $descriptor(Teuchos::RCP< CLASS > *),
                                 SWIG_POINTER_OWN));
}

%typemap(varin) Teuchos::RCP< CONST CLASS > &
%{
#error "varin typemap not implemented"
%}

%typemap(varout) Teuchos::RCP< CONST CLASS > &
%{
#error "varout typemap not implemented"
%}

// RCP by pointer
%typemap(in) Teuchos::RCP< CONST CLASS > * (Teuchos::RCP< CONST CLASS > tempshared)
{
  int newmem = 0;
  Teuchos::RCP< CLASS > * smartarg = CONVERTER($input, &newmem);
  if (!smartarg) SWIG_fail;
  tempshared = Teuchos::rcp_const_cast< CONST CLASS >(*smartarg);
  $1 = &tempshared;
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(out) Teuchos::RCP< CONST CLASS > *
{
  Teuchos::RCP< CONST CLASS > *smartresult =
    $1 && *$1 ? new Teuchos::RCP< CONST CLASS >(*$1) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                 $descriptor(Teuchos::RCP< CLASS > *),
                                 SWIG_POINTER_OWN));
  if ($owner) delete $1;
}

%typemap(varin) Teuchos::RCP< CONST CLASS > *
%{
#error "varin typemap not implemented"
%}

%typemap(varout) Teuchos::RCP< CONST CLASS > *
%{
#error "varout typemap not implemented"
%}

// RCP by pointer reference
%typemap(in) Teuchos::RCP< CONST CLASS > *& (Teuchos::RCP< CONST CLASS > tempshared,
                                            $*1_ltype temp = 0)
{
  int newmem = 0;
  Teuchos::RCP< CLASS > * smartarg = CONVERTER($input, &newmem);
  tempshared = %reinterpret_cast(*smartarg, Teuchos::RCP< CONST CLASS >);
  temp = &tempshared;
  $1 = &temp;
  if (newmem & SWIG_CAST_NEW_MEMORY) delete smartarg;
}

%typemap(out) Teuchos::RCP< CONST CLASS > *&
{
  Teuchos::RCP< CONST CLASS > *smartresult =
     *$1 && **$1 ? new Teuchos::RCP< CONST CLASS >(**$1) : 0;
  %set_output(SWIG_NewPointerObj(%as_voidptr(smartresult),
                                 $descriptor(Teuchos::RCP< CLASS > *),
                                 SWIG_POINTER_OWN));
}

%typemap(varin) Teuchos::RCP< CONST CLASS > *&
%{
#error "varin typemap not implemented"
%}

%typemap(varout) Teuchos::RCP< CONST CLASS > *&
%{
#error "varout typemap not implemented"
%}

// Typecheck typemaps
%typemap(typecheck,precedence=SWIG_TYPECHECK_POINTER,noblock=1)
  CONST CLASS,
  CONST CLASS &,
  CONST CLASS *,
  CONST CLASS *&,
  Teuchos::RCP< CONST CLASS >,
  Teuchos::RCP< CONST CLASS > &,
  Teuchos::RCP< CONST CLASS > *,
  Teuchos::RCP< CONST CLASS > *&
{
  $1 = 0;
  int res = SWIG_ConvertPtr($input, 0, $descriptor(Teuchos::RCP< CLASS > *), 0);
  if (SWIG_CheckState(res)) $1 = 1;
  else
    if (PyObject_HasAttrString($input, "__distarray__")) $1 = 1;
    else
      if ((Teuchos::DefaultComm<int>::getComm()->getSize() == 1) &&
          (PyArray_Check($input))) $1 = 1;
}

%template() Teuchos::RCP< CONST CLASS >;

%enddef

%define %teuchos_rcp_dap(CONVERTER, CLASS...)
  %teuchos_rcp_dap_typemaps(SWIGEMPTYHACK, CONVERTER, CLASS)
  %teuchos_rcp_dap_typemaps(const        , CONVERTER, CLASS)
  %teuchos_rcp_typemaps_overrides(SWIGEMPTYHACK, CLASS)
  %teuchos_rcp_typemaps_overrides(const        , CLASS)
%enddef
