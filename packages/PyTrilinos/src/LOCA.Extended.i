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

%define %loca_extended_docstring
"
PyTrilinos.LOCA.Extended is the python interface to namespace Extended
of the Trilinos continuation algorithm package LOCA:

    http://trilinos.sandia.gov/packages/nox

The purpose of LOCA.Extended is to provide classes that extend
NOX.Abstract classes to handle an arbitrary number of multi-vectors
and scalars.  The python version of LOCA.Extended supports the
following classes:

    * MultiVector         - Implemenatation of the NOX.Abstract.MultiVector
                            class for extended multi-vectors comprised of an
                            arbitrary number of multi-vectors and scalars
    * Vector              - Implemenatation of the NOX.Abstract.Vector class
                            for extended multi-vectors comprised of an
                            arbitrary number of multi-vectors and scalars
    * MultiAbstractGroup  - LOCA abstract interface for extended groups,
                            derived from the NOX.Abstract.Group, i.e., an
                            abstract interface for 'super' groups that have an
                            underlying group component
"
%enddef

%module(package   = "PyTrilinos.LOCA",
        directors = "1",
        docstring = %loca_extended_docstring) Extended

%{
// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "PyTrilinos_Teuchos_Util.hpp"

// LOCA includes
#include "LOCA.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
%}

// Configuration and optional includes
%include "PyTrilinos_config.h"
#ifdef HAVE_NOX_EPETRA
%{
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#include "Epetra_NumPyVector.hpp"
%}
#endif

// Standard exception handling
%include "exception.i"

// Include LOCA documentation
%feature("autodoc", "1");
%include "LOCA_dox.i"

// Ignore/renames
%ignore *::operator=;

// Trilinos module imports
%import "Teuchos.i"
%import "NOX.Abstract.i"

// Teuchos::RCP support
%teuchos_rcp(LOCA::Extended::MultiAbstractGroup)

//////////////////////////////////////
// LOCA::Extended MultiVector class //
//////////////////////////////////////
//%feature("director") LOCA::Extended::MultiVector;
%include "LOCA_Extended_MultiVector.H"

/////////////////////////////////
// LOCA::Extended Vector class //
/////////////////////////////////
//%feature("director") LOCA::Extended::Vector;
%ignore LOCA::Extended::Vector::getScalar(int);
%include "LOCA_Extended_Vector.H"

/////////////////////////////////////////////
// LOCA::Extended MultiAbstractGroup class //
/////////////////////////////////////////////
//%feature("director") LOCA::Extended::MultiAbstractGroup;
%include "LOCA_Extended_MultiAbstractGroup.H"
