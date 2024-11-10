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

#ifndef PYTRILINOS_TEUCHOS_HEADERS_
#define PYTRILINOS_TEUCHOS_HEADERS_

// Configuration
#include "Python3Compat.hpp"
#include "PyTrilinos_config.h"
#include "PyTrilinos_PythonException.hpp"
#include "Teuchos_config.h"
#include "Teuchos_DLLExportMacro.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif
#ifdef HAVE_STDINT_H
#undef HAVE_STDINT_H
#endif

// Teuchos include files
#include "PyTrilinos_Teuchos_Util.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_DefaultComm.hpp"
#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_FILEstream.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_LabeledObject.hpp"
#include "Teuchos_NullIteratorTraits.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_ParameterListExceptions.hpp"
#include "Teuchos_ParameterListModifier.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ReductionOp.hpp"
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerializationTraits.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_XMLInputSource.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLObjectImplem.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_any.hpp"
#include "Teuchos_iostream_helpers.hpp"

#endif

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

