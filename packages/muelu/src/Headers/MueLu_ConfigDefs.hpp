// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_CONFIGDEFS_HPP
#define MUELU_CONFIGDEFS_HPP

#include "MueLu_config.hpp"

#include <Teuchos_ConfigDefs.hpp>

// Tpetra
#include <Tpetra_KokkosCompat_DefaultNode.hpp>  // default template parameter of many MueLu classes

// Memory management
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_RCP.hpp>

// Verbose levels
#include <Teuchos_Describable.hpp>

// Misc
#include <Teuchos_ParameterList.hpp>

// Default S,L,G,N types, for use as default template arguments
// Available as MueLu::DefaultScalar, MueLu::DefaultLocalOrdinal, etc.
#include <MueLu_Details_DefaultTypes.hpp>

// Special macro for exception testing
// MUELU_TEST_FOR_EXCEPTION is only active if MueLu is configured with MueLu_ENABLE_DEBUG:BOOL=ON
// If you want an exception test both in the release and debug version of MueLu you still can use directly
// TEUCHOS_TEST_FOR_EXCEPTION
#ifdef HAVE_MUELU_DEBUG
#define MUELU_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg) \
  TEUCHOS_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg);
#else
#define MUELU_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg)
#endif

//! Namespace for MueLu classes and methods
namespace MueLu {

// import Teuchos memory management classes into MueLu
using Teuchos::arcp;
using Teuchos::arcp_reinterpret_cast;
using Teuchos::arcpFromArrayView;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::null;
using Teuchos::ParameterList;
using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_implicit_cast;
using Teuchos::rcp_static_cast;
using Teuchos::rcpFromRef;

// verbose levels
using Teuchos::VERB_DEFAULT;
using Teuchos::VERB_EXTREME;
using Teuchos::VERB_HIGH;
using Teuchos::VERB_LOW;
using Teuchos::VERB_MEDIUM;
using Teuchos::VERB_NONE;

}  // namespace MueLu

// This include file defines macros to avoid warnings under CUDA.  See github issue #1133.
#include "Teuchos_CompilerCodeTweakMacros.hpp"

#endif /* MUELU_CONFIGDEFS_H */
