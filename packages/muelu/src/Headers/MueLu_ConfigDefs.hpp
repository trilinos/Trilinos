// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
