// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TSQR_DLLEXPORTMACRO_H
#define TSQR_DLLEXPORTMACRO_H

// Variables lack import thunks; without dllimport the compiler emits a direct
// reference the linker cannot resolve from an import library (LNK2019).
//
// CMake injects -Dkokkostsqr_EXPORTS only when compiling kokkostsqr itself,
// and propagates -DKOKKOSTSQR_SHARED_BUILD to consumers via INTERFACE compile
// definitions for shared builds only.  Static builds see neither, so the macro
// expands to empty.  BUILD_SHARED_LIBS is a CMake variable never passed as -D.

#ifdef _WIN32
#  ifdef kokkostsqr_EXPORTS
#    define KOKKOSTSQR_LIB_DLL_EXPORT __declspec(dllexport)
#  elif defined(KOKKOSTSQR_SHARED_BUILD)
#    define KOKKOSTSQR_LIB_DLL_EXPORT __declspec(dllimport)
#  else
#    define KOKKOSTSQR_LIB_DLL_EXPORT
#  endif
#else
#  define KOKKOSTSQR_LIB_DLL_EXPORT
#endif

#endif // TSQR_DLLEXPORTMACRO_H
