// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DLLEXPORTMACRO_H
#define TPETRA_DLLEXPORTMACRO_H

// Variables lack import thunks; without dllimport the compiler emits a direct
// reference the linker cannot resolve from an import library (LNK2019).
//
// CMake injects -Dtpetra_EXPORTS only when compiling tpetra itself, and
// propagates -DTPETRA_SHARED_BUILD to consumers via INTERFACE compile
// definitions for shared builds only.  Static builds see neither, so the
// macro expands to empty.  BUILD_SHARED_LIBS is a CMake variable and is
// NOT propagated as a -D compiler flag.

#ifdef _WIN32
#ifdef tpetra_EXPORTS
#define TPETRACORE_LIB_DLL_EXPORT __declspec(dllexport)
#elif defined(TPETRA_SHARED_BUILD)
#define TPETRACORE_LIB_DLL_EXPORT __declspec(dllimport)
#else
#define TPETRACORE_LIB_DLL_EXPORT
#endif
#else
#define TPETRACORE_LIB_DLL_EXPORT
#endif

#endif  // TPETRA_DLLEXPORTMACRO_H
