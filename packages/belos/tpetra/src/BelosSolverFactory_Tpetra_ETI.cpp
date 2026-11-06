// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Explicit, DLL-exported instantiation of Belos::Impl::SolverFactoryParent
// for every Tpetra ETI combination. SolverFactoryParent's registry
// (get_solverManagers()) is a function-local static inside a template static
// method - vague linkage - so without this file every consumer of
// BelosSolverFactory_Tpetra.hpp (this DLL's own registrar .cpp files, and any
// downstream executable or shared library) would instantiate its own private,
// unsynchronized copy of the registry across the Windows DLL boundary. This
// is the one translation unit that owns the real, exported instantiation;
// BelosSolverFactory_Tpetra.hpp declares the matching extern template for
// every other translation unit to import instead of re-instantiating.

#include "BelosSolverFactory_Tpetra.hpp"
#include "BelosTpetra_DLLExportMacro.h"
#include "TpetraCore_ETIHelperMacros.h"

TPETRA_ETI_MANGLING_TYPEDEFS()

#define BELOS_TPETRA_SOLVERFACTORYPARENT_INSTANT( SC, LO, GO, NT ) \
  template class BELOSTPETRA_LIB_DLL_EXPORT \
    Belos::Impl::SolverFactoryParent<SC, ::Tpetra::MultiVector<SC,LO,GO,NT>, \
      ::Tpetra::Operator<SC,LO,GO,NT>, ::Belos::DefaultDenseMatrix<int,SC> >;

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( BELOS_TPETRA_SOLVERFACTORYPARENT_INSTANT )

#undef BELOS_TPETRA_SOLVERFACTORYPARENT_INSTANT
