// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_config.hpp"

#if defined(HAVE_MUELU_EXPLICIT_INSTANTIATION)

// We need both the _decl.hpp and _def.hpp files here, because if ETI
// is ON, then the .hpp file will only include the _decl.hpp file.
#include "MueLu_Details_LinearSolverFactory_decl.hpp"
#include "MueLu_Details_LinearSolverFactory_def.hpp"
// We need this whether or not ETI is on, in order to define typedefs
// for making Tpetra's macros work.
#include "TpetraCore_ETIHelperMacros.h"

#ifdef HAVE_MUELU_EPETRA
// Do explicit instantiation of MueLu::Details::LinearSolverFactory,
// for Epetra objects.
template class MueLu::Details::LinearSolverFactory<Epetra_MultiVector, Epetra_Operator, double>;
#endif  // HAVE_MUELU_EPETRA

// Define typedefs that make the Tpetra macros work
TPETRA_ETI_MANGLING_TYPEDEFS()

// Do explicit instantiation of MueLu::Details::LinearSolverFactory, for
// Tpetra objects, for all combinations of Tpetra template parameters
// for which Tpetra does explicit template instantiation (ETI).
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(MUELU_DETAILS_LINEARSOLVERFACTORY_INSTANT)

// TODO amk: do we also have to do this for Xpetra?

#endif  // HAVE_MUELU_EXPLICIT_INSTANTIATION