// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Amesos2_config.h"

#if defined(HAVE_AMESOS2_EXPLICIT_INSTANTIATION)

// We need both the _decl.hpp and _def.hpp files here, because if ETI
// is ON, then the .hpp file will only include the _decl.hpp file.
#include "Amesos2_Details_LinearSolverFactory_decl.hpp"
#include "Amesos2_Details_LinearSolverFactory_def.hpp"
// We need this whether or not ETI is on, in order to define typedefs
// for making Tpetra's macros work.
#include "TpetraCore_ETIHelperMacros.h"

#ifdef HAVE_AMESOS2_EPETRA
// Do explicit instantiation of Amesos2::Details::LinearSolverFactory,
// for Epetra objects.
template class Amesos2::Details::LinearSolverFactory<Epetra_MultiVector, Epetra_Operator, double>;
#endif // HAVE_AMESOS2_EPETRA

// Define typedefs that make the Tpetra macros work.
TPETRA_ETI_MANGLING_TYPEDEFS()

// Do explicit instantiation of Amesos2::Details::LinearSolverFactory, for
// Tpetra objects, for all combinations of Tpetra template parameters
// for which Tpetra does explicit template instantiation (ETI).
//
// NOTE (mfh 23 Jul 2015): Amesos2 has a required dependency on
// Tpetra, so we don't have to protect use of Tpetra with a macro.
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( AMESOS2_DETAILS_LINEARSOLVERFACTORY_INSTANT )

#endif // HAVE_AMESOS2_EXPLICIT_INSTANTIATION
