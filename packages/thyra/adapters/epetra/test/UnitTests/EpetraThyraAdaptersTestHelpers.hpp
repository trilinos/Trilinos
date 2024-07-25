// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_as.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


//
// Helper code and declarations
//

using Teuchos::as;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Ptr;
using Teuchos::outArg;
using Teuchos::Array;
using Teuchos::Comm;
typedef Teuchos_Ordinal Ordinal;


int g_localDim = 4;
bool g_dumpAll = false;
bool g_show_all_tests = false;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "local-dim", &g_localDim, "Local dimension of each vector." );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "show-all-tests", "no-show-all-tests", &g_show_all_tests,
    "Set if all tests are shown or not." );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "dump-all", "no-dump-all", &g_dumpAll,
    "Dump lots of data" );
}


RCP<const Epetra_Comm> getEpetraComm()
{
#ifdef HAVE_MPI
  return rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  return rcp(new Epetra_SerialComm());
#endif
}



} // namespace

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

