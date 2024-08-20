// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_EPETRA_UTILS_H
#define BELOS_EPETRA_UTILS_H

/*! \file BelosEpetraUtils.h
    \brief Provides utilities for Epetra examples and tests.
*/

#include "Epetra_config.h"
#include "Belos_config.h"

// RWH TODO FIXME: Trilinos_Util_distrib_msr_matrix available only if 32 bit GIDs available.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

#ifdef EPETRA_MPI
#include "mpi.h"
#endif

class Epetra_Map;
class Epetra_Comm;
class Epetra_CrsMatrix;
class Epetra_MultiVector;

#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
using Teuchos::RCP;
using Teuchos::rcp;

namespace Belos {

  namespace Util {

#ifdef HAVE_BELOS_TRIUTILS

  int createEpetraProblem(
      std::string              &filename
      ,RCP<Epetra_Map>         *rowMap
      ,RCP<Epetra_CrsMatrix>   *A
      ,RCP<Epetra_MultiVector> *B
      ,RCP<Epetra_MultiVector> *X
      ,int                     *MyPID
      ,int                     &numRHS
      );

  int createEpetraProblem(
      std::string              &filename
      ,RCP<Epetra_Map>         *rowMap
      ,RCP<Epetra_CrsMatrix>   *A
      ,RCP<Epetra_MultiVector> *B
      ,RCP<Epetra_MultiVector> *X
      ,int                     *MyPID
      );

#endif

  int rebalanceEpetraProblem(
      RCP<Epetra_Map>         &rowMap
      ,RCP<Epetra_CrsMatrix>   &A
      ,RCP<Epetra_MultiVector> &B
      ,RCP<Epetra_MultiVector> &X
      ,Epetra_Comm             &Comm
      );

  } // namespace Util

  namespace Test {
    class MPISession {
    public:
      MPISession (Teuchos::Ptr<int> argc, Teuchos::Ptr<char**> argv);
      ~MPISession ();
      Teuchos::RCP<const Epetra_Comm> getComm ();

    private:
      // Lazily initialized Epetra communicator wrapper.
      Teuchos::RCP<Epetra_Comm> comm_;
    };
  } // namespace Test

} // namespace Belos

#endif // EPETRA_NO_32BIT_GLOBAL_INDICES

#endif // BELOS_EPETRA_UTILS_H
