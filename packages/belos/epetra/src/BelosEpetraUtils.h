//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef BELOS_EPETRA_UTILS_H
#define BELOS_EPETRA_UTILS_H

/*! \file BelosEpetraUtils.h
    \brief Provides utilities for Epetra examples and tests.
*/

#include "Epetra_config.h"

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

#endif // BELOS_EPETRA_UTILS_H
