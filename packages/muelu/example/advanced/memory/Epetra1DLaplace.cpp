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
#include <vector>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_Comm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Version.h>

#include "MueLu_MemoryProfiler.hpp"

int main(int argc, char *argv[]) {
  int ierr, i;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool success = false;
  bool verbose = true;
  try {
    // int myRank = Comm.MyPID();

    // int numGlobalElements = 10000000;
    int numGlobalElements = 100;

    Teuchos::CommandLineProcessor cmdp(false, true);
    cmdp.setOption("numGlobalElements", &numGlobalElements, "Global problem size.");
    if (cmdp.parse(argc, argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      throw -1;
    }

    Epetra_Map Map(numGlobalElements, 0, Comm);

    int NumMyElements = Map.NumMyElements();

    std::vector<int> MyGlobalElements(NumMyElements);
    Map.MyGlobalElements(&MyGlobalElements[0]);

    int NumNz = 3;
    // std::vector<int> NumNz(NumMyElements);
    // for (i=0; i<NumMyElements; i++)
    //     if (MyGlobalElements[i]==0 || MyGlobalElements[i] == numGlobalElements-1)
    //       NumNz[i] = 2;
    //     else
    //       NumNz[i] = 3;
    //  Epetra_CrsMatrix A(Copy, Map, &NumNz[0]);

    MemoryUsageStart("Epetra");
    PrintMemoryUsage("Initial memory usage", "epetra-init.heap");

    Epetra_CrsMatrix A(Copy, Map, NumNz);

    PrintMemoryUsage("Memory after CrsMatrix constructor", "epetra-after-ctor.heap");

    std::vector<double> Values(2);
    Values[0] = -1.0;
    Values[1] = -1.0;
    std::vector<int> Indices(2);
    double two = 2.0;
    int NumEntries;

    for (i = 0; i < NumMyElements; i++) {
      if (MyGlobalElements[i] == 0) {
        Indices[0] = 1;
        NumEntries = 1;
      } else if (MyGlobalElements[i] == numGlobalElements - 1) {
        Indices[0] = numGlobalElements - 2;
        NumEntries = 1;
      } else {
        Indices[0] = MyGlobalElements[i] - 1;
        Indices[1] = MyGlobalElements[i] + 1;
        NumEntries = 2;
      }

      ierr = A.InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert(ierr == 0);

      // Put in the diagonal entry
      ierr = A.InsertGlobalValues(MyGlobalElements[i], 1, &two, &MyGlobalElements[i]);
      assert(ierr == 0);
    }

    PrintMemoryUsage("Memory after InsertGlobalValues()", "epetra-after-insert.heap");

    ierr = A.FillComplete();
    assert(ierr == 0);

    PrintMemoryUsage("Memory after FillComplete()", "epetra-after-fillcomplete.heap");

    MemoryUsageStop();

    if (ierr == 0)
      success = true;
    else
      success = false;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
