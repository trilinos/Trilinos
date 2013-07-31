/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#include "Ifpack.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "Epetra_MultiVector.h"

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_Array.hpp"
#include <string>
#include <stdio.h>
#include <map>

using Teuchos::RCP;
using Teuchos::rcp;

const int N = 10;
const int MatType = 3; //0 -> Unit diagonal, 1 -> Dense, val=col, 2 -> Random Dense, 3 -> Random Sparse
const double tol = 1E-6;

TEUCHOS_UNIT_TEST( Ifpack_Hypre, Construct ) {

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Map RowMap(N, 0, comm);
  Epetra_CrsMatrix Matrix(Copy, RowMap, 1);
  for(int i = 0; i < N; i++){
    int indices[1];
    double values[1];
    indices[0] = i;
    values[0] = 1.0;
    Matrix.InsertGlobalValues(i, 1, values, indices);
  }
  Matrix.FillComplete(); 
  //Ifpack_Hypre preconditioner(&Matrix);
  //preconditioner.Initialize();
}

TEUCHOS_UNIT_TEST( Ifpack_Hypre, ParameterList ){
  
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  Epetra_Map*            Map;
  // pointer to the matrix to be created
  Epetra_CrsMatrix*      Matrix;
  // container for parameters
  Teuchos::ParameterList GaleriList;
  // here we specify the global dimension of the problem
  int nx = 10 * Comm.NumProc();
  int ny = 10 * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);

  try
  {
    // Creates a simple linear map; for more details on the map creation
    // refer to the documentation
    Map = Galeri::CreateMap("Cartesian2D", Comm, GaleriList);

    // Creates a diagonal matrix with 1's on the diagonal
    Matrix   = Galeri::CreateCrsMatrix("Biharmonic2D", Map, GaleriList);

    // To created objects must be free'd using delete
    Ifpack_Hypre preconditioner(Matrix);
  
    
    int numVec = 2;
    Epetra_MultiVector X(Matrix->RowMatrixRowMap(), 2);
    Epetra_MultiVector KnownX(Matrix->RowMatrixRowMap(), 2);
    KnownX.Random();
    Epetra_MultiVector B(Matrix->RowMatrixRowMap(), 2);
    Matrix->Multiply(false, KnownX, B);

    //AztecOO problem(&matrix, &X, &B);
    //problem.SetPrecOperator(&matrix);
    //problem.Iterate(1000, 1e-7);
  
    Teuchos::ParameterList list("Preconditioner List");
  //RCP<FunctionParameter> functs[11];
  //functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 1000)); /* max iterations */
  //functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, 1e-7)); /* conv. tolerance */
  //functs[2] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTwoNorm, 1)); /* use the two norm as the stopping criteria */
  //functs[3] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetPrintLevel, 2)); /* print solve info */
  //functs[4] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetLogging, 1)); /* needed to get run info later */
  //functs[5] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, 1)); /* print amg solution info */
  //functs[6] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, 6));
  //functs[7] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, 6)); /* Sym G.S./Jacobi hybrid */ 
  //functs[8] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, 1));
  //functs[9] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, 0.0)); /* conv. tolerance zero */
  //functs[10] = rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, 1)); /* do only one iteration! */

  //list.set("Solver", PCG);
  //list.set("Preconditioner", BoomerAMG);
  //list.set("SolveOrPrecondition", Solver);
  //list.set("SetPreconditioner", true);
  //list.set("NumFunctions", 11);
  //list.set<RCP<FunctionParameter>*>("Functions", functs);

    preconditioner.SetParameters(list);
    preconditioner.Compute();
    //delete preconditioner;
    delete Map;
    delete Matrix;
  }
  catch (Galeri::Exception& rhs)
  {
    if (Comm.MyPID() == 0)
    {
      cerr << "Caught exception: ";
      rhs.Print();
    }
  }

}

