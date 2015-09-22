//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Flops.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Version.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "EpetraExt_SubCopy_CrsMatrix.h"
using namespace Trilinos_Util;

// prototypes


int main(int argc, char *argv[])
{
  int ierr = 0, i, forierr = 0;
  bool debug = false;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

// Generate Laplace matrix with 4 rows on each processor
 CrsMatrixGallery laplace_2d("laplace_2d", Comm);
 laplace_2d.Set("problem_size",Comm.NumProc()*Comm.NumProc()*4);
 Epetra_CrsMatrix * laplace_2d_matrix = laplace_2d.GetMatrix();
 if (verbose) cout << "Orig matrix = " << *laplace_2d_matrix << endl;
 const Epetra_Map * origRowMap = laplace_2d.GetMap();
 int * origGids = origRowMap->MyGlobalElements();
 
 // Next generate a sub map of the original problem, taking every other GID.
 Epetra_IntSerialDenseVector newRowMapGids(origRowMap->NumMyElements()/2 + 1);
 int numNewRowGids = 0;
 for (int i=0; i<origRowMap->NumMyElements(); i=i+2)
   newRowMapGids[numNewRowGids++] = origGids[i];
 Epetra_Map newRowMap(-1, numNewRowGids, newRowMapGids.Values(), 0, origRowMap->Comm());


 // Use this new map to create a SubCopy transform that we can use on-demand to create submatrices
 EpetraExt::CrsMatrix_SubCopy subMatrixTransform(newRowMap);

 // Use the transform we just created to create our submatrix
 Epetra_CrsMatrix & subA = subMatrixTransform(*laplace_2d_matrix);
 if (verbose) cout << "Sub matrix (every other row/column) = " << subA << endl;


 // Now test it using the fwd() method to see if changes in the original matrix are carried to the submatrix
 (*laplace_2d_matrix)[0][0] = 12.0;
 if (verbose) cout << "Orig matrix = " << *laplace_2d_matrix << endl;
 subMatrixTransform.fwd();
 assert(subA[0][0]==12.0);
  if (verbose) cout << "Sub matrix (every other row/column) = " << subA << endl;


  // Now change the submatrix and use the rvs() method to see if changes are carried to the original matrix
 subA[0][0] = 24.0;
 if (verbose) cout << "Sub matrix (every other row/column) = " << subA << endl;
 subMatrixTransform.rvs();
 assert((*laplace_2d_matrix)[0][0]==24.0);
 if (verbose) cout << "Orig matrix = " << *laplace_2d_matrix << endl;

 if (Comm.MyPID()==0) cout << "EpetraExt::CrsMatrix_SubCopy tests passed." << endl;
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0;
}


