/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
*/

#ifdef EPETRA_MPI
#include <mpi.h>
#endif
#include "Epetra_C_wrappers.h"

int main(int argc, char *argv[]) {
  int i;
  EPETRA_OBJECT_PTR Comm, Map, X, Y;
  int MyPID, NumProc;
  int NumGlobalElements;
  int NumMyElements;

#ifdef EPETRA_MPI
  /* Initialize MPI */
  MPI_Init(&argc,&argv);
  Comm = epetra_mpicomm_create2( MPI_COMM_WORLD );
#else
  Comm = epetra_serialcomm_create();
#endif

   MyPID = epetra_comm_mypid(Comm);
   NumProc = epetra_comm_numproc(Comm);


   /* Construct a Map that puts 2 elements on each PE */

   NumGlobalElements = 2*NumProc;
   Map = epetra_map_create1(NumGlobalElements, 0, Comm);
  
 
  X = epetra_vector_create1(Map);
  Y = epetra_vector_create1(Map);

  epetra_vector_random(X);
  epetra_vector_random(Y);
  printf("Contents of X vector\n");
  epetra_vector_print(X);


  printf("Contents of Y vector\n");
  epetra_vector_print(Y);

  /* Add X and Y (need to pass Y twice for now, since this is the only update 
     interface wrapped by C at this time) */
  epetra_vector_update(X, 1.0, Y, 0.0, Y, 1.0);

  printf("Sum of X and Y vectors\n");
  epetra_vector_print(X);

  epetra_vector_destroy(X);
  epetra_vector_destroy(Y);
  epetra_map_destroy(Map);
  epetra_comm_destroy(Comm);
  
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return 0 ;
}
