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


// Petra_Comm Test routine

#include <math.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#ifdef PETRA_MPI
#include <mpi.h>
#endif
#include "Petra_Comm.h"
#include "Petra_Time.h"

int main(int argc, char *argv[]) {

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // I'm alive !!!



  Petra_Comm & petracomm = *new Petra_Comm( MPI_COMM_WORLD );
  int MyPID =  petracomm.MyPID();
  int NumProc =  petracomm.NumProc();
  cout << "Processor "<<MyPID<<" of " << NumProc << " is alive."<<endl;

  if (NumProc!=2) {
    cout << " This special routine only works for 2 processors " << endl;
    abort();
  }


  Petra_Comm * other_comm;
  MPI_Status status;

  unsigned int icomm = &petracomm;

  if (MyPID==1) cout << "Address of Petra_Comm object on PE 1 = " << &petracomm << endl;

  if (MyPID==0) MPI_Recv((void *) &other_comm, 1, MPI_UNSIGNED, 1, 99, MPI_COMM_WORLD, &status);
  else MPI_Send ( (void *) &icomm, 1, MPI_UNSIGNED, 0, 99, MPI_COMM_WORLD);

  if (MyPID==0) cout << "Address of other Petra_Comm object on PE 0 = " << other_comm << endl;
  
  int otherPID = other_comm->MyPID();

  if (MyPID==0) cout << "Processor "<<MyPID<<" of " << NumProc
		     << " has a neighbor processor with ID "
		     << otherPID << " of " << other_comm->NumProc() <<endl;
 
  delete &petracomm;
  MPI_Finalize();

  return 0;
}
 
