//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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

#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_MsrMatrix.h"

#ifdef EPETRA_MPI
#define AZ_MPI
#define AZTEC_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#endif
#include "Trilinos_Util.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
int power_method(bool TransA, Epetra_MsrMatrix& A,
           Epetra_Vector& q,
           Epetra_Vector& z,
           Epetra_Vector& resid,
           double * lambda, int niters, double tolerance,
           bool verbose);


int main(int argc, char *argv[])
{
  int    proc_config[AZ_PROC_SIZE];// Processor information.

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  AZ_set_proc_config(proc_config,MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
  AZ_set_proc_config(proc_config,0);
#endif

  int temp; if (comm.MyPID()==0) {cout << "Type 1 and enter to continue" <<endl; cin >> temp; } comm.Barrier();

  if(argc != 2) cerr << "Error: enter name of data file on command line" << endl; 

  cout << comm << endl;

  /* Set exact solution to NULL */

  int    *update;                  /* vector elements updated on this node. */
  int    *bindx;
  double *val;
  double *xguess, *b, *xexact, *xsolve;
  int    N_update;           /* # of block unknowns updated on this node    */
  int    numLocalEquations, numGlobalEquations, n_nonzeros;

  xexact = NULL;

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

    Trilinos_Util_read_hb(argv[1], comm.MyPID(), &numGlobalEquations, &n_nonzeros,
             &val,  &bindx, &xguess, &b, &xexact);

  Trilinos_Util_distrib_msr_matrix(comm, &numGlobalEquations, &n_nonzeros, &N_update,
		  &update, &val, &bindx, &xguess, &b, &xexact);

#ifdef DEBUG
  for (i = 0; i<N_update; i++)
    if (val[i] == 0.0 ) cout << "Zero diagonal at row " << i << endl;
#endif

  /* convert matrix to a local distributed matrix */
  int * external = 0;
  int * update_index = 0;
  int * extern_index = 0;
  int * data_org = 0;

  AZ_transform(proc_config, &external, bindx, val, update,
	       &update_index, &extern_index, &data_org,
	       N_update, 0, 0, 0, 0,
               AZ_MSR_MATRIX);

  AZ_MATRIX * Amat = AZ_matrix_create(N_update);
  AZ_set_MSR(Amat, bindx, val, data_org, N_update, update, AZ_LOCAL);

  Epetra_MsrMatrix A(proc_config, Amat); // Create Epetra_MsrMatrix

  const Epetra_BlockMap & map = A.RowMatrixRowMap();


  // Create vectors for Power method

  Epetra_Vector q(map);
  Epetra_Vector z(map);
  Epetra_Vector resid(map);
  
  // variable needed for iteration
  double lambda = 0.0;
  // int niters = 10000;
  int niters = 10000;
  double tolerance = 1.0e-5;

  ////////////////////////////////////////////////////////////////////////////////////////////////

  // Iterate

  Epetra_Flops flopcounter;
  A.SetFlopCounter(flopcounter);
  q.SetFlopCounter(A);
  z.SetFlopCounter(A);
  resid.SetFlopCounter(A);

  bool verbose = (comm.MyPID()==0);
  Epetra_Time timer(comm);
  power_method(false, A, q, z, resid, &lambda, niters, tolerance, verbose);
  double elapsed_time = timer.ElapsedTime();
  double total_flops = A.Flops() + q.Flops() + z.Flops() + resid.Flops();
  double MFLOPs = total_flops/elapsed_time/1000000.0;
  
  if (verbose) cout << "\n\nTotal MFLOPs for solve = " << MFLOPs << endl<< endl;


  free ((void *) xguess);
  free ((void *) b);
  free ((void *) xexact);
  free ((void *) val);
  free ((void *) bindx);
  free ((void *) update);

  // Must delete
  if (external!=0) AZ_free ((void *) external);
  if (update_index!=0) AZ_free ((void *) update_index);
  if (extern_index!=0) AZ_free ((void *) extern_index);
  if (data_org!=0) AZ_free ((void *) data_org);

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
int power_method(bool TransA, Epetra_MsrMatrix& A,
           Epetra_Vector& q,
           Epetra_Vector& z,
           Epetra_Vector& resid,
           double * lambda, int niters, double tolerance,
           bool verbose) {
    
  // Fill z with random Numbers
  z.Random();

  // variable needed for iteration
  double normz, residual;

  int ierr = 1;

  for (int iter = 0; iter < niters; iter++)
    {
      z.Norm2(&normz); // Compute 2-norm of z
      q.Scale(1.0/normz, z);
      A.Multiply(TransA, q, z); // Compute z = A*q
      q.Dot(z, lambda); // Approximate maximum eigenvaluE
      if (iter%100==0 || iter+1==niters)
     {
       resid.Update(1.0, z, -(*lambda), q, 0.0); // Compute A*q - lambda*q
       resid.Norm2(&residual);
       if (verbose) cout << "Iter = " << iter << "  Lambda = " << *lambda 
                    << "  Residual of A*q - lambda*q = " << residual << endl;
     } 
      if (residual < tolerance) {
     ierr = 0;
     break;
      }
    }
  return(ierr);
}
