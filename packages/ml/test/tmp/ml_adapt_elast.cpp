//@HEADER
// ************************************************************************
//
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "AztecOO.h"
#include "Epetra_MsrMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
class Epetra_Map; 

// includes required by ML
#include "ml_MultiLevelPreconditioner.h"

#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Teuchos;
using namespace Trilinos_Util;

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Epetra_Time Time(Comm);

  // Create the linear problem using the class `Trilinos_Util::CrsMatrixGallery.'
  // Various matrix examples are supported; please refer to the
  // Trilinos tutorial for more details.
  
  // create Aztec stuff
  int    proc_config[AZ_PROC_SIZE], options[AZ_OPTIONS_SIZE];
#ifdef ML_MPI
  /* get number of processors and the name of this processor */
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
  int proc   = proc_config[AZ_node];
  int nprocs = proc_config[AZ_N_procs];
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
  int proc   = 0;
  int nprocs = 1;
#endif
  // read in the matrix size
  FILE *fp = fopen("ExampleMatrices/cantilever2D/data_matrix.txt","r");
  int leng;
  fscanf(fp,"%d",&leng);
  int num_PDE_eqns=2;
  int N_grid_pts = leng/num_PDE_eqns;

  // make a linear distribution of the matrix respecting the blocks size
  int leng1 = leng/nprocs;
  int leng2 = leng-leng1*nprocs;
  if (proc >= leng2)
  {
     leng2 += (proc*leng1);
  }
  else
  {
     leng1++;
     leng2 = proc*leng1;
  }
  int     N_update = leng1;
  int*    update  = new int[N_update+1];
  int     i;
  double *val=NULL;
  int    *bindx=NULL;
  for (i=0; i<N_update; i++) update[i] = i+leng2;
  
  // create the Epetra_CrSMatrix
  Epetra_Map*        StandardMap = new Epetra_Map(leng,N_update,update,0,Comm);
  Epetra_CrsMatrix*  A           = new Epetra_CrsMatrix(Copy,*StandardMap,1);
  
  AZ_input_msr_matrix("ExampleMatrices/cantilever2D/data_matrix.txt",
                      update, &val, &bindx, N_update, proc_config);

  
  for (i=0; i<leng; i++)
  {
    int row = update[i];
    A->SumIntoGlobalValues(row,1,&(val[i]),&row);
    A->SumIntoGlobalValues(row,bindx[i+1]-bindx[i],&(val[bindx[i]]),&(bindx[bindx[i]]));
  }
  A->TransformToLocal();
  
  // create solution and right-hand side (MultiVectors are fine as well)
  Epetra_Vector* LHS = new Epetra_Vector(A->OperatorDomainMap());
  Epetra_Vector* RHS = new Epetra_Vector(A->OperatorRangeMap());
  LHS->Random();
  RHS->Random();

  // build the epetra linear problem
  Epetra_LinearProblem Problem(A, LHS, RHS);
  
  // Construct a solver object for this problem
  AztecOO solver(Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("aggregation: damping factor", 0.0);

  // number of relaxation sweeps
  MLList.set("adaptive: max sweeps", 10);
  // number of additional null space vectors to compute
  MLList.set("adaptive: num vectors",2);

#if 1
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(dynamic_cast<Epetra_RowMatrix&>(*A), MLList, false);

  // need to allocate and fill the null space (also the
  // default one, as in this case). This vector is no longer
  // needed after a call to ComputeAdaptivePreconditioner().
  int NullSpaceSize = 2;
  vector<double> NullSpace((NullSpaceSize*A->NumMyRows()));
  for (i = 0 ; i < A->NumMyRows() ; ++i)
  {
    NullSpace[i] = 1.0;
    ++i;
    NullSpace[i] = 0.0;
  }
  for (i = A->NumMyRows() ; i < 2*A->NumMyRows() ; ++i)
  {
    NullSpace[i] = 0.0;
    ++i;
    NullSpace[i] = 1.0;
  }

  MLPrec->ComputeAdaptivePreconditioner(NullSpaceSize,&NullSpace[0]);
#else
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(dynamic_cast<Epetra_RowMatrix&>(*A), MLList);
#endif

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================
  
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 tolerance  
  solver.Iterate(1550, 1e-5);

  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
    cout << "Total Time = " << Time.ElapsedTime() << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(0);
  
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos");
  puts("--enable-aztecoo --enable-triutils");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
