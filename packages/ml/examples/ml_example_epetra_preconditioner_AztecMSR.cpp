
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

// Example of wrappers for Aztec MSR or VBR matrices to Epetra_RowMatrix.
//
// The code creates a distributed MSR matrix (in the same way of example
// ml_aztec_simple.cpp and ml_aztec_simple_METIS.cpp), then creates the
// Epetra objects required by the class ML_Epetra::MultiLevelPreconditioner.
//
// This example makes use of Aztec2Petra() (contained in the aztecoo package).
// This function creates shallow copies for vectors and for VBR matrices,
// while MSR matrices (as in this example) are deep copied. The user
// has to delete the created objects. For more details, please refer
// to the documentation in Trilinos/packages/aztecoo/src/Aztec2Petra.h.

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif
#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)


#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Aztec2Petra.h"
#include "Teuchos_ParameterList.hpp"

// includes required by ML
#include "ml_epetra_preconditioner.h"
#include "ml_include.h"

// data required by this example to build the Aztec matrix
struct partition_data {               
  int *my_global_ids;      /* my_global_ids[i]: id of ith local unknown.     */
  int *needed_external_ids;/* global ids of ghost unknowns.                  */
  int Nlocal;              /* Number of local unknowns.                      */
  int Nglobal;             /* Number of global unknowns.                     */
  int *my_local_ids;       /* my_local_ids[i]: local id of ith global unknown*/
  int mypid;               /* processor id                                   */
  int nprocs;              /* total number of processors.                    */
  int Nghost;              /* number of ghost variables on processor.        */
};

extern void        BuildPartition(struct partition_data *Partition);
extern AZ_MATRIX   *BuildMatrix(struct partition_data *);

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  int    Nnodes=128*128;          /* Total number of nodes in the problem.*/
                                    /* 'Nnodes' must be a perfect square.   */
  struct       partition_data Partition = {NULL, NULL,0,0,NULL,0,0,0};

  AZ_MATRIX    *AztecMatrix;
  int          proc_config[AZ_PROC_SIZE];

  // init for MPI and Aztec (old style)
#ifdef ML_MPI
  MPI_Init(&argc,&argv);
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif

  // create a linear decomposition (in Aztec format)
  Partition.Nglobal = Nnodes;
  BuildPartition(&Partition);

  // create a matrix as Aztec matrix (in this case MSR format)
  AztecMatrix = BuildMatrix( &Partition);

  // create room for solution and rhs. Random numbers will be inserted
  // later (using Epetra_Vector methods)
  double * lhs = new double[Partition.Nlocal+Partition.Nghost];
  double * rhs = new double[Partition.Nlocal];

  // =========================================== //
  // conversion from MSR/VBR to Epetra_RowMatrix //
  // =========================================== //
  
  Epetra_Comm * Comm;                      // communicator for Epetra ( == proc_config)
  Epetra_RowMatrix * EpetraMatrix;         // matrix ( == AztecMatrix)
  Epetra_Vector * EpetraLHS, * EpetraRHS;  // solution and rhs
  Epetra_BlockMap * Map;                   // row decomposition ( == update) 
  int * global_indices = 0;                // used for VBR
  
  Aztec2Petra(proc_config, AztecMatrix, lhs, rhs, Comm, Map,
	      EpetraMatrix, EpetraLHS, EpetraRHS, &global_indices);

  // at this point AztecMatrix is no longer required if it is an MSR
  // matrix (deep copied into EpetraMatrix); AztecMatrix is still reqruired
  // for VBR matrices (only shallow copied into EpetraMatrix).
  
  EpetraLHS->Random();
  EpetraRHS->Random();

  // ======================== //
  // build the linear problem //
  // ======================== //

  // build the epetra linear problem
  Epetra_LinearProblem Problem(EpetraMatrix, EpetraLHS, EpetraRHS);

  // build the AztecOO linear problem (used by AztecOO)
  AztecOO solver(Problem);

  // ==================================================== //
  // build MultiLevelPreconditioner based on EpetraMatrix //
  // ==================================================== //

  // create a parameter list for ML options
  Teuchos::ParameterList MLList;

  // Set defaults for classic smoothed aggregation
  // Refer to the user's guide for the available parameters,
  // and also to the example files ml_example_epetra_preconditioner.cpp
  // and ml_example_epetra_preconditioner_2level.cpp for more
  // examples.
  ML_Epetra::SetDefaults("DD",MLList);

  // solve with symmetric Gauss-Seidel
  MLList.set("smoother: type", "symmetric Gauss-Seidel");

  // create the preconditioner object and compute hierarchy
  // (see comments contained in file ml_example_epetra_preconditioner.cpp)
  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*EpetraMatrix, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // ================== //
  // solve with AztecOO //
  // ================== //

  // GMRES is required if the preconditioner is non-symmetric (this
  // depends on user's options). AZ_cg can be used for symmetric
  // preconditioners

  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-5 tolerance  
  solver.Iterate(500, 1e-5);

  // =========== //
  // free memory //
  // =========== //
  
  // these are the objects of this example
  delete [] lhs;
  delete [] rhs;
  if (AztecMatrix != NULL) {
    AZ_free(AztecMatrix->bindx);
    AZ_free(AztecMatrix->val);
    AZ_free(AztecMatrix->data_org);
    AZ_matrix_destroy(&AztecMatrix);
  }
    
  // these are the objects created by Aztec2Petra
  delete EpetraMatrix;
  delete EpetraLHS;
  delete EpetraRHS;
  delete Map;
  if (global_indices!=0) AZ_free((void *) global_indices);  
  delete MLPrec;
  
#ifdef ML_MPI
  MPI_Finalize();
#endif
  return 0;
}

/* Assign unknowns to processors */
void BuildPartition(struct partition_data *Partition)
{
  int    proc_config[AZ_PROC_SIZE];

#ifdef ML_MPI
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif
  
  AZ_input_update(NULL,&(Partition->Nlocal), &(Partition->my_global_ids),
		  proc_config, Partition->Nglobal, 1, AZ_linear);
  Partition->Nghost = 0;   /* will be computed later */
}

/********************************************************************/
/* Set up Kn_mat                                                    */
/*  1) First set it up as an Aztec DMSR matrix                      */
/*     This stores the diagonal of row j in Ke_val[j] and stores a  */
/*     pointer to row j's off-diagonal nonzeros in Ke_bindx[j] with */
/*     column/nonzero entries in Ke_bindx[Ke_bindx[j]:Ke_bindx[j+1]-1] */
/*     and Ke_val[Ke_bindx[j]:Ke_bindx[j+1]-1].                     */
/*  2) call AZ_transform to convert global column indices to local  */
/*     indices and to set up Aztec's communication structure.       */
/*  3) Stuff the arrays into an Aztec matrix.                       */
/*------------------------------------------------------------------*/

AZ_MATRIX *BuildMatrix(struct partition_data *Partition)
{
  int *Kn_bindx;
  double *Kn_val;
  int    proc_config[AZ_PROC_SIZE];
  AZ_MATRIX *Kn_mat;
  int    *reordered_glob = NULL, *cpntr = NULL, *Kn_data_org = NULL;
  int i, ii, jj, nx, gid, Nlocal, nz_ptr;
  int *reordered_externs = NULL;  /* Aztec thing */
  size_t size;

  Nlocal = Partition->Nlocal;
  Kn_bindx = (int    *) malloc((27*Nlocal+5)*sizeof(int));
  Kn_val   = (double *) malloc((27*Nlocal+5)*sizeof(double));
  Kn_bindx[0] = Nlocal+1;

  nx = (int) sqrt( ((double) Partition->Nglobal) + .00001);

  for (i = 0; i < Nlocal; i++) {
    gid = (Partition->my_global_ids)[i];

    nz_ptr = Kn_bindx[i];
    ii = gid%nx;
    jj = (gid - ii)/nx;


    if (ii != nx-1) { Kn_bindx[nz_ptr] = gid+ 1; Kn_val[nz_ptr++] = -1.;}
    if (jj != nx-1) { Kn_bindx[nz_ptr] = gid+nx; Kn_val[nz_ptr++] = -1.;}
    if (jj !=    0) { Kn_bindx[nz_ptr] = gid-nx; Kn_val[nz_ptr++] = -1.;}
    if (ii !=    0) { Kn_bindx[nz_ptr] = gid- 1; Kn_val[nz_ptr++] = -1.;}

    Kn_val[i] = 4.;
    Kn_bindx[i+1] = nz_ptr;
  }

  /* Transform the global Aztec matrix into a local Aztec matrix. That is,   */
  /* replace global column indices by local indices and set up communication */
  /* data structure 'Ke_data_org' that will be used for matvec's.            */

#ifdef ML_MPI
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif

  AZ_transform(proc_config,&(Partition->needed_external_ids),
			       Kn_bindx, Kn_val, Partition->my_global_ids,
			       &reordered_glob, &reordered_externs, 
			       &Kn_data_org, Nlocal, 0, 0, 0, 
			       &cpntr, AZ_MSR_MATRIX);
  Partition->Nghost = Kn_data_org[AZ_N_external];

  /* Convert old style Aztec matrix to newer style Aztec matrix */

  Kn_mat = AZ_matrix_create( Nlocal );
  AZ_set_MSR(Kn_mat, Kn_bindx, Kn_val, Kn_data_org, 0, NULL, AZ_LOCAL);

  size = sizeof(double)*(Nlocal+Partition->Nghost);

  AZ_free(reordered_glob);
  AZ_free(reordered_externs);
  
  return(Kn_mat);
} 

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --enable-epetra --enable-teuchos\nto run this example");

  return 0;
}

#endif
