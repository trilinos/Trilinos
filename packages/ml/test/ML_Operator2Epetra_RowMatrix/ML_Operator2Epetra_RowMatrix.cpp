/* 
 * Goal of this test:
 * - compare the two converters from ML_Operator to Epetra_RowMatrix.
 *   ML offers two ways to convert from ML_Operator to an Epetra_RowMatrix
 *   derived class. The first way is ML_Operator2EpetraCrsMatrix, which
 *   creates a new Epetra_CrsMatrix, and fills it with the elements of the
 *   input ML_Operator. This is an expensive conversion.
 *   The other conversion is given by class ML_Epetra::RowMatrix, that defines
 *   suitable wraps from the ML data format, to the Epetra data format.
 *
 * This test will:
 * - create an Aztec matrix;
 * - convert this matrix into ML_Operator;
 * - create an Epetra_CrsMatrix;
 * - create an ML_Epetra::RowMatrix;
 * - multiply those two matrices by a random vector, and compare
 *   the results.
 *
 * \date 29-Aug-04
 *
 * \author Marzio Sala, SNL 9214
 *
 */

#include <iostream>
#include <math.h>
#include "az_aztec.h"
#include "ml_include.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "ml_epetra_utils.h"
#include "ml_RowMatrix.h"

// The matrix construction part is copied from ml_aztec_simple.c.
struct user_partition_data {               
  int *my_global_ids;      /* my_global_ids[i]: id of ith local unknown.     */
  int *needed_external_ids;/* global ids of ghost unknowns.                  */
  int Nlocal;              /* Number of local unknowns.                      */
  int Nglobal;             /* Number of global unknowns.                     */
  int *my_local_ids;       /* my_local_ids[i]: local id of ith global unknown*/
  int mypid;               /* processor id                                   */
  int nprocs;              /* total number of processors.                    */
  int Nghost;              /* number of ghost variables on processor.        */
};

extern void        user_partition(struct user_partition_data *Partition);
extern AZ_MATRIX   *user_Kn_build(struct user_partition_data *);

#ifdef ML_MPI
#define COMMUNICATOR   MPI_COMM_WORLD
#else
#define COMMUNICATOR   AZ_NOT_MPI
#endif


int main(int argc, char *argv[])
{
  int ierr;
  int NumVectors = 16;
  int NumGlobalNodes = 1024;
  struct user_partition_data Partition = {NULL, NULL,0,0,NULL,0,0,0};

  /* See Aztec User's Guide for information on these variables */

  AZ_MATRIX    *Kn_mat;
  int          proc_config[AZ_PROC_SIZE];

  /* get processor information (id & # of procs) and set ML's printlevel. */

#ifdef ML_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  Partition.Nglobal = NumGlobalNodes;
  user_partition(&Partition);

  // FIXME: this is just to create the ML_Comm object
  ML* ml_handle;
  ML_Create(&ml_handle,1);

  // create an AZ_MATRIX, convert it to ML_Operator
  ML_Operator* ML_Mat;
  AZ_MATRIX* AZ_Mat;

  AZ_Mat = user_Kn_build( &Partition);
  ML_Mat = ML_Operator_Create(ml_handle->comm);
  AZ_convert_aztec_matrix_2ml_matrix(AZ_Mat, ML_Mat, proc_config);

  // 1.- convert ML_Operator to a newly created Epetra_CrsMatrix.
  //     This conversion is expensive in that a new matrix
  //     is stored, that contains the elements of the ML_Operator.
  //     This ML_Operator is no longer required by the Epetra_CrsMatrix

  Epetra_CrsMatrix* CrsMatrix;

  int MaxNumNonzeros;
  double CPUTime;

  ML_Operator2EpetraCrsMatrix(ML_Mat, CrsMatrix, MaxNumNonzeros,
			      true, CPUTime);

  // 2.- convert ML_Operator to ML_Epetra::RowMatrix
  ML_Epetra::RowMatrix* MLMatrix = new ML_Epetra::RowMatrix(ML_Mat,Comm);

  // 3.- be sure that the two matrices have the same row map

  if (!CrsMatrix->RowMatrixRowMap().SameAs(MLMatrix->RowMatrixRowMap())) {
    ML_EXIT(-1); // not the same map
  }

  if (!CrsMatrix->RowMatrixColMap().SameAs(MLMatrix->RowMatrixColMap())) {
    ML_EXIT(-1); // not the same map
  }

  // 3.- be sure that the two approaches give the same results
  //     if applied to the same vector (in this case a random one)

  const Epetra_Map& RowMap = CrsMatrix->RowMatrixRowMap();

  Epetra_MultiVector X1(RowMap,NumVectors);

  X1.Random();
  Epetra_MultiVector X2(X1);
  
  Epetra_MultiVector Y1(RowMap,NumVectors);
  Epetra_MultiVector Y2(RowMap,NumVectors);

  CrsMatrix->Apply(X1,Y1);
  
  MLMatrix->Apply(X2,Y2);
  
  Y1.Update(-1.0,Y2,1.0);
  
  double* Norm2 = new double[NumVectors];
  Y1.Norm2(Norm2);

  double TotalNorm = 0.0;
  for (int i = 0 ; i < NumVectors ; ++i) {
    TotalNorm += Norm2[i];
    if (Norm2[i] > 1e-15) {
      ML_EXIT(-2); 
    }
  }
    
  // at this point the test is passed. Some fancy (??) output,
  // and I give up.
  if (Comm.MyPID() == 0)
    cout << "Total norm = " << TotalNorm << endl;

  // free memory

  ML_Destroy(&ml_handle);
  ML_Operator_Destroy(&ML_Mat);
  
  delete CrsMatrix;
  delete MLMatrix;

  if (Partition.my_local_ids != NULL) free(Partition.my_local_ids);
  if (Partition.my_global_ids != NULL) free(Partition.my_global_ids);
  if (Partition.needed_external_ids != NULL)
    free(Partition.needed_external_ids);

  if (AZ_Mat != NULL) {
    AZ_free(AZ_Mat->bindx);
    AZ_free(AZ_Mat->val);
    AZ_free(AZ_Mat->data_org);
    AZ_matrix_destroy(&AZ_Mat);
  }

  delete [] Norm2;
  
#ifdef ML_MPI
  MPI_Finalize();
#endif
  return(0);

} // main driver 

/*------------------------------------------------------------------*/

void user_partition(struct user_partition_data *Partition)
{
  int    proc_config[AZ_PROC_SIZE];

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_input_update(NULL,&(Partition->Nlocal), &(Partition->my_global_ids),
		  proc_config, Partition->Nglobal, 1, AZ_linear);
  Partition->Nghost = 0;   /* will be computed later */
}

/*------------------------------------------------------------------*/

AZ_MATRIX *user_Kn_build(struct user_partition_data *Partition)

{
  int *Kn_bindx;
  double *Kn_val;
  int    proc_config[AZ_PROC_SIZE];
  AZ_MATRIX *Kn_mat;
  int    *reordered_glob = NULL, *cpntr = NULL, *Kn_data_org = NULL;
  int i, ii, jj, nx, gid, Nlocal, nz_ptr;
  int *reordered_externs = NULL;  /* Aztec thing */


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

  AZ_set_proc_config(proc_config, COMMUNICATOR);

  AZ_transform(proc_config,&(Partition->needed_external_ids),
			       Kn_bindx, Kn_val, Partition->my_global_ids,
			       &reordered_glob, &reordered_externs, 
			       &Kn_data_org, Nlocal, 0, 0, 0, 
			       &cpntr, AZ_MSR_MATRIX);
  Partition->Nghost = Kn_data_org[AZ_N_external];
  AZ_free(reordered_glob);
  AZ_free(reordered_externs);

  /* Convert old style Aztec matrix to newer style Aztec matrix */

  Kn_mat = AZ_matrix_create( Nlocal );
  AZ_set_MSR(Kn_mat, Kn_bindx, Kn_val, Kn_data_org, 0, NULL, AZ_LOCAL);

  return(Kn_mat);
}
