/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/************************************************************************/
/*          Utilities for Trilinos/ML users                             */
/*----------------------------------------------------------------------*/
/* Authors : Mike Heroux (SNL)                                          */
/*           Jonathan Hu  (SNL)                                         */
/*           Ray Tuminaro (SNL)                                         */
/************************************************************************/


#include "ml_common.h"

#ifdef ML_WITH_EPETRA
#include "ml_epetra_utils.h"
#include "Epetra_FECrsMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


/******************************************************************************/

int Epetra_ML_matvec(void *data, int in, double *p, int out, double *ap)
{
  /* ML matvec wrapper for Epetra matrices. */

  Epetra_RowMatrix *A = (Epetra_RowMatrix *) data;
/*MS*/
  Epetra_Vector X(View, A->OperatorDomainMap(), p);
  Epetra_Vector Y(View, A->OperatorRangeMap(), ap);
/*ms*/
  
  A->Multiply(false, X, Y);
  
  return 1;
}

/******************************************************************************/

int Epetra_ML_getrow(void *data, int N_requested_rows, int requested_rows[], 
		    int allocated_space, int columns[], double values[],
		    int row_lengths[])
{
/*
 * GetRow function for matrix of type Epetra_RowMatrix.
 * Supply local matrix (without ghost node columns) for rows given by
 * requested_rows[0 ... N_requested_rows-1].  Return this information in
 * 'row_lengths, columns, values'.  If there is not enough space to complete
 * this operation, return 0. Otherwise, return 1.
 *
 * Parameters
 * ==========
 * data             On input, points to user's data containing matrix values.
 * N_requested_rows On input, number of rows for which nonzero are to be
 *                  returned.
 * requested_rows   On input, requested_rows[0...N_requested_rows-1] give the
 *                  row indices of the rows for which nonzero values are
 *                  returned.
 * row_lengths      On output, row_lengths[i] is the number of nonzeros in the
 *                  row 'requested_rows[i]'
 * columns,values   On output, columns[k] and values[k] contains the column
 *                  number and value of a matrix nonzero where all nonzeros for
 *                  requested_rows[i] appear before requested_rows[i+1]'s
 *                  nonzeros.  NOTE: Arrays are of size 'allocated_space'.
 * allocated_space  On input, indicates the space available in 'columns' and
 *                  'values' for storing nonzeros. If more space is needed,
 *                  return 0.
 */
  Epetra_RowMatrix *Abase = (Epetra_RowMatrix *) data;
  int nz_ptr = 0;
  int NumEntries;
  int MaxPerRow;

  Epetra_CrsMatrix * Acrs = dynamic_cast<Epetra_CrsMatrix *>(Abase);
  int MatrixIsCrsMatrix = (Acrs!=0); // If this pointer is non-zero,
                                  // the cast to Epetra_CrsMatrix worked
  int *Indices;
  double *Values;
  if (!MatrixIsCrsMatrix) {
    MaxPerRow = Abase->MaxNumEntries();
    Values = new double [MaxPerRow]; 
    Indices = new int [MaxPerRow]; 
  }

  for (int i = 0; i < N_requested_rows; i++)
  {
    int ierr;
    int LocalRow = requested_rows[i];
    if (MatrixIsCrsMatrix)
      ierr = Acrs->ExtractMyRowView(LocalRow, NumEntries, Values, Indices);
    else
      ierr = Abase->ExtractMyRowCopy(LocalRow, MaxPerRow, NumEntries,
                                      Values, Indices);
    if (ierr) return(0); //JJH I think this is the correct thing to return if
                         //    A->ExtractMyRowCopy returns something nonzero ..
    row_lengths[i] = NumEntries;
    if (nz_ptr + NumEntries > allocated_space) return(0);
    for (int j=0; j<NumEntries; j++) {
      columns[nz_ptr] = Indices[j];
      values[nz_ptr++] = Values[j];
    }
  }

  return(1);
}

/******************************************************************************/

int Epetra_ML_comm_wrapper(double vec[], void *data)
{
  /*
   * Update vec's ghost node via communication. Note: the length of vec is
   * given by N_local + N_ghost where Amat was created via
   *                 AZ_matrix_create(N_local);
   * and a 'getrow' function was supplied via
   *                 AZ_set_MATFREE_getrow(Amat, , , , N_ghost, );
 *
 * Parameters
 * ==========
 * vec              On input, vec contains data. On output, ghost values
 *                  are updated.
 *
 * data             On input, points to user's data containing matrix values.
 *                  and communication information.
 */


  Epetra_RowMatrix *A = (Epetra_RowMatrix *) data;

  if (A->Comm().NumProc()==1) return(1); // Nothing to do in serial mode.

//  Epetra_Vector X_target(View, A->RowMatrixImporter()->TargetMap(), vec); //ghosted
//  Epetra_Vector X_source(View, A->RowMatrixImporter()->SourceMap(), vec); //loc only

  if( A->RowMatrixImporter() != 0 ) {
    Epetra_Vector X_target(View, A->RowMatrixImporter()->TargetMap(),
			   vec); //ghosted
    Epetra_Vector X_source(View, A->RowMatrixImporter()->SourceMap(),
			   vec); //loc only
  
//  assert(X_target.Import(X_source, *(A->RowMatrixImporter()),Insert)==0);
    X_target.Import(X_source, *(A->RowMatrixImporter()), Insert);
  }
  
  return(1);
}
/******************************************************************************/

int Epetra2MLMatrix(Epetra_RowMatrix * A, ML_Operator *newMatrix)
{
  int isize, osize;

  osize = A->NumMyRows();
  isize = A->OperatorDomainMap().NumMyElements();
  //  isize = A->NumMyCols();
  int N_ghost = A->RowMatrixColMap().NumMyElements() - isize;

  if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

  ML_Operator_Set_ApplyFuncData(newMatrix, isize, osize,
                              ML_EMPTY, (void*) A, osize,
                              NULL, 0);

  ML_CommInfoOP_Generate( &(newMatrix->getrow->pre_comm), 
                        Epetra_ML_comm_wrapper, (void *) A, 
                        newMatrix->comm, isize, N_ghost);

  ML_Operator_Set_Getrow(newMatrix, ML_EXTERNAL, newMatrix->outvec_leng,
                       Epetra_ML_getrow);

  ML_Operator_Set_ApplyFunc (newMatrix, ML_EXTERNAL, Epetra_ML_matvec);

  return 0;
}

int EpetraMatrix2MLMatrix(ML *ml_handle, int level,
                         Epetra_RowMatrix * A)
{
  int isize, osize;

  osize = A->NumMyRows();
  isize = osize;
  int N_ghost = A->NumMyCols() - A->NumMyRows();

  if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

  ML_Init_Amatrix(ml_handle, level,isize, osize, (void *) A);
  ML_Set_Amatrix_Getrow(ml_handle, level, Epetra_ML_getrow,
            Epetra_ML_comm_wrapper, isize+N_ghost);

  ML_Set_Amatrix_Matvec(ml_handle,  level, Epetra_ML_matvec);

  return 1;
}
int ML_back_to_epetraCrs(ML_Operator *Mat1Mat2, ML_Operator *Result, 
			 ML_Operator *Mat1, ML_Operator *Mat2)
{
  int *global_rows;

  Epetra_RowMatrix *Mat1_epet = (Epetra_RowMatrix *) Mat1->data;
  Epetra_RowMatrix *Mat2_epet = (Epetra_RowMatrix *) Mat2->data;

  Epetra_CrsMatrix *Result_epet = new Epetra_CrsMatrix(Copy, 
				            Mat1_epet->RowMatrixRowMap(),
					    Mat2_epet->RowMatrixColMap(), 0);
  int allocated = 0, row_length;
  int *bindx = NULL;
  double *val = NULL;

  global_rows = Mat1_epet->RowMatrixRowMap().MyGlobalElements();
  for (int i = 0; i < Mat1Mat2->getrow->Nrows; i++) {
    ML_get_matrix_row(Mat1Mat2, 1, &i, &allocated, &bindx, &val,
		      &row_length, 0);

    Result_epet->InsertGlobalValues(global_rows[i],
					       row_length, val,
					       bindx);
  }
  int ierr=Result_epet->TransformToLocal(&(Mat1_epet->OperatorRangeMap()),
				    &(Mat2_epet->OperatorDomainMap()));

  if (bindx != NULL) ML_free(bindx);
  if (val != NULL) ML_free(val);
  if (ierr!=0) {
    cerr <<"Error in Epetra_VbrMatrix TransformToLocal" << ierr << endl;
    EPETRA_CHK_ERR(ierr);
  }

  Epetra2MLMatrix((Epetra_RowMatrix *) Result_epet, Result);

  return 1;
}


Epetra_CrsMatrix *Epetra_MatrixMult(Epetra_RowMatrix *B_crs, Epetra_RowMatrix *Bt_crs)
{
  ML_Comm *comm, *temp;
  Epetra_RowMatrix *result;

  ML_Comm_Create(&comm);
  ML_Operator *B_ml, *Bt_ml, *BBt_ml;
  B_ml  = ML_Operator_Create(comm);
  Bt_ml = ML_Operator_Create(comm);
  BBt_ml  = ML_Operator_Create(comm);
  Epetra2MLMatrix(B_crs, B_ml);
  Epetra2MLMatrix(Bt_crs, Bt_ml);
  ML_2matmult(B_ml, Bt_ml, BBt_ml, ML_EpetraCRS_MATRIX);

  ML_Comm_Destroy(&comm);
  global_comm = temp;

  /* Need to blow about BBt_ml but keep epetra stuff */

  result = (Epetra_RowMatrix *) BBt_ml->data;
  ML_Operator_Destroy(&B_ml);
  ML_Operator_Destroy(&Bt_ml);
  ML_Operator_Destroy(&BBt_ml);

  return dynamic_cast<Epetra_CrsMatrix*>(result);
   
}


Epetra_CrsMatrix *Epetra_MatrixAdd(Epetra_RowMatrix *B_crs, Epetra_RowMatrix *Bt_crs, double scalar)
{
  ML_Comm *comm, *temp;

  temp = global_comm;
  ML_Comm_Create(&comm);
  ML_Operator *B_ml, *Bt_ml, *BBt_ml;
  B_ml  = ML_Operator_Create(comm);
  Bt_ml = ML_Operator_Create(comm);
  BBt_ml  = ML_Operator_Create(comm);
  Epetra2MLMatrix(B_crs, B_ml);
  Epetra2MLMatrix(Bt_crs, Bt_ml);
  Epetra_CrsMatrix *BBt_crs = new Epetra_CrsMatrix(Copy,
				            B_crs->RowMatrixRowMap(),
					    B_crs->RowMatrixColMap(), 0);
  BBt_ml->data = (void *) BBt_crs;
  ML_Operator_Add(B_ml, Bt_ml, BBt_ml, ML_EpetraCRS_MATRIX, scalar);
  BBt_crs->TransformToLocal(&(B_crs->OperatorRangeMap()),
				     &(B_crs->OperatorDomainMap()));

  ML_Comm_Destroy(&comm);
  global_comm = temp;

  /* Need to blow about BBt_ml but keep epetra stuff */

  ML_Operator_Destroy(&B_ml);
  ML_Operator_Destroy(&Bt_ml);
  ML_Operator_Destroy(&BBt_ml);

  return BBt_crs;
   
}

int ML_Epetra_CRSinsert(ML_Operator *A, int row, int *cols, double *vals, int length)
{
  int *global_rows;
  Epetra_CrsMatrix *A_crs = (Epetra_CrsMatrix *) A->data;

  global_rows = A_crs->RowMatrixRowMap().MyGlobalElements();
  A_crs->InsertGlobalValues(global_rows[row],length, vals, cols);

  return 0;
}

extern "C" {

Epetra_CrsMatrix * Q  = NULL;
Epetra_FECrsMatrix *  Qt = NULL; 
  
ML_Operator * ML_BuildQ( int StartingNumElements,
			 int ReorderedNumElements,
			 int reordered_decomposition[],
			 double StartingNullSpace[],
			 double ReorderedNullSpace[],
			 int ComputeNewNullSpace,
			 double StartingBdry[], double ReorderedBdry[],
			 MPI_Comm mpi_communicator,
			 ML_Comm *ml_communicator ) 
{
  
  ML_Operator * ML_Q2;
  
#ifndef ML_MPI
  /* ********************************************************************** */
  /* ONe should not call this function with one processor only (as he has   */
  /* nothing to redistributed. I simply return (this is also checked later) */
  /* ********************************************************************** */

  return NULL;
#endif

  Epetra_MpiComm Comm( mpi_communicator );

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  Epetra_Map StartingMap(-1,StartingNumElements,0,Comm);
  Epetra_Map ReorderedMap(-1,ReorderedNumElements,0,Comm);

  Q = new Epetra_CrsMatrix(Copy,StartingMap,1);

  int * MyGlobalElements = StartingMap.MyGlobalElements();

  // fill Q
  for( int i=0 ; i<StartingNumElements ; ++i ) {
    double one = 1.0;
    int indices = reordered_decomposition[i];
    Q->InsertGlobalValues(MyGlobalElements[i], 1, &one, &indices );
  }

  assert(Q->FillComplete(ReorderedMap,StartingMap)==0);

  ML_Q2 = ML_Operator_Create( ml_communicator );

  if( ComputeNewNullSpace == 1 ) {
    Epetra_Vector xxx(View, StartingMap, StartingNullSpace);
    Epetra_Vector yyy(View, ReorderedMap, ReorderedNullSpace);
    Q->Multiply(true,xxx,yyy);
  }

  Epetra_Vector xxx2(View, StartingMap, StartingBdry);
  Epetra_Vector yyy2(View, ReorderedMap, ReorderedBdry);
  Q->Multiply(true,xxx2,yyy2);

  Epetra2MLMatrix( Q, ML_Q2);

  return ML_Q2;
  
}

void ML_DestroyQ(void) 
{

  delete Q;
  Q = NULL;

  return;
  
} /* ML_DestroyQ */

ML_Operator * ML_BuildQt( int StartingNumElements,
			  int ReorderedNumElements,
			  int reordered_decomposition[],
			  MPI_Comm mpi_communicator,
			  ML_Comm *ml_communicator ) 
{
  
  ML_Operator * ML_Qt2;
  
#ifndef ML_MPI
  /* ********************************************************************** */
  /* ONe should not call this function with one processor only (as he has   */
  /* nothing to redistributed. I simply return (this is also checked later) */
  /* ********************************************************************** */

  return NULL;
#endif

  Epetra_MpiComm Comm( mpi_communicator );

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if( NumProc == 1 ) return NULL;

  Epetra_Map StartingMap(-1,StartingNumElements,0,Comm);
  Epetra_Map ReorderedMap(-1,ReorderedNumElements,0,Comm);

  Qt = new Epetra_FECrsMatrix(Copy,ReorderedMap,1);

  int * MyGlobalElements = StartingMap.MyGlobalElements();

  // fill Q
  for( int i=0 ; i<StartingNumElements ; ++i ) {
    int row = reordered_decomposition[i];
    double one = 1.0;
    int indices = MyGlobalElements[i];
    Qt->SumIntoGlobalValues(1, &row, 1, &indices, &one );
  }

  Qt->GlobalAssemble(false);

  // Q will be applied to vectors defined on StartingMap,
  // and the output vector will be defined on ReorderdMap
  assert(Qt->FillComplete(ReorderedMap,StartingMap)==0);

  ML_Qt2 = ML_Operator_Create( ml_communicator );

  Epetra2MLMatrix( Qt, ML_Qt2);

  return ML_Qt2;
  
} /* ML_BuildQt */

void ML_DestroyQt( void ) 
{

  delete Qt;
  Qt = NULL;

  return;
  
} /* ML_DestroyQt */

} /* extern "C" */

#ifdef WKC
int Epetra_ML_matvec_WKC (void *data, int in, double *p, int out, double *ap)
{
  /* ML matvec wrapper for Epetra matrices. */

  Epetra_RowMatrix *A = (Epetra_RowMatrix *) data;
  Epetra_MultiVector &X(*(Epetra_MultiVector *)p);
  Epetra_MultiVector &Y(*(Epetra_MultiVector *)ap);

  A->Multiply(false, X, Y);

  return 1;
}
#endif



#else

  /*noop for certain compilers*/
  int ML_EPETRA_EMPTY;

#endif /*ifdef ML_WITH_EPETRA*/

