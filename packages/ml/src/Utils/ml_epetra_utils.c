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

#ifdef ML_WITH_EPETRA

#include "ml_epetra_utils.h"

/******************************************************************************/

int Epetra_ML_matvec(void *data, int in, double *p, int out, double *ap)
{
  /* ML matvec wrapper for Epetra matrices. */

  Epetra_CrsMatrix *A = (Epetra_CrsMatrix *) data;
  Epetra_Vector X(View, A->DomainMap(), p);
  Epetra_Vector Y(View, A->RangeMap(), ap);
  
  A->Multiply(false, X, Y);
  
  return 1;
}

/******************************************************************************/

int Epetra_ML_getrow(void *data, int N_requested_rows, int requested_rows[], 
		    int allocated_space, int columns[], double values[],
		    int row_lengths[])
{
/*
 * GetRow function for matrix of type Epetra_CrsMatrix.
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
  Epetra_CrsMatrix  *A = (Epetra_CrsMatrix *) data;
  int nz_ptr = 0;
  int NumRows = A->NumMyRows();
  int NumEntries;
  double *Values;
  int *Indices;

  for (int i = 0; i < N_requested_rows; i++)
  {
    int LocalRow = requested_rows[i];
    int ierr = A->ExtractMyRowView(LocalRow, NumEntries, Values, Indices);
    if (ierr) return(1);
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


  Epetra_CrsMatrix *A = (Epetra_CrsMatrix *) data;

  if (A->Comm().NumProc()==1) return(1); // Nothing to do in serial mode.

  Epetra_Vector X_target(View, A->Importer()->TargetMap(), vec); //ghosted
  Epetra_Vector X_source(View, A->Importer()->SourceMap(), vec); //loc only

  assert(X_target.Import(X_source, *(A->Importer()),Insert)==0);

  return(1);
}
/******************************************************************************/

int EpetraMatrix2MLMatrix(ML *ml_handle, int level,
                         Epetra_CrsMatrix * A)
{
  int isize, osize;

  osize = A->NumMyRows();
  isize = osize;
  int N_ghost = A->NumMyCols() - A->NumMyRows();


  ML_Init_Amatrix(ml_handle, level,isize, osize, (void *) A);
  ML_Set_Amatrix_Getrow(ml_handle, level, Epetra_ML_getrow,
            Epetra_ML_comm_wrapper, isize+N_ghost);

  ML_Set_Amatrix_Matvec(ml_handle,  level, Epetra_ML_matvec);

  return 1;
}

#endif /*ifdef ML_WITH_EPETRA*/
