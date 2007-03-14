/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/************************************************************************/
/*          Utilities for Trilinos/ML users                             */
/*----------------------------------------------------------------------*/
/* Authors : Mike Heroux   (SNL)                                        */
/*           Jonathan Hu   (SNL)                                        */
/*           Ray Tuminaro  (SNL)                                        */
/*           Marzio Sala   (SNL)                                        */
/*           Michael Gee   (SNL)                                        */
/*           Chris Siefert (SNL)                                        */
/************************************************************************/

#include "ml_common.h"

#ifdef ML_WITH_EPETRA
#include <vector>
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Time.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "ml_FilterType.h"
#ifdef HAVE_ML_TEUCHOS
#include "Teuchos_ParameterList.hpp"
#endif

using namespace std;


// ====================================================================== 

typedef struct {
  ML_Epetra::FilterType Type;
  double AThresh;
  double RThresh;
  double FirstDivider;
  double SecondDivider;
  int Eqns;
  double* Mask;
} ML_Filter_Data;

static ML_Filter_Data Filter_;

// ====================================================================== 
int Epetra_ML_GetCrsDataptrs(ML_Operator *data, double **values, int **cols, int **rowptr)
{
  ML_Operator *mat_in;

  *values = NULL;
  *cols   = NULL;
  mat_in = (ML_Operator *) data;

  if ( (mat_in->matvec->func_ptr != ML_Epetra_matvec) && 
       (mat_in->matvec->func_ptr != ML_Epetra_CrsMatrix_matvec)) return 0;

  Epetra_RowMatrix *A = (Epetra_RowMatrix *) ML_Get_MyMatvecData(mat_in);

  Epetra_CrsMatrix * CrsA = NULL;

  CrsA = dynamic_cast<Epetra_CrsMatrix *>(A);
  if( CrsA != NULL ) 
    CrsA->ExtractCrsDataPointers(*rowptr, *cols, *values);

  return 0;
}

int ML_Epetra_matvec(ML_Operator *data, int in, double *p, int out, double *ap)
{
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;
  /* ML matvec wrapper for Epetra matrices. */

  // general case
  Epetra_RowMatrix *A = (Epetra_RowMatrix *) ML_Get_MyMatvecData(mat_in);

  // for VBR matrices first
  Epetra_VbrMatrix * VbrA = NULL;
  VbrA = dynamic_cast<Epetra_VbrMatrix *>(A);

  if( VbrA != NULL ) {
    Epetra_Vector X(View, VbrA->DomainMap(), p);
    Epetra_Vector Y(View, VbrA->RangeMap(), ap);
    VbrA->Multiply(false, X, Y);
  } else {   
    Epetra_Vector X(View, A->OperatorDomainMap(), p);
    Epetra_Vector Y(View, A->OperatorRangeMap(), ap);
  
    A->Multiply(false, X, Y);
  }

  return 1;
}

// ====================================================================== 

int ML_Epetra_CrsMatrix_matvec(ML_Operator *data, int in, double *p,
                                                  int out, double *ap)
{
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;
  /* ML matvec wrapper for Epetra matrices. */

  Epetra_CrsMatrix *A = (Epetra_CrsMatrix *) ML_Get_MyMatvecData(mat_in);

  Epetra_Vector X(View, A->OperatorDomainMap(), p);
  Epetra_Vector Y(View, A->OperatorRangeMap(), ap);
  
  A->Multiply(false, X, Y);

  return 1;
}

// ====================================================================== 

int ML_Epetra_VbrMatrix_matvec(ML_Operator *data, int in, double *p,
                                                  int out, double *ap)
{
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;
  /* ML matvec wrapper for Epetra matrices. */

  Epetra_VbrMatrix *A = (Epetra_VbrMatrix *) ML_Get_MyMatvecData(mat_in);

  Epetra_Vector X(View, A->DomainMap(), p);
  Epetra_Vector Y(View, A->RangeMap(), ap);
  A->Multiply(false, X, Y);

  return 1;
}

// ====================================================================== 

int ML_Epetra_matvec_Filter(ML_Operator *mat_in, int in, double *p, 
                            int out, double *ap)
{
  Epetra_RowMatrix *A = (Epetra_RowMatrix *) ML_Get_MyMatvecData(mat_in);
  int NumMyRows = A->NumMyRows();

  int row_lengths = 0;
  int allocated_space = A->MaxNumEntries();
  vector<int> columns(allocated_space + 1);
  vector<double> values(allocated_space + 1);

  // FIXME: DOES NOT WORK IN PARALLEL!
  assert (A->Comm().NumProc() == 1);
  
  for (int i = 0 ; i < NumMyRows ; ++i) {
    ap[i] = 0.0;
    int ierr;
    ierr = ML_Epetra_getrow_Filter(mat_in, 1, &i, allocated_space, 
                                   &columns[0], &values[0], &row_lengths);
    assert (ierr == 1);

    for (int j = 0 ; j < row_lengths ; ++j)
      ap[i] += values[j] * p[columns[j]];
  }

  return 1;
}

// ====================================================================== 
// General getrow for Epetra matrix classes.
// This function is deprecated, use one of the following instead:
// - ML_Epetra_RowMatrix_getrow
// - ML_Epetra_CrsMatrix_getrow
// - ML_Epetra_VbrMatrix_getrow
// ====================================================================== 

int ML_Epetra_getrow(ML_Operator *data, int N_requested_rows, int requested_rows[], 
		    int allocated_space, int columns[], double values[],
		    int row_lengths[])
{

  cout << "Fuction ML_Epetra_getrow() is no longer supported." << endl;
  cout << "You should use one of the following instead:" << endl;
  cout << "- ML_Epetra_RowMatrix_getrow();" << endl;
  cout << "- ML_Epetra_CrsMatrix_getrow();" << endl;
  cout << "- ML_Epetra_VbrMatrix_getrow()." << endl;
  cout << "If you don't know what is your matrix type, then use" << endl;
  cout << "the generic function for Epetra_RowMatrix's." << endl;
  cout << "You may need to update your Epetra wrapper and set the" << endl;
  cout << "appropriete function instead if ML_Epetra_getrow()" << endl;

  ML_RETURN(-1);

#if 0

  int nz_ptr = 0;
  int NumEntries;
  int MaxPerRow = 0;
  int NumPDEEqns=1;
  int * BlockIndices;
  Epetra_SerialDenseMatrix ** Entries;
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;

  Epetra_RowMatrix *Abase = (Epetra_RowMatrix *) ML_Get_MyGetrowData(mat_in);
  
  Epetra_CrsMatrix * Acrs = dynamic_cast<Epetra_CrsMatrix *>(Abase);
  int MatrixIsCrsMatrix = (Acrs!=0); // If this pointer is non-zero,
                                  // the cast to Epetra_CrsMatrix worked

  Epetra_VbrMatrix * Avbr = dynamic_cast<Epetra_VbrMatrix *>(Abase);
  int MatrixIsVbrMatrix = (Avbr!=0); // If this pointer is non-zero,
                                  // the cast to Epetra_VbrMatrix worked

  int *Indices;
  double *Values;
  int MatrixIsRowMatrix = false;
  if (MatrixIsCrsMatrix) {
    // do nothing for Crs
  } else  if( MatrixIsVbrMatrix ) {
    // for Vbr we need to know the number of PDE for each row
    if( Avbr->NumMyRows() % Avbr->NumMyBlockRows() != 0 ){
      cerr << "Error : NumPDEEqns does not seem to be constant\n";
      exit( EXIT_FAILURE );
    }
    NumPDEEqns = (Avbr->NumMyRows())/(Avbr->NumMyBlockRows());
  } else {
    // general RowMatrix case
    MatrixIsRowMatrix = true;
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
    else if (MatrixIsVbrMatrix) {
      // for vbr, we recover the local number of the BlockRow
      // (dividing by NumPDEEqns). In this way, we can get a view
      // of the local row (no memory allocation occurs).
      int PDEEqn = LocalRow % NumPDEEqns;
      int LocalBlockRow = LocalRow/NumPDEEqns;
      
      int RowDim;
      int NumBlockEntries;    
      ierr = Avbr->ExtractMyBlockRowView(LocalBlockRow,RowDim,
                     NumBlockEntries, BlockIndices, Entries);
      // I do here some stuff because Vbr matrices must
      // be treated differently.
      if (ierr) {
        if (MatrixIsRowMatrix) {
          delete [] Indices;
          delete [] Values;
        }
        return(0); 
      }
      NumEntries = NumBlockEntries*NumPDEEqns;
      if (nz_ptr + NumEntries > allocated_space) {
        if (MatrixIsRowMatrix) {
          delete [] Indices;
          delete [] Values;
        }
        return(0);
      }
      
      for( int j=0 ; j<NumBlockEntries ; ++j ) {
    for( int k=0 ; k<NumPDEEqns ; ++k ) {
      columns[nz_ptr] = BlockIndices[j]*NumPDEEqns+k;
      values[nz_ptr++] = (*Entries[j])(PDEEqn,k);
    }
      }
      row_lengths[i] = NumBlockEntries*NumPDEEqns;      
    }
    else
      ierr = Abase->ExtractMyRowCopy(LocalRow, MaxPerRow, NumEntries,
                                      Values, Indices);
    if (ierr) {
      if (MatrixIsRowMatrix) {
        delete [] Indices;
        delete [] Values;
      }
      return(0); //JJH I think this is the correct thing to return if
                 //    A->ExtractMyRowCopy returns something nonzero ..
    }

    if( !MatrixIsVbrMatrix ) {
      row_lengths[i] = NumEntries;
      if (nz_ptr + NumEntries > allocated_space) {
        if (MatrixIsRowMatrix) {
          delete [] Indices;
          delete [] Values;
        }
         return(0);
      }
      
      for (int j=0; j<NumEntries; j++) {
        columns[nz_ptr] = Indices[j];
        values[nz_ptr++] = Values[j];
      }
    }
  }

  if (MatrixIsRowMatrix) {
    delete [] Indices;
    delete [] Values;
  }
  
  return(1);
#endif
}

// ====================================================================== 
// Getrow for RowMatrix that are not Epetra_CrsMatrix or Epetra_VbrMatrix
// ====================================================================== 

int ML_Epetra_RowMatrix_getrow(ML_Operator *data, int N_requested_rows, 
                               int requested_rows[], int allocated_space, 
                               int columns[], double values[],
                               int row_lengths[])
{
  int nz_ptr = 0;
  int NumEntries;
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;

  Epetra_RowMatrix* A = (Epetra_RowMatrix *) ML_Get_MyGetrowData(mat_in);
  
  for (int i = 0; i < N_requested_rows; i++)
  {
    int ierr;
    int LocalRow = requested_rows[i];
    A->NumMyRowEntries(LocalRow, NumEntries);
    if (allocated_space < NumEntries)
      return(0); // to avoid Epetra print something on cout
    ierr = A->ExtractMyRowCopy(LocalRow, allocated_space, NumEntries,
                               values + nz_ptr, columns + nz_ptr);
    if (ierr) 
      return(0); //JJH I think this is the correct thing to return if
                 //    A->ExtractMyRowCopy returns something nonzero ..

    row_lengths[i] = NumEntries;
    // increase count of already used space...
    nz_ptr += NumEntries;
    // and decrease amount of available space
    allocated_space -= NumEntries;
    if (allocated_space < 0)
      return(0); // something was wrong here
  }

  return(1);
}

// ====================================================================== 
// Specialized getrow for Epetra_CrsMatrix class.
// ====================================================================== 

int ML_Epetra_CrsMatrix_getrow(ML_Operator *data, int N_requested_rows,
            int requested_rows[], 
		    int allocated_space, int columns[], double values[],
		    int row_lengths[])
{
  int nz_ptr = 0;
  int NumEntries;
  //int MaxPerRow = 0;
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;

  Epetra_CrsMatrix *Acrs =  (Epetra_CrsMatrix *) ML_Get_MyGetrowData(mat_in);
  
  for (int i = 0; i < N_requested_rows; i++)
  {
    int LocalRow = requested_rows[i];
    int *Indices;
    double *Values;

    int ierr = Acrs->ExtractMyRowView(LocalRow, NumEntries, Values, Indices);
    if (ierr)
      return(0); //JJH I think this is the correct thing to return if
                 //    A->ExtractMyRowCopy returns something nonzero ..

    row_lengths[i] = NumEntries;
    if (nz_ptr + NumEntries > allocated_space)
      return(0);
      
    for (int j=0; j<NumEntries; j++) {
      columns[nz_ptr] = Indices[j];
      values[nz_ptr++] = Values[j];
    }
  }

  return(1);
} //ML_Epetra_CrsMatrix_getrow

// ====================================================================== 
// Specialized getrow for Epetra_CrsMatrix class.
// ====================================================================== 

int ML_Epetra_CrsMatrix_get_one_row(ML_Operator *data, int N_requested_rows,
                                    int requested_rows[], 
                                    int allocated_space, int columns[], double values[],
                                    int row_lengths[])
{
  int nz_ptr = 0;
  int NumEntries;
  //int MaxPerRow = 0;
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;

  Epetra_CrsMatrix *Acrs =  (Epetra_CrsMatrix *) ML_Get_MyGetrowData(mat_in);
  
  for (int i = 0; i < N_requested_rows; i++)
  {
    int LocalRow = requested_rows[i];
    int *Indices;
    double *Values;

    int ierr = Acrs->ExtractMyRowView(LocalRow, NumEntries, Values, Indices);
    if (ierr)
      return(0); //JJH I think this is the correct thing to return if
                 //    A->ExtractMyRowCopy returns something nonzero ..

    row_lengths[i] = NumEntries;
    if (nz_ptr + NumEntries > allocated_space)
      return(0);
      
    for (int j=0; j<NumEntries; j++) {
      columns[nz_ptr] = Indices[j];
      values[nz_ptr++] = 1.0;
    }
  }

  return(1);
} //ML_Epetra_CrsMatrix_getrow


// ====================================================================== 
// Specialized getrow for Epetra_VbrMatrix class.
// ====================================================================== 

int ML_Epetra_VbrMatrix_getrow(ML_Operator *data, int N_requested_rows, 
                               int requested_rows[], int allocated_space, 
                               int columns[], double values[],
                               int row_lengths[])
{
  int nz_ptr = 0;
  int NumEntries;
  int * BlockIndices;
  Epetra_SerialDenseMatrix ** Entries;
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;
  Epetra_VbrMatrix * Avbr = (Epetra_VbrMatrix *) ML_Get_MyGetrowData(mat_in);

  /* moved into MultiLevelPreconditioner
  // for Vbr we need to know the number of PDE for each row
  if( Avbr->NumMyRows() % Avbr->NumMyBlockRows() != 0 ){
    cerr << "Error : NumPDEEqns does not seem to be constant\n";
    exit( EXIT_FAILURE );
  }
  */
  // NO! // int NumPDEEqns = mat_in->num_PDEs;
  // The above expression does not hold because of AmalgamateAndDropWeak,
  // which changes the value of mat_in->num_PDEs. Therefore we need
  // to use the one below, which might be slightly slower.
  int NumPDEEqns = (Avbr->NumMyRows())/(Avbr->NumMyBlockRows());

  for (int i = 0; i < N_requested_rows; i++)
  {
    int LocalRow = requested_rows[i];
    // for vbr, we recover the local number of the BlockRow
    // (dividing by NumPDEEqns). In this way, we can get a view
    // of the local row (no memory allocation occurs).
    int PDEEqn = LocalRow % NumPDEEqns;
    int LocalBlockRow = LocalRow/NumPDEEqns;
    
    int RowDim;
    int NumBlockEntries;    
    int ierr = Avbr->ExtractMyBlockRowView(LocalBlockRow,RowDim,
                   NumBlockEntries, BlockIndices, Entries);
    if (ierr) return(0); 
    NumEntries = NumBlockEntries*NumPDEEqns;
    if (nz_ptr + NumEntries > allocated_space)
      return(0);

    for( int j=0 ; j<NumBlockEntries ; ++j ) {
      for( int k=0 ; k<NumPDEEqns ; ++k ) {
        columns[nz_ptr] = BlockIndices[j]*NumPDEEqns+k;
        values[nz_ptr++] = (*Entries[j])(PDEEqn,k);
      }
    }
    row_lengths[i] = NumBlockEntries*NumPDEEqns;      
  }

  return(1);
} //ML_Epetra_VbrMatrix_getrow

// ====================================================================== 

#ifdef HAVE_ML_TEUCHOS
void ML_Set_Filter(Teuchos::ParameterList& List)
{
  Filter_.Type    = List.get("filter: type", ML_Epetra::ML_NO_FILTER);
  Filter_.AThresh = List.get("filter: absolute threshold", 0.0);
  Filter_.RThresh = List.get("filter: relative threshold", 1.0);
  Filter_.Eqns    = List.get("filter: equations", 1);
  Filter_.FirstDivider  = List.get("filter: first divider", 0);
  Filter_.SecondDivider = List.get("filter: second divider", 0);
  Filter_.Mask          = List.get("filter: mask", (double*)0);
}
#endif

// ====================================================================== 

int ML_Epetra_getrow_Filter(ML_Operator *data, int N_requested_rows, 
                            int requested_rows[], int allocated_space, 
                            int columns[], double values[], int row_lengths[])
{
  int ierr, eqn;
  ierr = ML_Epetra_getrow(data, N_requested_rows, requested_rows, allocated_space,
                          columns, values, row_lengths);

  if (ierr == 0)
    return(0);

  if (N_requested_rows != 1) {
    cerr << "Only N_requested_rows == 1 currently implemented..." << endl;
    exit(EXIT_FAILURE);
  }

  switch (Filter_.Type) {

  case ML_Epetra::ML_NO_FILTER:
    return(1);
    break;

  case ML_Epetra::ML_EQN_FILTER:

    for (int i = 0; i < row_lengths[0]; i++) {
      if (columns[i] % Filter_.Eqns != requested_rows[0] % Filter_.Eqns)
        values[i] = 0.0;
    }
    break;

  case ML_Epetra::ML_TWO_BLOCKS_FILTER:

    eqn = (requested_rows[0] % Filter_.Eqns);

    if (eqn < Filter_.FirstDivider) {
      // upper block
      for (int i = 0; i < row_lengths[0]; i++) {
        if (columns[i] % Filter_.Eqns >= Filter_.FirstDivider)
          // upper-right block
          values[i] = 0.0;
      }
    }
    else {
      // lower block
      for (int i = 0; i < row_lengths[0]; i++) {
        if (columns[i] % Filter_.Eqns < Filter_.FirstDivider)
          // lower-left block
          values[i] = 0.0;
      }
    }

    break;

  case ML_Epetra::ML_THREE_BLOCKS_FILTER:

    eqn = (requested_rows[0] % Filter_.Eqns);

    if (eqn < Filter_.FirstDivider) {
      for (int i = 0; i < row_lengths[0]; i++) {
        if (columns[i] % Filter_.Eqns >= Filter_.FirstDivider)
          values[i] = 0.0;
      }
    }
    else if (eqn < Filter_.SecondDivider) {
      for (int i = 0; i < row_lengths[0]; i++) {
        if (columns[i] % Filter_.Eqns <  Filter_.FirstDivider ||
            columns[i] % Filter_.Eqns >= Filter_.SecondDivider)
          values[i] = 0.0;
      }
    }
    else {
      for (int i = 0; i < row_lengths[0]; i++) {
        if (columns[i] % Filter_.Eqns <  Filter_.SecondDivider)
          values[i] = 0.0;
      }
    }

    break;

  case ML_Epetra::ML_MASK_FILTER:

    eqn = (requested_rows[0] % Filter_.Eqns);
    for (int i = 0; i < row_lengths[0]; i++) {
      values[i] *= Filter_.Mask[eqn * Filter_.Eqns + columns[i] % Filter_.Eqns];
    }
    break;

  default:

    cerr << "Error, file " << __FILE__ << ", line " << __LINE__ << endl;
    exit(EXIT_FAILURE);
  }

  if (Filter_.RThresh != 1.00 && Filter_.AThresh != 0.0) {
    for (int i = 0; i < row_lengths[0]; i++) {
      if (columns[i] == requested_rows[0]) {
        values[i] = Filter_.RThresh * values[i] + Filter_.AThresh * fabs(values[i]);
        break;
      }
    }
  }

  return(1);
}

// ====================================================================== 

int ML_Epetra_comm_wrapper(double vec[], void *data)
{
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

// ====================================================================== 

int ML_Epetra_CrsMatrix_comm_wrapper(double vec[], void *data)
{
  Epetra_CrsMatrix *A = (Epetra_CrsMatrix *) data;

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

// ====================================================================== 

int ML_Epetra_VbrMatrix_comm_wrapper(double vec[], void *data)
{
  Epetra_VbrMatrix *A = (Epetra_VbrMatrix *) data;

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

// ======================================================================

int ML_Operator_WrapEpetraMatrix(Epetra_RowMatrix * A, ML_Operator *newMatrix)
{
  int isize, osize;

  // FIXME ?? 
  //osize = A->NumMyRows();
  osize = A->OperatorRangeMap().NumMyElements();
  isize = A->OperatorDomainMap().NumMyElements();
  //  isize = A->NumMyCols();
  int N_ghost = A->RowMatrixColMap().NumMyElements() - isize;
  newMatrix->N_nonzeros = A->NumGlobalNonzeros();

  if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

  Epetra_CrsMatrix *Acrs = dynamic_cast<Epetra_CrsMatrix*>(A);

  if (Acrs) { // Epetra_CrsMatrix
    ML_Operator_Set_ApplyFuncData(newMatrix, isize, osize,
                                (void*) Acrs, osize,
                                NULL, 0);

    ML_CommInfoOP_Generate( &(newMatrix->getrow->pre_comm), 
                          ML_Epetra_CrsMatrix_comm_wrapper, (void *) Acrs, 
                          newMatrix->comm, isize, N_ghost);

    ML_Operator_Set_Getrow(newMatrix, newMatrix->outvec_leng,
                           ML_Epetra_CrsMatrix_getrow);

    ML_Operator_Set_ApplyFunc (newMatrix, ML_Epetra_CrsMatrix_matvec);

    newMatrix->type = ML_TYPE_CRS_MATRIX;
  }
  // TODO implement functionality for Epetra_VbrMatrix
  else { // RowMatrix
    ML_Operator_Set_ApplyFuncData(newMatrix, isize, osize,
                                (void*) A, osize,
                                NULL, 0);

    ML_CommInfoOP_Generate( &(newMatrix->getrow->pre_comm), 
                          ML_Epetra_comm_wrapper, (void *) A, 
                          newMatrix->comm, isize, N_ghost);

    ML_Operator_Set_Getrow(newMatrix, newMatrix->outvec_leng,
                           ML_Epetra_RowMatrix_getrow);

    ML_Operator_Set_ApplyFunc (newMatrix, ML_Epetra_VbrMatrix_matvec);

    newMatrix->type = ML_TYPE_ROW_MATRIX;
  }

  return 0;
}



/* This should (correctly) build the epetra maps from the ML_Operator object*/
void ML_Build_Epetra_Maps(ML_Operator* Amat,Epetra_Map **domainmap, Epetra_Map **rangemap){
  int    isize_offset, osize_offset;
  int Nghost;
  ML_Comm *comm;

  comm = Amat->comm;
#ifdef ML_MPI
  MPI_Comm mpi_comm ;
  mpi_comm = comm->USR_comm; 
  Epetra_MpiComm EpetraComm( mpi_comm ) ; 
#else
  Epetra_SerialComm EpetraComm ; 
#endif  

  Epetra_Time Time(EpetraComm);

  if (Amat->getrow->post_comm != NULL)  {
    if (Amat->comm->ML_mypid == 0)
      pr_error("Error: Please transpose matrix with ML_Operator_Transpose_byrow()\n       before calling ML_Operator2EpetraCrsMatrix().\n");
  }

  if (Amat->getrow->pre_comm == NULL)
    Nghost = 0;
  else {
    if (Amat->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Amat->getrow->pre_comm);
    Nghost = Amat->getrow->pre_comm->total_rcv_length;
  }

  int isize = Amat->invec_leng;
  int osize = Amat->outvec_leng;

  EpetraComm.ScanSum(&isize,&isize_offset,1); isize_offset-=isize;
  EpetraComm.ScanSum(&osize,&osize_offset,1); osize_offset-=osize;

  vector<double> global_isize; global_isize.resize(isize+Nghost+1);
  vector<int>    global_isize_as_int; global_isize_as_int.resize(isize+Nghost+1);
  vector<int>    global_osize_as_int(osize);
  
  for (int i = 0 ; i < isize; i++) {
          global_isize[i] = (double) (isize_offset + i);
          global_isize_as_int[i] = isize_offset + i;
  }
          
  for (int i = 0 ; i < osize; i++)
    global_osize_as_int[i] = osize_offset + i;
  
  for (int i = 0 ; i < Nghost; i++) global_isize[i+isize] = -1;
  
  *rangemap=new Epetra_Map( -1, osize, &global_osize_as_int[0], 0, EpetraComm ) ; 
  *domainmap=new Epetra_Map( -1, isize, &global_isize_as_int[0], 0, EpetraComm ) ; 
}/*end ML_Build_Epetra_Maps*/


// ================================================ ====== ==== ==== == = 
// This is a *ultra* lightweight wrap of an Epetra_CrsMatrix in ML.  This uses a
// "may change for experts only" function ExtractCrsDataPointers.
//
//  You need to have remapped the Epetra Matrix to include all the columns
//  before this routine gets called or else this won't work in parallel.
//
// -Chris Siefert 11/28/2006.
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
int ML_Operator_WrapEpetraCrsMatrix(Epetra_CrsMatrix * A, ML_Operator *newMatrix)
{
  int isize, osize,rv=0;
  osize = A->OperatorRangeMap().NumMyElements();
  isize = A->OperatorDomainMap().NumMyElements();
  
  int N_ghost = A->RowMatrixColMap().NumMyElements() - isize; 
  if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

  //  printf("[%d] Epetra->ML Local Size %dx%d\n",A->Comm().MyPID(),A->RowMap().NumMyElements(),A->ColMap().NumMyElements());
  
  /* Do the "View" Wrap */
  struct ML_CSR_MSRdata *epetra_csr= (struct ML_CSR_MSRdata*)malloc(sizeof(struct ML_CSR_MSRdata));
  epetra_csr->Nnz=newMatrix->N_nonzeros = A->NumGlobalNonzeros();
  epetra_csr->Nrows=osize;
  epetra_csr->Ncols=isize;  
  A->ExtractCrsDataPointers(epetra_csr->rowptr,epetra_csr->columns,epetra_csr->values);

  /* Sanity Check */
  if(!epetra_csr->rowptr || !epetra_csr->columns || !epetra_csr->values) rv=-1;
  
  /* Set the appropriate function pointers + data */
  ML_Operator_Set_ApplyFuncData(newMatrix, isize, osize,(void*) epetra_csr, osize,NULL,0);  
  ML_CommInfoOP_Generate(&(newMatrix->getrow->pre_comm),ML_Epetra_CrsMatrix_comm_wrapper, (void *) A, 
                          newMatrix->comm, isize, N_ghost);
  ML_Operator_Set_Getrow(newMatrix, newMatrix->outvec_leng,CSR_getrow);
  ML_Operator_Set_ApplyFunc (newMatrix, CSR_matvec);  
  newMatrix->data_destroy=free;
  newMatrix->type = ML_TYPE_CRS_MATRIX;  
  return rv;
}/*end ML_Operator_WrapEpetraCrsMatrix*/

// ================================================ ====== ==== ==== == 
// Thie provides a lightweight wrap of an ML_Operator in Epetra.  The Epetra
// object needs to be setup correctly beforehand or else this will have
// disasterous consequences.  We assume that the ML_Operator will persist until
// after the Epetra_CrsMatrix is destroyed, if this is set in View mode.
// -Chris Siefert 11/20/2006.
void Epetra_CrsMatrix_Wrap_ML_Operator(ML_Operator * A, const Epetra_Comm &Comm, const Epetra_Map &RowMap,Epetra_CrsMatrix **Result,Epetra_DataAccess CV){ 

  double bob;
  int mnz=10000;
  ML_Operator2EpetraCrsMatrix(A,*Result,mnz,false,bob);


  //  printf("[%d] ML->Epetra Wrap [global] %dx%d\n",Comm.MyPID(),(*Result)->NumGlobalRows(),(*Result)->NumGlobalCols());
  //  printf("[%d] ML->Epetra Wrap [local ] %dx%d\n",Comm.MyPID(),(*Result)->NumMyRows(),(*Result)->NumMyCols());
  
#ifdef THIS_CODE_DOESNT_WORK
  /* This is a very dangerous way to do this.  Live on the edge. */
  int *cols, *gcols;
  double* vals;
  struct ML_CSR_MSRdata* M_= (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(A);

  /* Build the Column Map */
  int *global_colmap;
  Epetra_Map *RowMap2, *DomainMap;
  ML_Build_Epetra_Maps(A,&DomainMap,&RowMap2);
  //  int num_local_cols=ML_build_global_numbering(A,&global_colmap,"cols");
  //  Epetra_Map ColMap(-1,num_local_cols,global_colmap,0,Comm);
  //  printf("[%d] ML->Epetra Local (R/C) %dx%d\n",Comm.MyPID(),RowMap2->NumMyElements(),ColMap->NumMyElements());
  
  /* Allocate the Epetra_CrsMatrix Object */
  //  *Result=new Epetra_CrsMatrix(CV,*RowMap2,*ColMap,0);
  //  *Result=new Epetra_CrsMatrix(CV,RowMap,ColMap,0);
  *Result=new Epetra_CrsMatrix(CV,*RowMap2,0);
  
  /* Fill the matrix. */
  for(int row=0;row<A->outvec_leng;row++){
    cols=&(M_->columns[M_->rowptr[row]]);
    vals=&(M_->values[M_->rowptr[row]]);
    (*Result)->InsertMyValues(row,M_->rowptr[row+1]-M_->rowptr[row],vals,cols);
  }/*end for*/
  //  (*Result)->FillComplete(*ColMap,*RowMap2);//hax
  (*Result)->FillComplete(*DomainMap,*RowMap2);//hax

  printf("[%d] ML->Epetra Wrap [global] %dx%d\n",Comm.MyPID(),(*Result)->NumGlobalRows(),(*Result)->NumGlobalCols());
  printf("[%d] ML->Epetra Wrap [local ] %dx%d\n",Comm.MyPID(),(*Result)->NumMyRows(),(*Result)->NumMyCols());

    
  /* Cleanup */
  delete DomainMap; delete RowMap2;
  //  free(global_colmap);
#endif
}/*end Epetra_CrsMatrix_Wrap_ML_Operator*/


// ============================================================================
int* ML_Epetra::FindLocalDiricheltRowsFromOnesAndZeros(const Epetra_CrsMatrix & Matrix, int &numBCRows){
  int *dirichletRows = new int[Matrix.NumMyRows()];
  numBCRows = 0;
  for (int i=0; i<Matrix.NumMyRows(); i++) {
    int numEntries, *cols;
    double *vals;
    int ierr = Matrix.ExtractMyRowView(i,numEntries,vals,cols);
    if (ierr == 0) {
      int nz=0;
      for (int j=0; j<numEntries; j++) if (vals[j] != 0.0) nz++;
      if (nz == 1) dirichletRows[numBCRows++] = i;      
    }/*end if*/
  }/*end fpr*/
  return dirichletRows;
}/*end FindLocalDiricheltRowsFromOnesAndZeros*/


// ====================================================================== 
 //! Finds Dirichlet the local Dirichlet columns, given the local Dirichlet rows
Epetra_IntVector * ML_Epetra::FindLocalDirichletColumnsFromRows(const int *dirichletRows, int numBCRows,const Epetra_CrsMatrix & Matrix){
  const Epetra_Map & ColMap = Matrix.ColMap();
  int indexBase = ColMap.IndexBase();
  Epetra_Map globalMap(Matrix.NumGlobalCols(),indexBase,Matrix.Comm());

  // create the exporter from this proc's column map to global 1-1 column map
  Epetra_Export Exporter(ColMap,globalMap);

  // create a vector of global column indices that we will export to
  Epetra_IntVector globColsToZero(globalMap);
  // create a vector of local column indices that we will export from
  Epetra_IntVector *myColsToZero= new Epetra_IntVector(ColMap);
  //  myColsToZero->PutScalar(0);
  myColsToZero->PutValue(0);  

  // for each local column j in a local dirichlet row, set myColsToZero[j]=1
  for (int i=0; i < numBCRows; i++) {
    int numEntries;
    double *vals;
    int *cols;
    Matrix.ExtractMyRowView(dirichletRows[i],numEntries,vals,cols);
    for (int j=0; j < numEntries; j++)
      (*myColsToZero)[ cols[j] ] = 1;
  }/*end for*/

  // export to the global column map
  globColsToZero.Export(*myColsToZero,Exporter,Add);
  // now import from the global column map to the local column map
  myColsToZero->Import(globColsToZero,Exporter,Insert);

  return myColsToZero;
}/*end FindLocalDirichletColumnsFromRows*/


  // ====================================================================== 
Epetra_IntVector * ML_Epetra::LocalRowstoColumns(int *Rows, int numRows,const Epetra_CrsMatrix & Matrix){
  const Epetra_Map & ColMap = Matrix.ColMap();
  int indexBase = ColMap.IndexBase();
  Epetra_Map globalMap(Matrix.NumGlobalCols(),indexBase,Matrix.Comm());

  // create the exporter from this proc's column map to global 1-1 column map
  Epetra_Export Exporter(ColMap,globalMap);

  // create a vector of global column indices that we will export to
  Epetra_IntVector globColsToZero(globalMap);
  // create a vector of local column indices that we will export from
  Epetra_IntVector *myColsToZero= new Epetra_IntVector(ColMap);
  myColsToZero->PutValue(0);  

  // flag all local columns corresponding to the local rows specified
  for (int i=0; i < numRows; i++) 
    (*myColsToZero)[Matrix.LCID(Matrix.GRID(Rows[i]))]=1;

  // export to the global column map
  globColsToZero.Export(*myColsToZero,Exporter,Add);
  // now import from the global column map to the local column map
  myColsToZero->Import(globColsToZero,Exporter,Insert);

  return myColsToZero;
}/*end LocalRowstoColumns*/

  // ====================================================================== 
void ML_Epetra::Apply_BCsToMatrixRows(const int *dirichletRows, int numBCRows, const Epetra_CrsMatrix & Matrix)
{
  /* This function zeros out *rows* of Matrix that correspond to Dirichlet rows.
     Input:
         Matrix             matrix
     Output:
         Grad               matrix with Dirichlet *rows* zeroed out

     Comments: The graph of Matrix is unchanged.
  */
  
  // -------------------------
  // now zero out the rows
  // -------------------------
  for (int i=0; i < numBCRows; i++) {
    int numEntries;
    double *vals;
    int *cols;
    Matrix.ExtractMyRowView(dirichletRows[i],numEntries,vals,cols);
    for (int j=0; j < numEntries; j++)
      vals[j] = 0.0;
  }/*end for*/
}/*end Apply_BCsToMatrixRows*/



// ====================================================================== 
void ML_Epetra::Apply_BCsToMatrixColumns(const Epetra_IntVector &dirichletColumns,const Epetra_CrsMatrix & Matrix){
  /* This function zeros out columns of Matrix.
     Input:
         dirichletColumns   outpuy from FindLocalDirichletColumnsFromRow
         Matrix             matrix to nuke columns of 
     Output:
         Matrix             matrix with columns zeroed out

     Comments: The graph of Matrix is unchanged.
  */
  for (int i=0; i < Matrix.NumMyRows(); i++) {
    int numEntries;
    double *vals;
    int *cols;
    Matrix.ExtractMyRowView(i,numEntries,vals,cols);
    for (int j=0; j < numEntries; j++) {
      if (dirichletColumns[ cols[j] ] > 0)  vals[j] = 0.0;
    }/*end for*/
  }/*end for*/
}/* end Apply_BCsToMatrixColumns */


// ====================================================================== 
void ML_Epetra::Apply_BCsToMatrixColumns(const int *dirichletRows, int numBCRows, const Epetra_CrsMatrix & Matrix){
  /* This function zeros out columns of Matrix.
     Input:
         dirichletRows      output from FindLocalDirichletRowsFromOnesAndZeros
         numBCRows          output from FindLocalDirichletRowsFromOnesAndZeros
         Matrix             matrix to nuke columns of 
     Output:
         Matrix             matrix with columns zeroed out

     Comments: The graph of Matrix is unchanged.
  */
  Epetra_IntVector* dirichletColumns=FindLocalDirichletColumnsFromRows(dirichletRows,numBCRows,Matrix);
  Apply_BCsToMatrixColumns(*dirichletColumns,Matrix);
  delete dirichletColumns;  
}/* end Apply_BCsToMatrixColumns */

// ====================================================================== 
void ML_Epetra::Apply_BCsToMatrixColumns(const Epetra_RowMatrix & iBoundaryMatrix, const Epetra_RowMatrix & iMatrix){
  const Epetra_CrsMatrix *BoundaryMatrix = dynamic_cast<const Epetra_CrsMatrix*> (&iBoundaryMatrix);
  const Epetra_CrsMatrix *Matrix = dynamic_cast<const Epetra_CrsMatrix*>(&iMatrix);

  if (BoundaryMatrix == 0 || Matrix == 0) {
    cout << "Not applying Dirichlet boundary conditions to gradient "
         << "because cast failed." << endl;
    return;
  }

  // locate Dirichlet edges
  int numBCRows;
  int *dirichletRows = FindLocalDiricheltRowsFromOnesAndZeros(*Matrix,numBCRows);
  Apply_BCsToMatrixColumns(dirichletRows,numBCRows,*Matrix);

  delete [] dirichletRows;
}/* end Apply_BCsToMatrixColumns */



// ====================================================================== 
void ML_Epetra::Remove_Zeroed_Rows(const Epetra_CrsMatrix & Matrix)
{
  /* Finds zeroed out rows and plops a 1 on the diagonal
     rows/columns to nuke.
     Input:
         Matrix             matrix
     Output:
         Grad               matrix with zerod out rows getting a one

     Comments: The graph of Matrix is unchanged.
  */
  int i,j,N=Matrix.NumMyRows();
  int numEntries, gridx,cidx;
  double *vals;
  int *cols;

  /* Zero the rows, add ones to diagonal */
  for (i=0; i < N; i++) {
    Matrix.ExtractMyRowView(i,numEntries,vals,cols);
    gridx=Matrix.GRID(i);
    for (j=0, cidx=-1; j < numEntries; j++){
      if(vals[j]!=0) break;
      if(gridx == Matrix.GCID(cols[j])) cidx=j;      
    }/*end for*/
    if(cidx!=-1) vals[cidx]=1.0;
  }/*end for*/
  
}/*end Apply_OAZToMatrix*/




// ====================================================================== 
void ML_Epetra::Apply_OAZToMatrix(int *dirichletRows, int numBCRows, const Epetra_CrsMatrix & Matrix)
{
  /* This function does row/column ones-and-zeros on a matrix, given the
     rows/columns to nuke.
     Input:
         Matrix             matrix
     Output:
         Grad               matrix with Dirichlet rows/columns OAZ'd.

     Comments: The graph of Matrix is unchanged.
  */

  int numEntries;
  double *vals;
  int *cols;

  /* Find the local column numbers to nuke */
  Epetra_IntVector *dirichletColumns=LocalRowstoColumns(dirichletRows,numBCRows,Matrix);

  /* Zero the columns */
  FILE *f=fopen("dcols.dat","w");
  for (int i=0; i < Matrix.NumMyRows(); i++) {
    Matrix.ExtractMyRowView(i,numEntries,vals,cols);
    for (int j=0; j < numEntries; j++) 
      if ((*dirichletColumns)[ cols[j] ] > 0){
        vals[j] = 0.0;
        fprintf(f,"%d\n",cols[j]);
      }
    
  }/*end for*/
  
  /* Zero the rows, add ones to diagonal */
  for (int i=0; i < numBCRows; i++) {
    Matrix.ExtractMyRowView(dirichletRows[i],numEntries,vals,cols);
    for (int j=0; j < numEntries; j++)
      if(Matrix.GRID(dirichletRows[i])==Matrix.GCID(cols[j])) vals[j]=1.0;
      else vals[j] = 0.0;
  }/*end for*/


  delete dirichletColumns;  
}/*end Apply_OAZToMatrix*/





// ====================================================================== 
void ML_Epetra::Apply_BCsToGradient(
             const Epetra_RowMatrix & iEdgeMatrix,
             const Epetra_RowMatrix & iGrad)
{
  /* This function zeros out *rows* of T that correspond to Dirichlet rows in
     the curl-curl matrix.  It mimics what was done previously in
     ML_Tmat_applyDirichletBC().
     Input:
         EdgeMatrix         curl-curl matrix
         Grad               gradient matrix
     Output:
         Grad               gradient matrix with *rows* zeroed out

     Comments: The graph of Grad is unchanged.
  */

  const Epetra_CrsMatrix *EdgeMatrix = dynamic_cast<const Epetra_CrsMatrix*>
                                           (&iEdgeMatrix );
  const Epetra_CrsMatrix *Grad = dynamic_cast<const Epetra_CrsMatrix*>(&iGrad );

  if (EdgeMatrix == 0 || Grad == 0) {
    cout << "Not applying Dirichlet boundary conditions to gradient "
         << "because cast failed." << endl;
    return;
  }

  // locate Dirichlet edges
  int *dirichletEdges = new int[EdgeMatrix->NumMyRows()];
  int numBCEdges = 0;
  for (int i=0; i<EdgeMatrix->NumMyRows(); i++) {
    int numEntries, *cols;
    double *vals;
    int ierr = EdgeMatrix->ExtractMyRowView(i,numEntries,vals,cols);
    if (ierr == 0) {
      int nz=0;
      for (int j=0; j<numEntries; j++) if (vals[j] != 0.0) nz++;
      if (nz == 1) {
        dirichletEdges[numBCEdges++] = i;
      }
    }
  }
  printf("Picking up %d Dirichlet rows\n",numBCEdges);

  
  // -------------------------
  // now zero out the rows
  // -------------------------
  for (int i=0; i < numBCEdges; i++) {
    int numEntries;
    double *vals;
    int *cols;
    Grad->ExtractMyRowView(dirichletEdges[i],numEntries,vals,cols);
    for (int j=0; j < numEntries; j++)
      vals[j] = 0.0;
  }
  delete [] dirichletEdges;
} //Apply_BCsToGradient

// ====================================================================== 

int ML_Epetra_CrsGraph_matvec(ML_Operator *data, int in, double *p,
                              int out, double *ap)
{
  cerr << "ML_Epetra_CrsGraph_matvec() not implemented." << endl;
  ML_RETURN(-1);
}

// ====================================================================== 

int ML_Epetra_CrsGraph_getrow(ML_Operator *data, int N_requested_rows,
                              int requested_rows[], int allocated_space, 
                              int columns[], double values[],
                              int row_lengths[])
{
  int nz_ptr = 0;
  int NumEntries;
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;

  Epetra_CrsGraph *Graph =  (Epetra_CrsGraph *) ML_Get_MyGetrowData(mat_in);
  
  for (int i = 0; i < N_requested_rows; i++)
  {
    int LocalRow = requested_rows[i];
    int *Indices;

    int ierr = Graph->ExtractMyRowView(LocalRow, NumEntries, Indices);
    if (ierr)
      return(0); //JJH I think this is the correct thing to return if
                 //    A->ExtractMyRowCopy returns something nonzero ..

    row_lengths[i] = NumEntries;
    if (nz_ptr + NumEntries > allocated_space)
      return(0);
      
    for (int j=0; j<NumEntries; j++) {
      columns[nz_ptr] = Indices[j];
      values[nz_ptr++] = 1.0; // simply set each entry to 1
    }
  }

  return(1);
} //ML_Epetra_CrsGraph_getrow

// ====================================================================== 
int ML_Epetra_CrsGraph_comm_wrapper(double vec[], void *data)
{
  Epetra_CrsGraph*A = (Epetra_CrsGraph*) data;

  if (A->Comm().NumProc()==1) return(1); // Nothing to do in serial mode.

  if( A->Importer() != 0 ) {
    // this is SLOW
    const Epetra_BlockMap& RowMap = A->RowMap(); // this is a block map
    const Epetra_BlockMap& ColMap = A->ColMap(); // this is a block map

    Epetra_Map RowMap2(-1, RowMap.NumMyElements(), RowMap.MyGlobalElements(), ColMap.IndexBase(), RowMap.Comm());
    Epetra_Map ColMap2(-1, ColMap.NumMyElements(), ColMap.MyGlobalElements(), ColMap.IndexBase(), ColMap.Comm());
    Epetra_Import Importer(ColMap2, RowMap2);

    Epetra_Vector X_target(View, 
                           ColMap2,
                           //A->Importer()->TargetMap(),
			   vec); //ghosted
    Epetra_Vector X_source(View, 
                           RowMap2, 
                           //A->Importer()->SourceMap(),
			   vec); //loc only
  
    X_target.Import(X_source, 
                     Importer,
                    //*(A->Importer()), 
                    Insert);
  }
  
  return(1);
}
// ======================================================================

int ML_Operator_WrapEpetraCrsGraph(Epetra_CrsGraph* Graph, ML_Operator *newMatrix)
{
  int isize, osize;

  osize = Graph->RangeMap().NumMyElements();
  isize = Graph->DomainMap().NumMyElements();
  assert (Graph->HaveColMap() == true);
  int N_ghost = Graph->NumMyBlockCols() - isize;

  if (N_ghost < 0) N_ghost = 0;

  ML_Operator_Set_ApplyFuncData(newMatrix, isize, osize,
                                (void*) Graph, osize,
                                NULL, 0);

  ML_CommInfoOP_Generate( &(newMatrix->getrow->pre_comm), 
                         ML_Epetra_CrsGraph_comm_wrapper, (void *) Graph, 
                         newMatrix->comm, isize, N_ghost);

  ML_Operator_Set_Getrow(newMatrix, newMatrix->outvec_leng,
                         ML_Epetra_CrsGraph_getrow);

  ML_Operator_Set_ApplyFunc (newMatrix, ML_Epetra_CrsGraph_matvec);

  return 0;
}

// ======================================================================

int EpetraMatrix2MLMatrix(ML *ml_handle, int level,
                         Epetra_RowMatrix * A)
{
  int isize, osize;

  osize = A->NumMyRows();
  isize = osize;
  int N_ghost = A->NumMyCols() - A->NumMyRows();

  if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

  ML_Init_Amatrix(ml_handle, level,isize, osize, (void *) A);
  ML_Set_Amatrix_Getrow(ml_handle, level, ML_Epetra_RowMatrix_getrow,
                        ML_Epetra_comm_wrapper, isize+N_ghost);

  ML_Set_Amatrix_Matvec(ml_handle,  level, ML_Epetra_matvec);

  return 1;
}

// ======================================================================
int ML_back_to_epetraCrs(ML_Operator *Mat1Mat2, ML_Operator *Result, 
			 ML_Operator *Mat1, ML_Operator *Mat2)
{
  //---------------------------------------------------------------------------
  Epetra_CrsMatrix *Mat1_epet = (Epetra_CrsMatrix *) Mat1->data;
  Epetra_CrsMatrix *Mat2_epet = (Epetra_CrsMatrix *) Mat2->data;
  
  //---------------------------------------------------------------------------
  // for temporary use create a linear row map, range map, domain map and a matrix
  Epetra_Map* linrowmap = new Epetra_Map(Mat1_epet->RowMap().NumGlobalElements(),
                                         Mat1_epet->RowMap().NumMyElements(),0,
                                         Mat1_epet->Comm());
  Epetra_Map* linrangemap = new Epetra_Map(Mat1_epet->OperatorRangeMap().NumGlobalElements(),
                                           Mat1_epet->OperatorRangeMap().NumMyElements(),0,
                                           Mat1_epet->Comm());
  Epetra_Map* lindomainmap = new Epetra_Map(Mat2_epet->OperatorDomainMap().NumGlobalElements(),
                                            Mat2_epet->OperatorDomainMap().NumMyElements(),0,
                                            Mat2_epet->Comm());
  Epetra_CrsMatrix* tmpresult = new Epetra_CrsMatrix(Copy,*linrowmap,500,false);
  

  // see results as ML_Operator
  //ML_Operator_Print(Mat2,"Mat2");
  //ML_Operator_Print(Mat1Mat2,"Mat1Mat2");
  
  // warning
  // when either Mat1 or Mat2 contain empty columns, this routine fails.
  // This is due to difering philosophies in ML and Epetra w.r.t
  // local column indices. In case of empty columns, ML keeps those
  // local column indices while Epetra strips them out.
  // This appears in serial and in parallel.
  // I currently do not see any elegant way to fix this easily.
  // (One fix proposed by Ray would be to add another layer of indirect addressing to
  // the ML_Epetra_CrsMatrix_getrow to fix this. But this will come at
  // some price and in most cases will not be needed)
  
  //---------------------------------------------------------------------------
  // The result matrix Mat1Mat2:
  // - ML created a new row numbering that is a linear map for the rows 
  //   no matter what the input maps were
  // - row indices are local
  // - ML created a column numbering matching the new row map
  // - col indices are global
  
  //---------------------------------------------------------------------------
  // fill the temporary matrix tmpresult which has linear maps as well
  int allocated = 0, row_length;
  int *bindx = NULL;
  double *val = NULL;
  int* global_rows = linrowmap->MyGlobalElements();
  if (Mat1Mat2->getrow->Nrows != Mat1_epet->RowMap().NumMyElements())
  {
    cout << "Rowmap of ML_Operator and Epetra_CrsMatrix are different!\n";
    exit(-1);
  }
  for (int i=0; i<Mat1Mat2->getrow->Nrows; ++i) 
  {
    // get the row
    ML_get_matrix_row(Mat1Mat2, 1, &i, &allocated, &bindx, &val,&row_length, 0);
    // ML pads empty rows with a zero, take these out again in the result
    if (row_length==1 && val[0]==0.0) continue;
    // the row index i is an ml linear map
    // the row index global_rows[i] is the true Epetra grid
    // we have columns bindx which are global but refer to MLs linear map
    // this matches the map of tmpresult
    int err = tmpresult->InsertGlobalValues(global_rows[i],row_length, 
                                            val,bindx);
    if (err!=0 && err != 1) cout << "tmpresult->InsertGlobalValues returned " << err << endl;
  }
  if (bindx != NULL) ML_free(bindx);
  if (val != NULL) ML_free(val);

  int err = tmpresult->FillComplete(*lindomainmap,*linrangemap);
  if (err) 
  {
    cerr <<"Error in Epetra_CrsMatrix FillComplete" << err << endl;
    EPETRA_CHK_ERR(err);
  }
  delete linrowmap;
  delete linrangemap;
  delete lindomainmap;

  //---------------------------------------------------------------------------
  // compute the global column lookup of the final result matrix
  // the unknown column map is an overlapping version of the Mat2_epet->OperatorDomainMap()
  const Epetra_Comm& comm = Mat2_epet->Comm();
  const Epetra_Map& dommap = Mat2_epet->OperatorDomainMap();
  vector<int> gcolumns(dommap.NumGlobalElements());
  int countold=0;
  int count=0;
  for (int proc=0; proc<comm.NumProc(); ++proc)
  {
    if (proc==comm.MyPID())
      for (int i=0; i<dommap.NumMyElements(); ++i)
      {
        //cout << "Proc " << proc << " gcolumns[" << countold << "+" << count << "] = " << dommap.GID(i) << endl;
        gcolumns[countold+count] = dommap.GID(i);
        if (gcolumns[countold+count]<0) cout << "Cannot find gcid for lcid\n";
        ++count;
      }
    comm.Broadcast(&count,1,proc);
    comm.Broadcast(&gcolumns[countold],count,proc);
    countold += count;
    count=0;
  }
  
  //if (comm.MyPID()==0)
  //for (int i=0; i<(int)gcolumns.size(); ++i) cout << "gcolumns[ " << i << "] = " << gcolumns[i] << endl;
  //---------------------------------------------------------------------------
  // create the final result matrix with the correct row map
  Epetra_CrsMatrix *Result_epet = new Epetra_CrsMatrix(Copy,Mat1_epet->RowMap(),
                                                       20,false);
  //---------------------------------------------------------------------------
  // fill the final result from the tmpresult
  vector<int> gcid(50);
  for (int i=0; i<tmpresult->NumMyRows(); ++i)
  {
    int lrid = i; // holds for both matrices
    int grid = Result_epet->GRID(i); // holds for Result_epet
    if (grid<0) cout << "Cannot find grid for lrid\n";
    int numindices;
    int* indices;
    double* values;
    int err = tmpresult->ExtractMyRowView(lrid,numindices,values,indices);
    if (err) cout << "tmpresult->ExtractMyRowView returned " << err << endl;
    // indices[j] are lcid which is what we wanted
    if (numindices>(int)gcid.size()) gcid.resize(numindices);
    for (int j=0; j<numindices; ++j)
    {
      // get the gcid in the tmpresult matrix
      // gcid in tmpresult is from the linear column map
      int tmpgcid = tmpresult->GCID(indices[j]);
      if (tmpgcid<0 || tmpgcid>=(int)gcolumns.size()) 
        cout << "Cannot find tmpgcid " << tmpgcid << " for lcid (out of range)\n";
      // get the gcid from the lookup vector
      gcid[j] = gcolumns[tmpgcid];
    }
    // insert row into final result matrix
    err = Result_epet->InsertGlobalValues(grid,numindices,values,&(gcid[0]));
    if (err != 0 && err != 1) cout << "Result_epet->InsertGlobalValues returned " << err << endl;
  }
  int ierr=Result_epet->FillComplete(Mat2_epet->OperatorDomainMap(),
                                     Mat1_epet->OperatorRangeMap());
  if (ierr!=0) {
    cerr <<"Error in Epetra_CrsMatrix FillComplete" << ierr << endl;
    EPETRA_CHK_ERR(ierr);
  }

  // tidy up
  delete tmpresult;
  gcolumns.clear();
  gcid.clear();
  
  // wrap the result
  ML_Operator_WrapEpetraMatrix((Epetra_RowMatrix *) Result_epet, Result);

  return 1;
}

// ======================================================================

Epetra_CrsMatrix *Epetra_MatrixMult(Epetra_RowMatrix *B_crs, Epetra_RowMatrix *Bt_crs)
{
  ML_Comm *comm, *temp;
  Epetra_RowMatrix *result;

  temp = global_comm;
  ML_Comm_Create(&comm);
  ML_Operator *B_ml, *Bt_ml, *BBt_ml;
  B_ml  = ML_Operator_Create(comm);
  Bt_ml = ML_Operator_Create(comm);
  BBt_ml  = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraMatrix(B_crs, B_ml);
  ML_Operator_WrapEpetraMatrix(Bt_crs, Bt_ml);
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

// ====================================================================== 

Epetra_CrsMatrix *Epetra_MatrixAdd(Epetra_RowMatrix *B_crs, Epetra_RowMatrix *Bt_crs, double scalar)
{
  ML_Comm *comm, *temp;

  temp = global_comm;
  ML_Comm_Create(&comm);
  ML_Operator *B_ml, *Bt_ml, *BBt_ml;
  B_ml  = ML_Operator_Create(comm);
  Bt_ml = ML_Operator_Create(comm);
  BBt_ml  = ML_Operator_Create(comm);
  ML_Operator_WrapEpetraMatrix(B_crs, B_ml);
  ML_Operator_WrapEpetraMatrix(Bt_crs, Bt_ml);
  Epetra_CrsMatrix *BBt_crs = new Epetra_CrsMatrix(Copy,
				            B_crs->RowMatrixRowMap(),
					    B_crs->RowMatrixColMap(), 0);
  BBt_ml->data = (void *) BBt_crs;
  ML_Operator_Add(B_ml, Bt_ml, BBt_ml, ML_EpetraCRS_MATRIX, scalar);

  BBt_crs->FillComplete(B_crs->OperatorRangeMap(),
                        B_crs->OperatorDomainMap());

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


#ifndef ML_CPP
#ifdef __cplusplus
extern "C" 
{
#endif
#endif

Epetra_CrsMatrix * Q  = NULL;
Epetra_FECrsMatrix *  Qt = NULL; 

// ======================================================================
ML_Operator * ML_BuildQ( int StartingNumElements,
			 int ReorderedNumElements,
			 int NumPDEEqns, int NullSpaceDim,
			 int * reordered_decomposition,
			 double * StartingNullSpace,
			 double * ReorderedNullSpace,
			 int ComputeNewNullSpace,
			 double * StartingBdry, double * ReorderedBdry,
			 USR_COMM mpi_communicator,
			 ML_Comm *ml_communicator ) 
{
  
  ML_Operator * ML_Q2;
  
#ifdef ML_MPI
  Epetra_MpiComm Comm( mpi_communicator );
#else
  Epetra_SerialComm Comm;
#endif
  
  Epetra_Map StartingMap(-1,StartingNumElements*NumPDEEqns,0,Comm);
  Epetra_Map ReorderedMap(-1,ReorderedNumElements*NumPDEEqns,0,Comm);
  
  Q = new Epetra_CrsMatrix(Copy,StartingMap,1);

  int * MyGlobalElements = StartingMap.MyGlobalElements();

  // fill Q
  for( int i=0 ; i<StartingNumElements ; i++ ) {
    // i and PointCol are for the amalagamated configuration
    double one = 1.0;
    int PointCol = reordered_decomposition[i];
    for( int j=0 ; j<NumPDEEqns ; ++j ) {
      // GlobalRow and GlobalCol are for the amalgamated conf
      int GlobalRow = MyGlobalElements[i*NumPDEEqns] + j;
      int GlobalCol = PointCol*NumPDEEqns + j;
      Q->InsertGlobalValues(GlobalRow, 1, &one, &GlobalCol );
      // It appears that this is the safest way to code
      // the Q operator. If we skip the diagonal values, then
      // the ML-epetra conversion generally crashes with
      // Zoltan aggregation. Appearantly, ParMETIS without
      // Zoltan does not require the diagonal element...
      // This is just slightly more expensive....
      GlobalRow = MyGlobalElements[i*NumPDEEqns] + j;
      GlobalCol = GlobalRow;
      double zero = 0.0;
      // NOTE: this function may return a warning
      // (if the element has already been inserted)
      Q->InsertGlobalValues(GlobalRow, 1, &zero, &GlobalCol );
    }
  }
  
  Q->FillComplete(ReorderedMap,StartingMap);
  
  {int itemp;
  Comm.MaxAll(&ComputeNewNullSpace,&itemp,1);
  if( itemp == 1 ) ComputeNewNullSpace = 1;
  }
  
  if( ComputeNewNullSpace == 1 ) {

    if( NumPDEEqns == NullSpaceDim ) {

      double ** StartArrayOfPointers = new double * [NullSpaceDim];
      double ** ReordArrayOfPointers = new double * [NullSpaceDim];
      
      for( int k=0 ; k<NullSpaceDim ; ++k ) {
	StartArrayOfPointers[k] = StartingNullSpace+k*StartingNumElements*NumPDEEqns;
	ReordArrayOfPointers[k] = ReorderedNullSpace+k*ReorderedNumElements*NumPDEEqns;
      }
      
      Epetra_MultiVector startNS(View,StartingMap,StartArrayOfPointers,NullSpaceDim);
      Epetra_MultiVector reordNS(View,ReorderedMap,ReordArrayOfPointers,NullSpaceDim);
      
      Q->Multiply(true,startNS,reordNS);
      
      delete [] StartArrayOfPointers;
      delete [] ReordArrayOfPointers;

    } else {
      
      Epetra_Vector startNS2(StartingMap);
      Epetra_Vector reordNS2(ReorderedMap);

      for( int i=0 ; i<NullSpaceDim ; ++i ) {
	startNS2.PutScalar(0.0);
	for( int j=0 ; j<StartingNumElements ; ++j ) {
	  startNS2[j] = StartingNullSpace[j*NullSpaceDim+i];
	}
	Q->Multiply(true,startNS2,reordNS2);
	for( int j=0 ; j<ReorderedNumElements ; ++j ) {
	  ReorderedNullSpace[j*NullSpaceDim+i] = reordNS2[j];
	}
      }
    }
  }

  double * Start = NULL;
  double * Reord = NULL;
  
  if( StartingNumElements != 0 ) Start = new double[StartingNumElements*NumPDEEqns];
  if( ReorderedNumElements != 0 ) Reord = new double[ReorderedNumElements*NumPDEEqns];

  Epetra_Vector xxx(View,StartingMap,Start);
  Epetra_Vector yyy(View,ReorderedMap,Reord);

  xxx.PutScalar(0.0);
  yyy.PutScalar(0.0);

  for( int i=0 ; i<StartingNumElements ; ++i ) {
    xxx[i*NumPDEEqns] = StartingBdry[i];
  }

  Q->Multiply(true,xxx,yyy);

  ML_Q2 = ML_Operator_Create( ml_communicator );  
  
  ML_Operator_WrapEpetraMatrix(Q, ML_Q2);

  for( int i=0 ; i<ReorderedNumElements ; ++i ) {
    ReorderedBdry[i] = yyy[i*NumPDEEqns];
  }
  
  if( Start != NULL ) delete [] Start;
  if( Reord != NULL ) delete [] Reord;

  return ML_Q2;
}

// ======================================================================
int ML_ApplyQ(int StartingNumElements,
	      int ReorderedNumElements,
	      int NumVectors,
	      double* StartingVectors,
	      double* ReorderedVectors)
{

  int NumPDEEqns = Q->OperatorRangeMap().NumMyElements() / StartingNumElements;

  if (NumPDEEqns == 1) {
 
    // in this case I can go on with pointers
    double** StartArrayOfPointers = new double * [NumVectors];
    double** ReordArrayOfPointers = new double * [NumVectors];

    for (int k = 0 ; k < NumVectors ; ++k) {
      StartArrayOfPointers[k] = StartingVectors + k * StartingNumElements;
      ReordArrayOfPointers[k] = ReorderedVectors + k * ReorderedNumElements;
    }

    Epetra_MultiVector startNS(View,Q->OperatorRangeMap(),
			       StartArrayOfPointers,NumVectors);
    Epetra_MultiVector reordNS(View,Q->OperatorDomainMap(),
			       ReordArrayOfPointers,NumVectors);
    Q->Multiply(true,startNS,reordNS);

    delete [] StartArrayOfPointers;
    delete [] ReordArrayOfPointers;

  }
  else {
    // here instead I must allocate, can be coded better
    assert (Q->OperatorRangeMap().NumMyElements() == StartingNumElements * NumPDEEqns);
    assert (Q->OperatorDomainMap().NumMyElements() == ReorderedNumElements * NumPDEEqns);

    Epetra_MultiVector startNS(Q->OperatorRangeMap(), NumVectors);
    Epetra_MultiVector reordNS(Q->OperatorDomainMap(), NumVectors);
    startNS.PutScalar(0.0);
    reordNS.PutScalar(0.0);

    for (int k = 0 ; k < NumVectors ; ++k) {
      for (int i = 0 ; i < StartingNumElements ; ++i) {
	startNS[k][i * NumPDEEqns] = StartingVectors[i + k * StartingNumElements];
      }
    }
    for (int k = 0 ; k < NumVectors ; ++k) {
      for (int i = 0 ; i < ReorderedNumElements ; ++i) {
	reordNS[k][i * NumPDEEqns] = ReorderedVectors[i + k * ReorderedNumElements];
      }
    }

    Q->Multiply(true,startNS,reordNS);

    for (int k = 0 ; k < NumVectors ; ++k) {
      for (int i = 0 ; i < ReorderedNumElements ; ++i) {
	ReorderedVectors[i + k * ReorderedNumElements] = reordNS[k][i * NumPDEEqns];
      }
    }

  }

  return 0;
}

void ML_DestroyQ(void) 
{

  delete Q;
  Q = NULL;

  return;
  
} /* ML_DestroyQ */

#if 0
// NOTE: this works ONLY if NumPDEEqns == 1. To be changed as
// done with ML_BuildQ for the general case

ML_Operator * ML_BuildQt( int StartingNumElements,
			  int ReorderedNumElements,
			  int reordered_decomposition[],
			  USR_COMM mpi_communicator,
			  ML_Comm *ml_communicator ) 
{
  
  ML_Operator * ML_Qt2;

  cout << "CHECK MEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE" << endl;
  exit( EXIT_FAILURE );
  
#ifndef ML_MPI
  /* ********************************************************************** */
  /* ONe should not call this function with one processor only (as he has   */
  /* nothing to redistributed. I simply return (this is also checked later) */
  /* ********************************************************************** */

  return NULL;
#else

  Epetra_MpiComm Comm( mpi_communicator );

  if( Comm.NumProc() == 1 ) return NULL;

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
  Qt->FillComplete(ReorderedMap,StartingMap);
  
  ML_Qt2 = ML_Operator_Create( ml_communicator );

  ML_Operator_WrapEpetraMatrix( Qt, ML_Qt2);

  return ML_Qt2;
#endif
  
} /* ML_BuildQt */

void ML_DestroyQt( void ) 
{

  delete Qt;
  Qt = NULL;

  return;
  
} /* ML_DestroyQt */
#endif

#ifndef ML_CPP
#ifdef __cplusplus
} /* extern "C" */
#endif
#endif

// ======================================================================
int ML_Operator2EpetraCrsMatrix(ML_Operator *Amat, Epetra_CrsMatrix * &
				CrsMatrix, int & MaxNumNonzeros,
				bool CheckNonzeroRow, double & CPUTime)
{
  int    isize_offset, osize_offset;
  int Nghost;
  ML_Comm *comm;

  comm = Amat->comm;
#ifdef ML_MPI
  MPI_Comm mpi_comm ;
  mpi_comm = comm->USR_comm; 
  Epetra_MpiComm EpetraComm( mpi_comm ) ; 
#else
  Epetra_SerialComm EpetraComm ; 
#endif  

  Epetra_Time Time(EpetraComm);

  if (Amat->getrow->post_comm != NULL)  {
    if (Amat->comm->ML_mypid == 0)
      pr_error("Error: Please transpose matrix with ML_Operator_Transpose_byrow()\n       before calling ML_Operator2EpetraCrsMatrix().\n");
  }

  if (Amat->getrow->pre_comm == NULL)
    Nghost = 0;
  else {
    if (Amat->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Amat->getrow->pre_comm);
    Nghost = Amat->getrow->pre_comm->total_rcv_length;
  }

  int isize = Amat->invec_leng;
  int osize = Amat->outvec_leng;

  EpetraComm.ScanSum(&isize,&isize_offset,1); isize_offset-=isize;
  EpetraComm.ScanSum(&osize,&osize_offset,1); osize_offset-=osize;

  vector<double> global_isize; global_isize.resize(isize+Nghost+1);
  vector<int>    global_isize_as_int; global_isize_as_int.resize(isize+Nghost+1);

  //vector<double> global_osize(osize);
  vector<int>    global_osize_as_int(osize);
  
  for (int i = 0 ; i < isize; i++) {
          global_isize[i] = (double) (isize_offset + i);
          global_isize_as_int[i] = isize_offset + i;
  }
          
  for (int i = 0 ; i < osize; i++) {
    //global_osize[i] = (double) (osize_offset + i);
    global_osize_as_int[i] = osize_offset + i;
  }
  for (int i = 0 ; i < Nghost; i++) global_isize[i+isize] = -1;
  
  Epetra_Map  rangemap( -1, osize, &global_osize_as_int[0], 0, EpetraComm ) ; 
  Epetra_Map  domainmap( -1, isize, &global_isize_as_int[0], 0, EpetraComm ) ; 
  
  CrsMatrix = new Epetra_CrsMatrix( Copy, rangemap, 0 ); 
  
  ML_exchange_bdry(&global_isize[0],Amat->getrow->pre_comm, 
 		 Amat->invec_leng,comm,ML_OVERWRITE,NULL);

  for ( int j = isize; j < isize+Nghost; j++ ) { 
    global_isize_as_int[j] = (int) global_isize[j];
  }

  // MS // introduced variable allocation for colInd and colVal
  // MS // improved efficiency in InsertGlobalValues

  int allocated = 128;
  vector<int> colInd(allocated);
  vector<double> colVal(allocated);
  int NumNonzeros;
  int ierr;
  int    ncnt;

  MaxNumNonzeros=0;
  
  for (int i = 0; i < osize; i++)
  {
    ierr = ML_Operator_Getrow(Amat,1,&i,allocated,&colInd[0],&colVal[0],&ncnt);

    if( ierr == 0 ) {
      do {
	allocated *= 2;
        if (allocated > 20000) // would look strange to have such a dense row
        {
          cerr << "Row " << i << " on processor " << comm->ML_mypid;
          cerr << " seems to have more than 20000 nonzeros." << endl;
          cerr << "This looks suspicious, so now I abort..." << endl;
          ML_EXIT(-1);
        }
        colInd.resize(allocated);
        colVal.resize(allocated);
	ierr = ML_Operator_Getrow(Amat,1,&i,allocated,&colInd[0],&colVal[0],&ncnt);
      } while( ierr == 0 );
    }

    // MS // check out how many nonzeros we have
    // MS // NOTE: this may result in a non-symmetric pattern for CrsMatrix

    NumNonzeros = 0;
    for (int j = 0; j < ncnt; j++) {
      if (colVal[j] != 0.0) {
	    colInd[NumNonzeros] = global_isize_as_int[colInd[j]];
	    colVal[NumNonzeros] = colVal[j];
	    NumNonzeros++;
      }
    }
    if( NumNonzeros == 0 && CheckNonzeroRow ) {
      cout << "*ML*WRN* in ML_Operator2EpetraCrsMatrix : \n*ML*WRN* Global row "
	   << global_osize_as_int[i]
	   << " has no nonzero elements (and " << ncnt
	   << " zero entries)" << endl
	   << "*ML*WRN* Now put 1 on the diagonal...\n";
      // insert a 1 on the diagonal
      colInd[NumNonzeros] = global_isize_as_int[i];
      colVal[NumNonzeros] = 1.0;
      NumNonzeros++;
    }
    MaxNumNonzeros = EPETRA_MAX(NumNonzeros,MaxNumNonzeros);
    
    CrsMatrix->InsertGlobalValues(global_osize_as_int[i], NumNonzeros, 
				  &colVal[0], &colInd[0]);
  }

  CrsMatrix->FillComplete(domainmap,rangemap);
  CrsMatrix->OptimizeStorage();

  CPUTime = Time.ElapsedTime();

  return 0;
  
} /* ML_Operator2EpetraCrsMatrix for rectangular matrices*/

#if 0
// FIXME: delete me ???
int ML_Operator2EpetraCrsMatrix_old(ML_Operator *Ke, Epetra_CrsMatrix * &
				CrsMatrix, int & MaxNumNonzeros,
				bool CheckNonzeroRow, double & CPUTime)
{

  double *global_nodes, *global_rows;
  int *global_rows_as_int, *global_nodes_as_int;
  int    Nnodes, node_offset, row_offset;
  int Nnodes_global, Nrows_global;
  int Nghost_nodes;
  int Nrows;
  ML_Comm *comm;

  comm = Ke->comm;
#ifdef ML_MPI
  MPI_Comm mpi_comm ;
  mpi_comm = comm->USR_comm; 
  Epetra_MpiComm EpetraComm( mpi_comm ) ; 
#else
  Epetra_SerialComm EpetraComm ; 
#endif  

  Epetra_Time Time(EpetraComm);

  if (Ke->getrow->pre_comm == NULL) 
    Nghost_nodes = 0;
  else {
    if (Ke->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Ke->getrow->pre_comm);
    Nghost_nodes = Ke->getrow->pre_comm->total_rcv_length;
  }

  Nnodes = Ke->invec_leng;
  Nrows = Ke->outvec_leng;

  assert( Nnodes == Nrows );

  // MS // moved to Epetra node_offset = ML_gpartialsum_int(Nnodes, comm);
  // Nnodes_global = Nnodes;
  // ML_gsum_scalar_int(&Nnodes_global, &dummy, comm);

  // MS moved to Epetra row_offset = ML_gpartialsum_int(Nrows, comm);
  //  Nrows_global = Nrows;
  // ML_gsum_scalar_int(&Nrows_global, &dummy, comm);

  EpetraComm.ScanSum(&Nnodes,&node_offset,1); node_offset-=Nnodes;
  EpetraComm.ScanSum(&Nrows,&row_offset,1); row_offset-=Nrows;

  EpetraComm.SumAll(&Nnodes,&Nnodes_global,1);
  EpetraComm.SumAll(&Nrows,&Nrows_global,1);  

  assert( Nnodes_global == Nrows_global ) ; 

  global_nodes = new double[Nnodes+Nghost_nodes+1];
  global_nodes_as_int = new int[Nnodes+Nghost_nodes+1];

  global_rows = new double[Nrows+1];
  global_rows_as_int = new int[Nrows+1];
  
  for (int i = 0 ; i < Nnodes; i++) global_nodes[i] = (double) (node_offset + i);
  for (int i = 0 ; i < Nrows; i++) {
    global_rows[i] = (double) (row_offset + i);
    global_rows_as_int[i] = row_offset + i;
  }
  for (int i = 0 ; i < Nghost_nodes; i++) global_nodes[i+Nnodes] = -1;
  
  Epetra_Map  EpetraMap( Nrows_global, Nrows, global_rows_as_int, 0, EpetraComm ) ; 
  
  CrsMatrix = new Epetra_CrsMatrix( Copy, EpetraMap, 0 ); 
  
  ML_exchange_bdry(global_nodes,Ke->getrow->pre_comm, 
 		 Ke->invec_leng,comm,ML_OVERWRITE,NULL);

  for ( int j = 0; j < Nnodes+Nghost_nodes; j++ ) { 
    global_nodes_as_int[j] = (int) global_nodes[j];
  }

  // MS // introduced variable allocation for colInd and colVal
  // MS // improved efficiency in InsertGlobalValues

  {
  int allocated = 1;
  int * colInd = new int[allocated];
  double * colVal = new double[allocated];
  int NumNonzeros;
  int ierr;
  int    ncnt;

  MaxNumNonzeros=0;
  
  for (int i = 0; i < Nrows; i++) {
    ierr = ML_Operator_Getrow(Ke,1,&i,allocated,colInd,colVal,&ncnt);

    if( ierr == 0 ) {
      do {
	delete [] colInd;
	delete [] colVal;
	allocated *= 2;
	colInd = new int[allocated];
	colVal = new double[allocated];
	ierr = ML_Operator_Getrow(Ke,1,&i,allocated,colInd,colVal,&ncnt);
      } while( ierr == 0 );
    }

    // MS // check out how many nonzeros we have
    // MS // NOTE: this may result in a non-symmetric patter for CrsMatrix

    NumNonzeros = 0;
    for (int j = 0; j < ncnt; j++) {
      if (colVal[j] != 0.0) {
	colInd[NumNonzeros] = global_nodes_as_int[colInd[j]];
	colVal[NumNonzeros] = colVal[j];
	NumNonzeros++;
      }
    }
    if( NumNonzeros == 0 && CheckNonzeroRow ) {
      cout << "*ML*WRN* in ML_Operator2EpetraCrsMatrix : \n*ML*WRN* Global row "
	   << global_rows_as_int[i]
	   << " has no nonzero elements (and " << ncnt
	   << " zero entries)" << endl
	   << "*ML*WRN* Now put 1 on the diagonal...\n";
      // insert a 1 on the diagonal
      colInd[NumNonzeros] = global_nodes_as_int[i];
      colVal[NumNonzeros] = 1.0;
      NumNonzeros++;
    }
    MaxNumNonzeros = EPETRA_MAX(NumNonzeros,MaxNumNonzeros);
    
    CrsMatrix->InsertGlobalValues( global_rows_as_int[i], NumNonzeros, 
				   colVal, colInd);
    
    //    CrsMatrix->InsertGlobalValues( global_rows_as_int[i], ncnt, 
    //				   colVal, colInd);
  }

  delete [] colInd;
  delete [] colVal;
  }
  
  delete [] global_nodes_as_int;
  delete [] global_rows_as_int;
  delete [] global_rows;
  delete [] global_nodes;

  CrsMatrix->FillComplete();

  CPUTime = Time.ElapsedTime();

  return 0;
  
} /* ML_Operator2EpetraCrsMatrix */
#endif

#ifdef WKC
int ML_Epetra_matvec_WKC (ML_Operator *data, int in, double *p, int out, double *ap)
{
  ML_Operator *mat_in;

  mat_in = data;
  Epetra_RowMatrix *A = (Epetra_RowMatrix *) ML_Get_MyMatvecData(mat_in);
  Epetra_MultiVector &X(*(Epetra_MultiVector *)p);
  Epetra_MultiVector &Y(*(Epetra_MultiVector *)ap);

  A->Multiply(false, X, Y);

  return 1;
}
#endif

// ============================================================================
#include "Epetra_FECrsMatrix.h"
// FIXME: change my name?
Epetra_FECrsMatrix* FakeMatrix = 0;


#ifndef ML_CPP
#ifdef __cplusplus
extern "C" 
{
#endif
#endif

int ML_Operator_DiscreteLaplacian(ML_Operator* Op, int SymmetricPattern,
				  double* x_coord, double* y_coord,
				  double* z_coord, double theta,
				  ML_Operator** NewOp)
{

  if (Op->invec_leng != Op->outvec_leng) 
    return(-1); // works only with square matrices

  int NumMyRows = Op->outvec_leng;
  int NumPDEEqns = 1; // FIXME or DELETEME
  int NumDimensions = 1;

  if (x_coord == 0)
    return(-2); // at least one dimension is required

  if (y_coord == 0 && z_coord != 0)
    return(-3); // cannot be this

  if (y_coord != 0)
    NumDimensions++;

  if (z_coord != 0)
    NumDimensions++;

  // need an Epetra_Comm object. This in general exists, but I
  // don't know an easy way to pass it up to here
#ifdef ML_MPI
  Epetra_MpiComm Comm(Op->comm->USR_comm);
#else
  Epetra_SerialComm Comm;
#endif

  // need to create a map. This may may already exist for the finest-level
  // matrix, but in general it does not

  Epetra_Map Map(-1,NumMyRows,0,Comm);
  int NumGlobalRows = Map.NumGlobalElements();
  
  // create the auxiliary matrix as VBR matrix. This should help
  // to save memory with respect to the creation of "pure" VBR
  // matrices.

  int MaxNnzRow = Op->max_nz_per_row;
  assert(MaxNnzRow > 0);
  // FIXME: I can compute it in this case

  FakeMatrix = new Epetra_FECrsMatrix(Copy,Map, MaxNnzRow);
  
  if (FakeMatrix == 0)
    return(-10); // problems occur

  if (ML_Get_PrintLevel() > 8) {
    cout << endl;
    cout << "Creating discrete Laplacian..." << endl;
    cout << "Number of dimensions = " << NumDimensions << endl;
    cout << "theta = " << theta;
    if (SymmetricPattern == 1)
      cout << ", using symmetric pattern" << endl;
    else
      cout << ", using original pattern" << endl;
    cout << endl;
  }

  // create the auxiliary matrix

  int allocated = 1;
  int * colInd = new int[allocated];
  double * colVal = new double[allocated];
  int NnzRow;

  double coord_i[3], coord_j[3];
  for( int i = 0; i<3 ; ++i ) {
    coord_i[i] = 0.0; coord_j[i] = 0.0;
  }

  // get global column number
 
  int Nghost;

  if (Op->getrow->pre_comm == NULL) 
    Nghost = 0;
  else {
    if (Op->getrow->pre_comm->total_rcv_length <= 0)
      ML_CommInfoOP_Compute_TotalRcvLength(Op->getrow->pre_comm);
    Nghost = Op->getrow->pre_comm->total_rcv_length;
  }

  vector<double> global;        
  vector<int>    global_as_int; 
  global.resize(NumMyRows+Nghost+1);
  global_as_int.resize(NumMyRows+Nghost+1);

  int offset;
  Comm.ScanSum(&NumMyRows,&offset,1); 
  offset -= NumMyRows;

  for (int i = 0 ; i < NumMyRows; i++) {
    global[i] = (double) (offset + i);
    global_as_int[i] = offset + i;
  }

  for (int i = 0 ; i < Nghost; i++) 
    global[i+NumMyRows] = -1;

  ML_exchange_bdry(&global[0],Op->getrow->pre_comm,
                  Op->invec_leng,Op->comm,ML_OVERWRITE,NULL);

  for ( int j = 0; j < NumMyRows+Nghost; j++ ) {
    global_as_int[j] = (int) global[j];
  }
  
  // =================== //
  // cycle over all rows //
  // =================== //

  for (int i = 0; i < NumMyRows ; i += NumPDEEqns ) {

    int GlobalRow = global_as_int[i];

    assert( GlobalRow != -1 );

    if( i%NumPDEEqns == 0 ) { // do it just once for each block row
      switch( NumDimensions ) {
      case 3:
	coord_i[2] = z_coord[i/NumPDEEqns];
      case 2:
	coord_i[1] = y_coord[i/NumPDEEqns];
      case 1:
	coord_i[0] = x_coord[i/NumPDEEqns];
      }
    }

    int ierr = ML_Operator_Getrow(Op,1,&i,allocated,colInd,colVal,&NnzRow);

    if( ierr == 0 ) {
      do {
	delete [] colInd;
	delete [] colVal;
	allocated *= 2;
	colInd = new int[allocated];
	colVal = new double[allocated];
	ierr = ML_Operator_Getrow(Op,1,&i,allocated,colInd,colVal,&NnzRow);
      } while( ierr == 0 );
    }

    // NOTE: for VBR matrices, the "real" value that will be used in
    // the subsequent part of the code is only the one for the first
    // equations. For each block, I replace values with the sum of
    // the abs of each block entry.

    for (int j = 0 ; j < NnzRow ; j += NumPDEEqns) {
      colVal[j] = fabs(colVal[j]);
      for (int k = 1 ; k < NumPDEEqns ; ++k) {
	colVal[j] += fabs(colVal[j+k]);
      }
    }

    // work only on the first equations. Theta will blend the
    // coordinate part with the sub of abs of row elements.

    int GlobalCol;
    double total = 0.0;

    for (int j = 0 ; j < NnzRow ; j += NumPDEEqns) {

     if (colInd[j]%NumPDEEqns == 0) { 

      // insert diagonal later
      if (colInd[j] != i) {
	if (colInd[j]%NumPDEEqns == 0) { // do it only once
	  switch(NumDimensions) {
	  case 3:
	    coord_j[2] = z_coord[colInd[j]/NumPDEEqns];
	  case 2:
	    coord_j[1] = y_coord[colInd[j]/NumPDEEqns];
	  case 1:
	    coord_j[0] = x_coord[colInd[j]/NumPDEEqns];
	  }
	}

	double tmp1=coord_i[0]-coord_j[0];
	double tmp2=coord_i[1]-coord_j[1];
	double tmp3=coord_i[2]-coord_j[2];
	double d2 = tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3;
	if( d2 == 0.0 ) {
	  cerr << endl;
	  cerr << "distance between node " << i/NumPDEEqns << " and node " 
	    << colInd[j]/NumPDEEqns << endl
	    << "is zero. Coordinates of these nodes are" << endl
	    << "x_i = " << coord_i[0] << ", x_j = " << coord_j[0] << endl  
	    << "y_i = " << coord_i[1] << ", y_j = " << coord_j[1] << endl  
	    << "z_i = " << coord_i[2] << ", z_j = " << coord_j[2] << endl  
	    << "Now proceeding with distance = 1.0" << endl;
	  cerr << endl;
	  d2 = 1.0;
	}

	double val = -(1.0-theta)*(1.0/d2) + theta*(colVal[j]);

	GlobalCol = global_as_int[colInd[j]];
	assert( GlobalCol != -1 );
	
	// insert this value on all rows
	for( int k=0 ; k<NumPDEEqns ; ++k ) {
	  int row = GlobalRow+k;
	  int col = GlobalCol+k;
	  if( row >= NumGlobalRows || col >= NumGlobalRows ) {
	    cerr << "trying to insert element (" << row 
	         << "," << col << "), " << endl
		 << "while NumGlobalRows = " << NumGlobalRows << endl
		 << "(GlobalRow = " << GlobalRow << ", GlobalCol = " << GlobalCol << ")" << endl
		 << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
	  }
	    
	  // FakeMatrix->InsertGlobalValues(row,1,&val,&col);
	  if( FakeMatrix->SumIntoGlobalValues(1,&row,1,&col,&val) != 0 ) {
	    int ierr = FakeMatrix->InsertGlobalValues(1,&row,1,&col,&val);
	    if( ierr ) {
	      cerr << "InsertGlobalValues return value = " << ierr << endl
		<< "for element (" << row << "," << col << ")" << endl
		<< "(file " << __FILE__ << ", line " << __LINE__
		<< ")" << endl;
	    }
	  }
	}

	total -= val;

	// put (j,i) element as well. this works also for
	// off-process elements. 
	if( SymmetricPattern == 1 ) {
	  for( int k=0 ; k<NumPDEEqns ; ++k ) {
	    int row = GlobalCol+k;
	    int col = GlobalRow+k;
	    if( row >= NumGlobalRows || col >= NumGlobalRows ) {
	      cerr << "trying to insert element (" << row 
		<< "," << col << "), " << endl
		<< "while NumGlobalRows = " << NumGlobalRows << endl
		<< "(GlobalRow = " << GlobalRow << ", GlobalCol = " << GlobalCol << ")" << endl
		<< "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
	    }
	    if( FakeMatrix->SumIntoGlobalValues(1,&row,1,&col,&val) != 0 ) { 
	      int ierr = FakeMatrix->InsertGlobalValues(1,&row,1,&col,&val);
	      if( ierr ) {
		cerr << "InsertGlobalValues return value = " << ierr << endl
		  << "for element (" << row << "," << col << ")" << endl
		  << "(file " << __FILE__ << ", line " << __LINE__
		  << ")" << endl;
	      }
	    }
	  }
	  total -= val;
	}
      } 
    }

    }

    // create lines with zero-row sum
    for( int k=0 ; k<NumPDEEqns ; ++k ) {
      int row = GlobalRow+k;
      assert( row < NumGlobalRows );
      if( FakeMatrix->SumIntoGlobalValues(1,&row,1,&row,&total) != 0) {
	int ierr = FakeMatrix->InsertGlobalValues(1,&row,1,&row,&total);
	if( ierr ) {
	  cerr << "InsertGlobalValues return value = " << ierr << endl
	    << "for element (" << row << "," << row << ")" << endl
	    << "(file " << __FILE__ << ", line " << __LINE__
	    << ")" << endl;
	}
      }
   }
  }

  FakeMatrix->GlobalAssemble();

  delete [] colInd;
  delete [] colVal;

  // create a new ML_Operator from this Epetra Matrix.

  *NewOp = ML_Operator_Create(Op->comm);
  ML_Operator_WrapEpetraMatrix(FakeMatrix,*NewOp);
  
  return(0);
}

int ML_Operator_Destroy_DiscreteLaplacian() 
{

  if (FakeMatrix != 0) {
    delete FakeMatrix;
    FakeMatrix = 0;
  }

  return 0;

}

#ifndef ML_CPP
#ifdef __cplusplus
} // extern "C"
#endif
#endif


string ML_toString(const int& x) {
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

string ML_toString(const double& x) {
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}


/*----------------------------------------------------------------------*
 |                                                           m.gee 03/05|
 |                                                                      |
 | reads an update vector from file and creates an Epetra_Map from it.  |
 |                                                                      |
 | the file has the following format:                                   |
 | first line: <number_global_elements>  <number_of_procs>              |
 | following lines: <rownumber> <proc_number> -1                        |
 |                                                                      |
 | the file has to have <number_global_elements> + 1 rows               |
 |                                                                      |
 | Input:  char* filename        name of file to read from              |
 |         Epetra_Comm& comm     a valid Epetra_Comm                    |
 |                                                                      |
 | Output: Epetra_Map*           an allocated Epetra_Map class          |
 |                               (the calling user is responsible       |
 |                                for destroying this)                  |
 |                                                                      |
 | returns Epetra_Map* on success, NULL otherwise                       |
 |                                                                      |
 |                                                                      |
 *----------------------------------------------------------------------*/
Epetra_Map* Epetra_ML_readupdatevector(char* filename, Epetra_Comm& comm)
{
  char  buffer[200];
  char* bptr      = 0;
  int numeq_total = 0;
  int numeq       = 0;
  Epetra_Map* map = 0;
  int proc        = comm.MyPID();
  int nproc       = comm.NumProc();
  
  FILE *fp = fopen(filename,"r");
  if (!fp) return 0;
  if (proc) 
  {
    fclose(fp);
    fp = 0;
  }

  int ok = 1;
  if (proc==0)
  {
     fgets(buffer,199,fp);
     numeq_total = strtol(buffer,&bptr,10); // read number of global rows
     int j = strtol(bptr,&bptr,10);
     if (j != nproc) ok = 0;
     else            ok = numeq_total;
     fgets(buffer,199,fp);
  }
  comm.Broadcast(&ok,1,0);
  if (!ok) return 0;
  else numeq_total = ok;
  
  int* gupdate = new int[numeq_total];
  if (proc==0)
  {
     for (int i=0; i<numeq_total; i++)
     {
        int row = strtol(buffer,&bptr,10);
        int thisproc = strtol(bptr,&bptr,10);
        gupdate[row] = thisproc;
        fgets(buffer,199,fp);
     }
     fclose(fp); fp = 0;
  }
   
  comm.Broadcast(gupdate,numeq_total,0);
  for (int i=0; i< numeq_total; i++)
     if (gupdate[i]==proc) numeq++;
     
  int* update = new int[numeq];
  
  int counter=0;
  for (int i=0; i<numeq_total; i++)
  {
     if (gupdate[i]==proc)
     {
        update[counter] = i;
        ++counter;
     }
  }   
  delete [] gupdate; gupdate = 0;
  
  map = new Epetra_Map(numeq_total,numeq,update,0,comm);
  
  return map;
}

/*----------------------------------------------------------------------*
 |                                                           m.gee 03/05|
 |                                                                      |
 | reads a matrix in aztec format                                       |
 |                                                                      |
 | the file has the following format:                                   |
 | first line: <number_global_rows>                                     |
 | following lines:                                                     |
 | <globalrownumb> <val_main_diag> <globalcolnumb> <value> ... -1       |
 |                                                                      |
 | the file has to have <number_global_rows> + 1 rows                   |
 |                                                                      |
 | Input:  char* filename        name of file to read from              |
 |         Epetra_Map& dmap       a valid Epetra_Map used as RowMap     |
 |         Epetra_Comm& comm     a valid Epetra_Comm                    |
 |                                                                      |
 | Output: Epetra_CrsMatrix*      an allocated Epetra_CrsMatrix class   |
 |                               (the calling user is responsible       |
 |                                for destroying this)                  |
 |                                                                      |
 | returns Epetra_CrsMatrix* on success, NULL otherwise                 |
 |                                                                      |
 *----------------------------------------------------------------------*/
Epetra_CrsMatrix* Epetra_ML_readaztecmatrix(char* filename,Epetra_Map& map,Epetra_Comm& comm)
{
   char  buffer[10000];
   char* bptr      = 0;

   int  numeq_total = map.NumGlobalElements();
   int  nproc       = comm.NumProc();
   int  proc        = comm.MyPID();
   
   Epetra_CrsMatrix* A = new Epetra_CrsMatrix(Copy,map,map,0);
   
   for (int activeproc=0; activeproc<nproc; activeproc++)
   {
      int ok=0;
      FILE* fp = 0;
      if (activeproc==proc)
      {
         cout << "Proc " << proc << " is reading the Epetra_CrsMatrix .."; fflush(stdout);
         fp = fopen(filename,"r");
         if (fp) 
         {
            ok=1;
            fgets(buffer,9999,fp);
            int readnumeq = strtol(buffer,&bptr,10);
            if (readnumeq != numeq_total)
               ok = 0;
         }
         else ok = 0;
      }
      comm.Broadcast(&ok,1,activeproc);
      if (!ok)
      {
         delete A;
         return 0;
      }
      if (activeproc==proc)
      {
         for (int i=0; i<numeq_total; i++)
         {
            fgets(buffer,9999,fp);
            int row = i;
            if (!map.MyGID(row)) // it's not one of my rows 
               continue;
            else                 // this row belongs to me, read it
            {
               cout << "."; fflush(stdout);
               // read row and insert them
               bptr = buffer;
               int column = 0;
               while (column != -1)
               {
                  column = strtol(bptr,&bptr,10);
                  if (column == -1) break;
                  double value = strtod(bptr,&bptr);
                  A->InsertGlobalValues(row,1,&value,&column);
               }
            }
         }
         cout << endl;   
         fclose(fp); fp = 0;
      }
      comm.Barrier();
   }
   A->FillComplete();
   
   return A;
}



/*----------------------------------------------------------------------*
 |                                                           m.gee 03/05|
 |                                                                      |
 | reads a vector in aztec format                                       |
 |                                                                      |
 | the file has the following format:                                   |
 | first line: <number_global_rows>                                     |
 | following lines:                                                     |
 | <globalrownumb> <val_> -1                                            |
 |                                                                      |
 | the file has to have <number_global_rows> + 1 rows                   |
 |                                                                      |
 | Input:  char* filename             name of file to read from         |
 | Output: Epetra_MultiVector& Vector valid Epetra_MultiVector          |
 |                                     matching the map                 |
 | Input:  Epetra_Map& map            a valid Epetra_Map                |
 |         Epetra_Comm& comm          a valid Epetra_Comm               |
 |         int ivec                   indice of vector to put values in |
 |                                                                      |
 | returns true on success, false otherwise                             |
 |                                                                      |
 *----------------------------------------------------------------------*/
bool Epetra_ML_readaztecvector(char* filename, Epetra_MultiVector& Vector, 
                               Epetra_Map& map,Epetra_Comm& comm, int ivec)
{
  char  buffer[200];
  char* bptr      = 0;

  int  numeq_total = map.NumGlobalElements();
  int  nproc       = comm.NumProc();
  int  proc        = comm.MyPID();
   
  FILE *fp = fopen(filename,"r");
  if (!fp) return false;
  if (proc) 
  {
    fclose(fp);
    fp = 0;
  }

  int ok = 0;
  if (proc==0)
  {
     fgets(buffer,199,fp);
     int tmp = strtol(buffer,&bptr,10); // read number of global rows
     if (tmp != numeq_total) ok = 0;
     else                    ok = 1;
     fclose(fp); fp = 0;
  }
  comm.Broadcast(&ok,1,0);
  if (!ok) return false;

  for (int activeproc=0; activeproc<nproc; activeproc++)
  {
     int ok = 0;
     FILE* fp = 0;
     if (activeproc==proc)
     {
        fp = fopen(filename,"r");
        if (fp)
        {
           ok = 1;
           fgets(buffer,199,fp);
        }
        else ok = 0;
     }
     comm.Broadcast(&ok,1,activeproc);
     if (!ok)
        return false;
     if (activeproc==proc)
     {
        for (int i=0; i<numeq_total; i++)
        {
           fgets(buffer,199,fp);
           int row = strtol(buffer,&bptr,10);
           if (!map.MyGID(row))
              continue;
           else
           {
              double value = strtod(bptr,&bptr);
              Vector.ReplaceGlobalValue(row,ivec,value);
           }
        }
        fclose(fp); fp = 0;
     }
     comm.Barrier();
  }
  

  return true;
}                               


/*----------------------------------------------------------------------*
 |                                                           m.gee 03/05|
 |                                                                      |
 | reads variable block information                                     |
 |                                                                      |
 | the file has the following format:                                   |
 | first line: <number_global_blocks>                                   |
 | following lines:                                                     |
 | <blocksize> <globalrownumber1> <pde_number1>  ...                    |
 |                                                                      |
 | the file has to have <number_global_blocks> + 1 rows                 |
 |                                                                      |
 | Input:  char* filename             name of file to read from         |
 |         Epetra_Map& map            a valid Epetra_Map                |
 |         Epetra_Comm& comm          a valid Epetra_Comm               |
 | Output  int** blocks               *blocks points to allocated       |
 |                                    vector matching map holding       |
 |                                    global block indices              |
 |         int** block_pde            *block_pde points to allocated    |
 |                                    vector matching map holding       |
 |                                    number of pde equation each entry |
 |                                    in *blocks belongs to             |
 |                                                                      |
 |                                                                      |
 | WARNING: The routine expects the map not to cut a single block onto  |
 |          several processors! It will return false if otherwise       |
 |                                                                      |
 | Returns true and allocated *blocks / *block_pde on success,          |
 | returns false and *blocks=NULL / *block_pde=NULL otherwise           |
 |                                                                      |
 *----------------------------------------------------------------------*/
bool Epetra_ML_readvariableblocks(char* filename, Epetra_Map& map,
                                  Epetra_Comm& comm, 
                                  int**blocks, int** block_pde)
{
  char  buffer[1000];
  char* bptr      = 0;

  int  numeq       = map.NumMyElements();
  int  nproc       = comm.NumProc();
  int  proc        = comm.MyPID();

  FILE *fp = fopen(filename,"r");
  if (!fp) return false;
  if (proc) 
  {
    fclose(fp);
    fp = 0;
  }

  int nblocks = 0;
  if (proc==0)
  {
     fgets(buffer,199,fp);
     nblocks = strtol(buffer,&bptr,10); // read number of global blocks
     fclose(fp); fp = 0;
  }
  comm.Broadcast(&nblocks,1,0);
  if (!nblocks) return false;

  *blocks    = new int[numeq];
  *block_pde = new int[numeq];
  
  int block_counter=0;
  int numeq_counter=0;
  for (int activeproc=0; activeproc<nproc; activeproc++)
  {
     int   ok = 0;
     FILE *fp = 0;
     if (activeproc==proc)
     {
        fp = fopen(filename,"r");
        if (fp)
        {
           ok = 1;
           fgets(buffer,999,fp);
        }
        else ok = 0;
     }
     comm.Broadcast(&ok,1,activeproc);
     if (!ok)
     {
        delete [] *blocks;    *blocks = 0;
        delete [] *block_pde; *block_pde = 0;
        return false;
     }
     ok = 1;
     if (activeproc==proc)
     {
        for (int i=0; i<nblocks; i++)
        {
           fgets(buffer,199,fp);
           int blocksize = strtol(buffer,&bptr,10);
           if (!blocksize)
           {
              ok = 0;
              break;
           }
           int myblock = 0;
           for (int j=0; j<blocksize; j++)
           {
              int row = strtol(bptr,&bptr,10);
              int pde = strtol(bptr,&bptr,10);
              if (map.MyGID(row)==true)
              {
                 ++myblock;
                 (*blocks)[numeq_counter]    = block_counter;
                 (*block_pde)[numeq_counter] = pde;
                 ++numeq_counter;
              }
              else if (j==0 && map.MyGID(row)==false)
                 break;
              else if (j>0 && map.MyGID(row)==false)
              {
                 cout << "**ERR** block split among several procs, abort reading\n";
                 ok = 0;
                 break;
              }
           }
           if (myblock) ++block_counter;
        if (!ok) break;
        }
        cout << "numeq " << numeq << endl;
        cout << "numeq_counter " << numeq_counter << endl;
     }
     comm.Broadcast(&ok,1,activeproc);
     if (!ok)
     {
        delete [] *blocks;    *blocks = 0;
        delete [] *block_pde; *block_pde = 0;
        return false;
     }
     comm.Broadcast(&block_counter,1,activeproc);
  }
  
  if (nblocks != block_counter)
  {
     cout << "**ERR**  Something went wrong, final number of blocks: " << block_counter << endl
          << "**ERR** not equal number of blocks from head of file : " << nblocks << endl;
     throw -1;
  }

  return true;
}                                  

/*----------------------------------------------------------------------*
 |                                                           m.gee 03/05|
 |                                                                      |
 | writes column of Epetra_Multivecotr to GID viz                       |
 |                                                                      |
 *----------------------------------------------------------------------*/
bool Epetra_ML_writegidviz(char* filename, int label, 
                           Epetra_MultiVector& vector, int ivec, 
                           Epetra_Map& map, Epetra_Comm& comm)
{
  char* bptr;
  char buffer[1000];
  char filename2[1000];
  
  int  numeq_total = map.NumGlobalElements();
  int  numeq       = map.NumMyElements();
  int  proc        = comm.MyPID();

  //----------------- reduce content of ivec Vector in vector to proc 0    
  double* values  = vector[ivec];
  double* send    = new double[numeq_total];
  double* gvalues = new double[numeq_total];
  for (int i=0; i<numeq_total; i++) send[i] = 0.0;
  for (int i=0; i<numeq; i++) 
  {
     int gID = map.GID(i);
     if (gID==-1) {
        cout << "**ERR Cannot find GID\n"; throw -1; }
     send[gID] = values[i];
  }
  comm.SumAll(send,gvalues,numeq_total);
  delete [] send;
  if (proc) delete [] gvalues;
  
  // ---------------------------------------------------open all files
  // copy filename not to modify it 
  strcpy(filename2,filename);
  int   ok    = 0;
  FILE* fin   = 0;
  FILE* foutr = 0;
  FILE* foutm = 0;
  if (proc==0)
  {
     fin = fopen(filename2,"r");
     if (fin) ok = 1;
  }
  comm.Broadcast(&ok,1,0);
  if (!ok)
  {
     delete [] gvalues;
     return false;
  } 
  bool newresfile=true;
  if (proc==0)
  {
     // try to open the mesh file for read to see whether it exists
     foutm = fopen("data.flavia.msh","r");
     if (foutm) // mesh exists, don't have to recreate
     {
        fclose(foutm); foutm = 0;
     }
     else // mesh file does not exist, create      
        foutm = fopen("data.flavia.msh","w");
     
     // try to open the mesh file for read to see whether it exists
     foutr = fopen("data.flavia.res","r");
     if (foutr) // result file exists, attach to it
     {
        fclose(foutr);
        foutr = fopen("data.flavia.res","a+w");
        newresfile=false;
     }
     else // result file does nopt exist yet, create it
     {
        foutr = fopen("data.flavia.res","w");
        newresfile=true;
     }
  }
  
  
  //----------------------------------- read the grid file
  int nnode      = 0;
  int dofpernode = 0; 
  int readnumeq  = 0;
  bool isshell=false;
  if (proc==0)
  {
     // read the grid file
     fgets(buffer,999,fin);
     while (strpbrk(buffer,"#"))
        fgets(buffer,999,fin);
     nnode      = strtol(buffer,&bptr,10);
     dofpernode = strtol(bptr,&bptr,10);
     readnumeq  = strtol(bptr,&bptr,10);
     if (strncmp(" SHELL",bptr,6)==0) 
        isshell=true;
     else                            
        isshell=false;
     if (readnumeq==numeq_total) ok=1;
     else                        ok=0;  
  }
  comm.Broadcast(&ok,1,0);
  if (!ok)
  {
     delete [] gvalues;
     return false;
  }
  
  //-------------------------- read nodal coordinates and dofs
  double* x   = 0;
  double* y   = 0;
  double* z   = 0;
  int**   dof = 0;
  if (proc==0)
  {
     // allocate vectors for nodal coordinates
     x = new double[nnode];
     y = new double[nnode];
     z = new double[nnode];
     // create array for dofs on each node
     dof    = new int* [nnode];
     dof[0] = new int[nnode*dofpernode];
     for (int i=1; i<nnode; i++)
        dof[i] = &(dof[0][i*dofpernode]);
     // read the nodes
     for (int i=0; i<nnode; i++)
     {
        fgets(buffer,999,fin);
        int node    = strtol(buffer,&bptr,10);
        x[node-1]   = strtod(bptr,&bptr);
        y[node-1]   = strtod(bptr,&bptr);
        z[node-1]   = strtod(bptr,&bptr);
        for (int j=0; j<dofpernode; j++)
           dof[node-1][j] = strtol(bptr,&bptr,10); 
     }
     // check whether we arrived at the line, were the elements begin
     fgets(buffer,999,fin);
     if (!(strpbrk(buffer,"#")))
     {
        ok = 0;
        delete [] x; delete [] y; delete [] z;
        delete [] dof[0]; delete [] dof;
     }
     else ok = 1;
  }     
  comm.Broadcast(&ok,1,0);
  if (!ok)
  {
     delete [] gvalues;
     return false;
  }

  //---------------------------- read the element topology
  int   nelement    = 0;
  int   nodesperele = 0;
  int** top         = 0;
  if (proc==0)
  {
     // read the elements
     fgets(buffer,999,fin);
     nelement    = strtol(buffer,&bptr,10);
     nodesperele = strtol(bptr,&bptr,10);
     // allocate array for topology
     top    = new int* [nelement];
     top[0] = new int[nelement*nodesperele];
     for (int i=1; i<nelement; i++)
        top[i] = &(top[0][i*nodesperele]);
     // read the elements
     for (int i=0; i<nelement; i++)
     {
        fgets(buffer,999,fin);
        int element    = strtol(buffer,&bptr,10);
        for (int j=0; j<nodesperele; j++)
          top[element-1][j] = strtol(bptr,&bptr,10);
     }
     // check for end of elements marker
     fgets(buffer,999,fin);
     if (!(strpbrk(buffer,"#")))
     {
        ok = 0;
        delete [] x; delete [] y; delete [] z;
        delete [] dof[0]; delete [] dof;
        delete [] top[0]; delete [] top;        
     }
     else ok = 1;
     fclose(fin);   fin   = 0;
  }  
  comm.Broadcast(&ok,1,0);
  if (!ok)
  {
     delete [] gvalues;
     return false;
  }

  //------------------------- printf the .flavia.msh file
  if (proc==0 && foutm)
  {
     // print nodal coordinates
     fprintf(foutm,"#-------------------------------------------------------------------------------\n");
     fprintf(foutm,"# visualization using GID\n");
     fprintf(foutm,"#-------------------------------------------------------------------------------\n");
     fprintf(foutm,"MESH datamesh DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n");
     fprintf(foutm,"COORDINATES\n");
     for (int i=0; i<nnode; i++)
     fprintf(foutm,"%6d   %20.10f   %20.10f   %20.10f\n",i+1,x[i],y[i],z[i]);
     fprintf(foutm,"END COORDINATES\n");
     
     // print elements
     fprintf(foutm,"ELEMENTS\n");
     for (int i=0; i<nelement; i++)
     {
        fprintf(foutm,"%6d   ",i+1);
        for (int j=0; j<nodesperele; j++)
           fprintf(foutm,"%6d   ",top[i][j]);
        fprintf(foutm,"\n");
     }
     fprintf(foutm,"END ELEMENTS\n");
     fflush(foutm);
     fclose(foutm); foutm = 0;
  }

  //------------- printf the .flavia.res file with the vector
  if (proc==0)
  {
     char sign='"';
     if (newresfile)
     {
     fprintf(foutr,"Gid Post Results File 1.0\n");
     fprintf(foutr,"#-------------------------------------------------------------------------------\n");
     fprintf(foutr,"# visualization using GID\n");
     fprintf(foutr,"#-------------------------------------------------------------------------------\n");
     fprintf(foutr,"RESULTRANGESTABLE %cstandard%c\n",sign,sign);
     fprintf(foutr,"            - -1000000.0 : %cvery small%c\n",sign,sign);
     fprintf(foutr," -1000000.0 -  1000000.0 : %cnormal%c\n",sign,sign);
     fprintf(foutr,"  1000000.0 -            : %cvery large%c\n",sign,sign);
     fprintf(foutr,"END RESULTRANGESTABLE\n");
     fprintf(foutr,"#-------------------------------------------------------------------------------\n");
     fprintf(foutr,"GAUSSPOINTS %cdatamesh%c ELEMTYPE Hexahedra %cdatamesh%c\n",sign,sign,sign,sign);
     fprintf(foutr,"NUMBER OF GAUSS POINTS: 8\n");
     fprintf(foutr,"NATURAL COORDINATES: Internal\n");
     fprintf(foutr,"END GAUSSPOINTS\n");
     fprintf(foutr,"#-------------------------------------------------------------------------------\n");
     }
     fprintf(foutr,"#===============================================================================\n");
     fprintf(foutr,"#===============================================================================\n");
     fprintf(foutr,"RESULT %cdisplacement%c %cML%c %d VECTOR ONNODES\n",sign,sign,sign,sign,label);
     fprintf(foutr,"RESULTRANGESTABLE %cstandard%c\n",sign,sign);
     fprintf(foutr,"COMPONENTNAMES %cx-displ%c,%cy-displ%c,%cz-displ%c\n",sign,sign,sign,sign,sign,sign);
     fprintf(foutr,"VALUES\n"); fflush(foutr);
     if (!isshell) // result does not come from a shell element
     {
        for (int i=0; i<nnode; i++)
        {
           fprintf(foutr," %6d   ",i+1);
           for (int j=0; j<dofpernode; j++)
           {
              int thisdof = dof[i][j];
              double val  = 0.0;
              if (thisdof<numeq_total)
                 val = gvalues[thisdof];
              else
                 val = 0.0;
              fprintf(foutr,"%20.10e   ",val); 
           }
           fprintf(foutr,"\n"); 
        }
     }
     else // results come from a shell element
     {
        int realnnode = nnode/2;
        for (int i=0; i<realnnode; i++)
        {
           int node_lower = i;
           int node_upper = i+realnnode;
           // print the lower surface node
           fprintf(foutr," %6d   ",node_lower+1);
           for (int j=0; j<dofpernode; j++)
           {
              int thisdof = dof[node_lower][j];
              double val = 0.0;
              if (thisdof<numeq_total)
                 val = gvalues[thisdof];
              else
                 val = 0.0;
              // this is mid surface displacement, subtract the relativ displacement
              int reldof  = dof[node_upper][j];
              double val2 = 0.0;
              if (reldof<numeq_total)
                 val2 = gvalues[reldof];
              else
                 val2 = 0.0;
              val -= val2;
              fprintf(foutr,"%20.10e   ",val); 
           }
           fprintf(foutr,"\n"); fflush(foutr);
           // print the upper surface node
           fprintf(foutr," %6d   ",node_upper+1);
           for (int j=0; j<dofpernode; j++)
           {
              int thisdof = dof[node_upper][j];
              double val = 0.0;
              if (thisdof<numeq_total)
                 val = gvalues[thisdof];
              else
                 val = 0.0;
              // this is a relativ displcement, add mid surface displacement to get upper total displ.
              int middof  = dof[node_lower][j];
              double val2 = 0.0;
              if (middof<numeq_total)
                 val2 = gvalues[middof];
              else
                 val2 = 0.0;
              val += val2;
              fprintf(foutr,"%20.10e   ",val); 
           }
           fprintf(foutr,"\n"); fflush(foutr);
        }
     }
     fprintf(foutr,"END VALUES\n");
  }
  // clean up
  if (proc==0)
  {
     delete [] x; delete [] y; delete [] z;
     delete [] dof[0]; delete [] dof;
     delete [] top[0]; delete [] top;        
     delete [] gvalues; 
     fflush(foutr); fclose(foutr);
  }
  return true;
}                           

#else

  /*noop for certain compilers*/
  int ML_EPETRA_EMPTY;

#endif /*ifdef ML_WITH_EPETRA*/

