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
/*           Marzio Sala (SNL)                                          */
/************************************************************************/

#include "ml_common.h"

#ifdef ML_WITH_EPETRA
#include <vector>
#include "ml_epetra_utils.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"

#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


// ====================================================================== 

int Epetra_ML_matvec(ML_Operator *data, int in, double *p, int out, double *ap)
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

int Epetra_ML_getrow(ML_Operator *data, int N_requested_rows, int requested_rows[], 
		    int allocated_space, int columns[], double values[],
		    int row_lengths[])
{
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
      if (ierr) return(0); 
      NumEntries = NumBlockEntries*NumPDEEqns;
      if (nz_ptr + NumEntries > allocated_space) return(0);
      
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
    if (ierr) return(0); //JJH I think this is the correct thing to return if
                         //    A->ExtractMyRowCopy returns something nonzero ..

    if( !MatrixIsVbrMatrix ) {
      row_lengths[i] = NumEntries;
      if (nz_ptr + NumEntries > allocated_space) return(0);
      
      for (int j=0; j<NumEntries; j++) {
	columns[nz_ptr] = Indices[j];
	values[nz_ptr++] = Values[j];
      }
    }
  }
  

  return(1);
}

// ====================================================================== 

int Epetra_ML_comm_wrapper(double vec[], void *data)
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

int Epetra2MLMatrix(Epetra_RowMatrix * A, ML_Operator *newMatrix)
{
  int isize, osize;

  osize = A->NumMyRows();
  isize = A->OperatorDomainMap().NumMyElements();
  //  isize = A->NumMyCols();
  int N_ghost = A->RowMatrixColMap().NumMyElements() - isize;

  if (N_ghost < 0) N_ghost = 0;  // A->NumMyCols() = 0 for an empty matrix

  ML_Operator_Set_ApplyFuncData(newMatrix, isize, osize,
                              (void*) A, osize,
                              NULL, 0);

  ML_CommInfoOP_Generate( &(newMatrix->getrow->pre_comm), 
                        Epetra_ML_comm_wrapper, (void *) A, 
                        newMatrix->comm, isize, N_ghost);

  ML_Operator_Set_Getrow(newMatrix, newMatrix->outvec_leng,
                       Epetra_ML_getrow);

  ML_Operator_Set_ApplyFunc (newMatrix, Epetra_ML_matvec);

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
    }
  }
  
  assert(Q->FillComplete(ReorderedMap,StartingMap)==0);
  
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

  for( int i=0 ; i<ReorderedNumElements ; ++i ) {
    ReorderedBdry[i] = yyy[i*NumPDEEqns];
  }
  
  ML_Q2 = ML_Operator_Create( ml_communicator );  
  
  Epetra2MLMatrix( Q, ML_Q2);

  if( Start != NULL ) delete [] Start;
  if( Reord != NULL ) delete [] Reord;

  return ML_Q2;

}


void ML_DestroyQ(void) 
{

  delete Q;
  Q = NULL;

  return;
  
} /* ML_DestroyQ */

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
  assert(Qt->FillComplete(ReorderedMap,StartingMap)==0);
  
  ML_Qt2 = ML_Operator_Create( ml_communicator );

  Epetra2MLMatrix( Qt, ML_Qt2);

  return ML_Qt2;
#endif
  
} /* ML_BuildQt */

void ML_DestroyQt( void ) 
{

  delete Qt;
  Qt = NULL;

  return;
  
} /* ML_DestroyQt */

} /* extern "C" */

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

  vector<double> global_osize(osize);
  vector<int>    global_osize_as_int(osize);
  
  for (int i = 0 ; i < isize; i++) {
          global_isize[i] = (double) (isize_offset + i);
          global_isize_as_int[i] = isize_offset + i;
  }
          
  for (int i = 0 ; i < osize; i++) {
    global_osize[i] = (double) (osize_offset + i);
    global_osize_as_int[i] = osize_offset + i;
  }
  for (int i = 0 ; i < Nghost; i++) global_isize[i+isize] = -1;
  
  Epetra_Map  rangemap( -1, osize, &global_osize_as_int[0], 0, EpetraComm ) ; 
  Epetra_Map  domainmap( -1, isize, &global_isize_as_int[0], 0, EpetraComm ) ; 
  
  CrsMatrix = new Epetra_CrsMatrix( Copy, rangemap, 0 ); 
  
  ML_exchange_bdry(&global_isize[0],Amat->getrow->pre_comm, 
 		 Amat->invec_leng,comm,ML_OVERWRITE,NULL);

  for ( int j = 0; j < isize+Nghost; j++ ) { 
    global_isize_as_int[j] = (int) global_isize[j];
  }

  // MS // introduced variable allocation for colInd and colVal
  // MS // improved efficiency in InsertGlobalValues

  int allocated = 1;
  int * colInd = new int[allocated];
  double * colVal = new double[allocated];
  int NumNonzeros;
  int ierr;
  int    ncnt;

  MaxNumNonzeros=0;
  
  for (int i = 0; i < osize; i++)
  {
    ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,&ncnt);

    if( ierr == 0 ) {
      do {
	delete [] colInd;
	delete [] colVal;
	allocated *= 2;
	colInd = new int[allocated];
	colVal = new double[allocated];
	ierr = ML_Operator_Getrow(Amat,1,&i,allocated,colInd,colVal,&ncnt);
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
    
    CrsMatrix->InsertGlobalValues( global_osize_as_int[i], NumNonzeros, 
				   colVal, colInd);
    
  }

  delete [] colInd;
  delete [] colVal;
  
  assert(CrsMatrix->FillComplete(domainmap,rangemap)==0);

  CPUTime = Time.ElapsedTime();

  return 0;
  
} /* ML_Operator2EpetraCrsMatrix for rectangular matrices*/

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

  assert(CrsMatrix->FillComplete()==0);

  CPUTime = Time.ElapsedTime();

  return 0;
  
} /* ML_Operator2EpetraCrsMatrix */

#ifdef WKC
int Epetra_ML_matvec_WKC (ML_Operator *data, int in, double *p, int out, double *ap)
{
  ML_Operator *mat_in;

  mat_in = (ML_Operator *) data;
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

extern "C" {

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
  // FIXME: I work with MPI_COMM_WORLD only, not with the comm of ML!
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
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
      colVal[j] = abs(colVal[j]);
      for (int k = 1 ; k < NumPDEEqns ; ++k) {
	colVal[j] += abs(colVal[j+k]);
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

	double d2 = pow(coord_i[0]-coord_j[0],2) + pow(coord_i[1]-coord_j[1],2) + pow(coord_i[2]-coord_j[2],2);
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
  Epetra2MLMatrix(FakeMatrix,*NewOp);
  
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

} // extern "C"

#else

  /*noop for certain compilers*/
  int ML_EPETRA_EMPTY;

#endif /*ifdef ML_WITH_EPETRA*/

