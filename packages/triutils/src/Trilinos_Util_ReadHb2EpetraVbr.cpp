#include <stdlib.h>
#include <stdio.h>
#include "Trilinos_Util.h"
#include "iohb.h"
#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"

void Trilinos_Util_ReadHb2EpetraVbr(char *data_file, char * partitioning,
				 const Epetra_Comm  &comm, 
				 Epetra_BlockMap *& map, 
				 Epetra_VbrMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact) {

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

  int NumGlobalEquations, NumMyNonzeros;
  double *val_msr, *x_in, *b_in, *xexact_in;
  int *bindx_msr;
  
  /* Set exact solution to NULL */
  xexact = NULL;
  Trilinos_Util_read_hb(data_file, comm.MyPID(), &NumGlobalEquations, &NumMyNonzeros,
			&val_msr,  &bindx_msr, &x_in, &b_in, &xexact_in);
  
  double *val;
  int NumGlobalElements, NumGlobalBlockEntries, *indx, *rpntr, *cpntr, *bpntr, *bindx;
  int NumMyBlockEntries, NumMyElements, * MyGlobalElements;
  
  Trilinos_Util_create_vbr(comm, partitioning,
			   &NumGlobalEquations, &NumGlobalElements, 
			   &NumMyNonzeros, &NumMyBlockEntries,
			   &NumMyElements, &MyGlobalElements,
			   bindx_msr, val_msr,
			   &val, &indx, &rpntr, &cpntr,
			   &bpntr, &bindx);
  
  if(comm.MyPID()==0)
    {
      free ((void *) val_msr);
      free ((void *) bindx_msr);
      free ((void *) cpntr);
    }
  
  int * ElementSizeList = 0;
  if (NumMyElements>0) ElementSizeList = new int[NumMyElements];
  
  for (int i=0; i<NumMyElements; i++) ElementSizeList[i] = rpntr[i+1] - rpntr[i];
  
  map = new Epetra_BlockMap(-1, NumMyElements, MyGlobalElements, 
			    ElementSizeList, 0, comm);
  
   A = new Epetra_VbrMatrix(Copy, *map, 0);
  
  /* Add block rows one-at-a-time */
  
  for (int i=0; i<NumMyElements; i++) {
    int BlockRow = MyGlobalElements[i];
    int NumBlockEntries = bpntr[i+1] - bpntr[i];
    int *BlockIndices = bindx + bpntr[i];
    int ierr = A->BeginInsertGlobalValues(BlockRow, NumBlockEntries, BlockIndices);
    if (ierr!=0) {
      cerr << "Error in BeginInsertGlobalValues(GlobalBlockRow = " << BlockRow 
	   << ") = " << ierr << endl; 
      abort();
    }
    int LDA = ElementSizeList[i];
    int NumRows = LDA;
    for (int j=bpntr[i]; j<bpntr[i+1]; j++) {
      int NumCols = (indx[j+1] - indx[j])/LDA;
      double * Values = val + indx[j];
      ierr = A->SubmitBlockEntry(Values, LDA, NumRows, NumCols);
      if (ierr!=0) {
	cerr << "Error in SubmitBlockEntry, GlobalBlockRow = " << BlockRow 
	     << "GlobalBlockCol = " << BlockIndices[j] << "Error = " << ierr << endl; 
	abort();
      }
    }
    ierr = A->EndSubmitEntries();
    if (ierr!=0) {
      cerr << "Error in EndSubmitEntries(GlobalBlockRow = " << BlockRow 
	   << ") = " << ierr << endl; 
      abort();
    }
  }  
  int ierr=A->TransformToLocal();    
  if (ierr!=0) cerr << "Error in Epetra_VbrMatrix TransformToLocal ierr = " << ierr << endl;
  
  xexact = new Epetra_Vector(Copy, *map, xexact_in);
  x = new Epetra_Vector(Copy, *map, x_in);
  b = new Epetra_Vector(Copy, *map, b_in);

  if(comm.MyPID()==0)
    {
      free ((void *) val);
      free ((void *) indx);
      free ((void *) rpntr);
      free ((void *) bpntr);
      free ((void *) bindx);
    }
  return;
}
