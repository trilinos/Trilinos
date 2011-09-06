
/*
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
*/

#include "Aztec2Petra.h"

int Aztec2Petra(int * proc_config,
		AZ_MATRIX * Amat, double * az_x, double * az_b,
		Epetra_Comm * & comm,
		Epetra_BlockMap * & map,
		Epetra_RowMatrix * &A,
		Epetra_Vector * & x,
		Epetra_Vector * & b,
		int ** global_indices) {

  // Build Epetra_Comm object

#ifdef AZTEC_MPI
    MPI_Comm * mpicomm = (MPI_Comm * ) AZ_get_comm(proc_config);
    comm = (Epetra_Comm *) new Epetra_MpiComm(*mpicomm);
#else
    comm = (Epetra_Comm *) new Epetra_SerialComm();
#endif  

  int * MyGlobalElements, *global_bindx, *update;
  
  if (!Amat->has_global_indices) {
    //create a global bindx
    AZ_revert_to_global(proc_config, Amat, &global_bindx, &update);
    MyGlobalElements = update;
  }
  else // Already have global ordering
    {
      global_bindx = Amat->bindx;
      MyGlobalElements = Amat->update;
      if (MyGlobalElements==0) EPETRA_CHK_ERR(-1);
    }

  // Get matrix information
  int NumMyElements = 0;
  if (Amat->data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
    NumMyElements = Amat->data_org[AZ_N_int_blk] + Amat->data_org[AZ_N_bord_blk];
  else
    NumMyElements = Amat->data_org[AZ_N_internal] + Amat->data_org[AZ_N_border];
  // int NumMyElements = Amat->N_update; // Note: This "official" way does not always work
  int * bpntr = Amat->bpntr;
  int * rpntr = Amat->rpntr;
  int * indx = Amat->indx;
  double * val = Amat->val;

  int NumGlobalElements;
  comm->SumAll(&NumMyElements, &NumGlobalElements, 1);


  // Make ElementSizeList (if VBR) - number of block entries in each block row

  int * ElementSizeList = 0;

  if (Amat->data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
  
    ElementSizeList = new int[NumMyElements];
    if (ElementSizeList==0) EPETRA_CHK_ERR(-1); // Ran out of memory
    
    for (int i=0; i<NumMyElements; i++) ElementSizeList[i] = rpntr[i+1] - rpntr[i];

    map = new Epetra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, 
			     ElementSizeList, 0, *comm);

    if (map==0) EPETRA_CHK_ERR(-2); // Ran out of memory

    delete [] ElementSizeList;
 
    Epetra_VbrMatrix * AA = new Epetra_VbrMatrix(View, *map, 0);
  
    if (AA==0) EPETRA_CHK_ERR(-3); // Ran out of memory

    /* Add block rows one-at-a-time */
    {for (int i=0; i<NumMyElements; i++) {
      int BlockRow = MyGlobalElements[i];
      int NumBlockEntries = bpntr[i+1] - bpntr[i];
      int *BlockIndices = global_bindx + bpntr[i];
      int ierr = AA->BeginInsertGlobalValues(BlockRow, NumBlockEntries, BlockIndices);
      if (ierr!=0) {
	cerr << "Error in BeginInsertGlobalValues(GlobalBlockRow = " << BlockRow 
	     << ") = " << ierr << endl; 
	EPETRA_CHK_ERR(ierr);
      }
      int LDA = rpntr[i+1] - rpntr[i];
      int NumRows = LDA;
      for (int j=bpntr[i]; j<bpntr[i+1]; j++) {
	int NumCols = (indx[j+1] - indx[j])/LDA;
	double * Values = val + indx[j];
	ierr = AA->SubmitBlockEntry(Values, LDA, NumRows, NumCols);
	if (ierr!=0) {
	  cerr << "Error in SubmitBlockEntry, GlobalBlockRow = " << BlockRow 
	       << "GlobalBlockCol = " << BlockIndices[j] << "Error = " << ierr << endl; 
	  EPETRA_CHK_ERR(ierr);
	}
      }
      ierr = AA->EndSubmitEntries();
      if (ierr!=0) {
	cerr << "Error in EndSubmitEntries(GlobalBlockRow = " << BlockRow 
	     << ") = " << ierr << endl; 
	EPETRA_CHK_ERR(ierr);
      }
    }}  
    int ierr=AA->FillComplete();    
    if (ierr!=0) {
      cerr <<"Error in Epetra_VbrMatrix FillComplete" << ierr << endl;
      EPETRA_CHK_ERR(ierr);
    }
    
    A = dynamic_cast<Epetra_RowMatrix *> (AA); // cast VBR pointer to RowMatrix pointer
  }
  else if  (Amat->data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
  
    /* Make numNzBlks - number of block entries in each block row */

    int * numNz = new int[NumMyElements];
    for (int i=0; i<NumMyElements; i++) numNz[i] = global_bindx[i+1] - global_bindx[i] + 1;

    Epetra_Map * map1 = new Epetra_Map(NumGlobalElements, NumMyElements,
				     MyGlobalElements, 0, *comm);

    Epetra_CrsMatrix * AA = new Epetra_CrsMatrix(Copy, *map1, numNz);

    map = (Epetra_BlockMap *) map1; // cast Epetra_Map to Epetra_BlockMap

    /* Add  rows one-at-a-time */

    for (int row=0; row<NumMyElements; row++) {
      double * row_vals = val + global_bindx[row];
      int * col_inds = global_bindx + global_bindx[row];
      int numEntries = global_bindx[row+1] - global_bindx[row];
      int ierr = AA->InsertGlobalValues(MyGlobalElements[row], numEntries, row_vals, col_inds);
      if (ierr!=0) {
	cerr << "Error puting row " << MyGlobalElements[row] << endl;
	EPETRA_CHK_ERR(ierr);
      }
      ierr = AA->InsertGlobalValues(MyGlobalElements[row], 1, val+row, MyGlobalElements+row);
      if (ierr!=0) {
	cerr << "Error putting  diagonal" << endl;
	EPETRA_CHK_ERR(ierr);
      }
    }

    int ierr=AA->FillComplete();
    if (ierr!=0) {
      cerr << "Error in Epetra_CrsMatrix_FillComplete" << endl;
      EPETRA_CHK_ERR(ierr);
    }
    A = dynamic_cast<Epetra_RowMatrix *> (AA); // cast CRS pointer to RowMatrix pointer
  }
  else cerr << "Not a supported AZ_MATRIX data type" << endl;


  // Create x vector
  x = new Epetra_Vector(View, *map,az_x);

  
  // RPP: Can not use the OperatorRangeMap in the ctor of the "b" vector 
  // below.  In MPSalsa, we delete the VbrMatrix yet still use the vector "b".
  // Deleting the matrix deletes the OperatorRangeMap that the b vector is 
  // based on.  Losing the map means "b" and all vectors that are created 
  // with the copy constructor of "b" break.  Mike has suggested 
  // using reference counting (Boost smart pointers) so the map is not 
  // deleted.  For now we will use the "map" variable as the base map for "b". 
  //b = new Epetra_Vector (View, A->OperatorRangeMap(), az_b);
  b = new Epetra_Vector (View, *map, az_b);

  *global_indices = 0; // Assume return array will be empty
  if (!Amat->has_global_indices) {
   AZ_free((void *) update);
   if (Amat->data_org[AZ_matrix_type] != AZ_VBR_MATRIX)
     AZ_free((void *) global_bindx);
   else
     global_indices = &global_bindx;
   }
  return 0;
}
