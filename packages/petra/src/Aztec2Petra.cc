#include "Aztec2Petra.h"

int Aztec2Petra(int * proc_config,
		AZ_MATRIX * Amat, double * az_x, double * az_b,
		Petra_Comm * & comm,
		Petra_BlockMap * & map,
		Petra_RDP_RowMatrix * &A,
		Petra_RDP_Vector * & x,
		Petra_RDP_Vector * & b) {

  // Build Petra_Comm object

#ifdef AZ_MPI
    MPI_Comm * mpicomm = (MPI_Comm * ) AZ_get_comm(proc_config);
    comm = new Petra_Comm(*mpicomm);
#else
    comm = new Petra_Comm();
#endif  

  int * MyGlobalElements, * global_bindx, *update;
  
  if (!Amat->has_global_indices) {
    //create a global bindx
    AZ_revert_to_global(proc_config, Amat, &global_bindx, &update);
    MyGlobalElements = update;
  }
  else // Already have global ordering
    {
      global_bindx = Amat->bindx;
      MyGlobalElements = Amat->update;
      if (MyGlobalElements==0) PETRA_CHK_ERR(-1);
    }

  // Get matrix information
  int NumMyElements = 0;
  if (Amat->data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
    NumMyElements = Amat->data_org[AZ_N_int_blk] + Amat->data_org[AZ_N_bord_blk];
  else
    NumMyElements = Amat->data_org[AZ_N_internal] + Amat->data_org[AZ_N_border];
  // int NumMyElements = Amat->N_update; // Note: This "official" way does not always work
  int * bpntr = Amat->bpntr;
  int * bindx = Amat->bindx;
  int * rpntr = Amat->rpntr;
  int * indx = Amat->indx;
  double * val = Amat->val;

  int NumGlobalElements;
  comm->SumAll(&NumMyElements, &NumGlobalElements, 1);


  // Make ElementSizeList (if VBR) - number of block entries in each block row

  int * ElementSizeList = 0;

  if (Amat->data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
  
    ElementSizeList = new int[NumMyElements];
    if (ElementSizeList==0) PETRA_CHK_ERR(-1); // Ran out of memory
    
    for (int i=0; i<NumMyElements; i++) ElementSizeList[i] = rpntr[i+1] - rpntr[i];

    map = new Petra_BlockMap(NumGlobalElements, NumMyElements, MyGlobalElements, 
			     ElementSizeList, 0, *comm);

    if (map==0) PETRA_CHK_ERR(-2); // Ran out of memory

    delete [] ElementSizeList;
 
    Petra_RDP_VBR_Matrix * AA = new Petra_RDP_VBR_Matrix(Copy, *map, 0);
  
    if (AA==0) PETRA_CHK_ERR(-3); // Ran out of memory

    /* Add block rows one-at-a-time */
    for (int i=0; i<NumMyElements; i++) {
      int BlockRow = MyGlobalElements[i];
      int NumBlockEntries = bpntr[i+1] - bpntr[i];
      int *BlockIndices = global_bindx + bpntr[i];
      int ierr = AA->BeginInsertGlobalValues(BlockRow, NumBlockEntries, BlockIndices);
      if (ierr!=0) {
	cerr << "Error in BeginInsertGlobalValues(GlobalBlockRow = " << BlockRow 
	     << ") = " << ierr << endl; 
	PETRA_CHK_ERR(ierr);
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
	  PETRA_CHK_ERR(ierr);
	}
      }
      ierr = AA->EndSubmitEntries();
      if (ierr!=0) {
	cerr << "Error in EndSubmitEntries(GlobalBlockRow = " << BlockRow 
	     << ") = " << ierr << endl; 
	PETRA_CHK_ERR(ierr);
      }
    }  
    int ierr=AA->TransformToLocal();    
    if (ierr!=0) {
      cerr <<"Error in Petra_RDP_VBR_Matrix TransformToLocal" << ierr << endl;
      PETRA_CHK_ERR(ierr);
    }

    A = (Petra_RDP_RowMatrix *) AA; // cast VBR pointer to RowMatrix pointer
  }
  else if  (Amat->data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
  
    /* Make numNzBlks - number of block entries in each block row */

    int * numNz = new int[NumMyElements];
    for (int i=0; i<NumMyElements; i++) numNz[i] = global_bindx[i+1] - global_bindx[i] + 1;

    /* Make ColInds - Exactly global_bindx, offset by diag (just copy pointer) */
    int * ColInds = global_bindx+NumMyElements+1;

    Petra_Map * map1 = new Petra_Map(NumGlobalElements, NumMyElements,
				     MyGlobalElements, 0, *comm);

    Petra_RDP_CRS_Matrix * AA = new Petra_RDP_CRS_Matrix(Copy, *map1, numNz);

    map = (Petra_BlockMap *) map1; // cast Petra_Map to Petra_BlockMap

    /* Add  rows one-at-a-time */

    for (int row=0; row<NumMyElements; row++) {
      double * row_vals = val + global_bindx[row];
      int * col_inds = global_bindx + global_bindx[row];
      int numEntries = global_bindx[row+1] - global_bindx[row];
      int ierr = AA->InsertGlobalValues(MyGlobalElements[row], numEntries, row_vals, col_inds);
      if (ierr!=0) {
	cerr << "Error puting row " << MyGlobalElements[row] << endl;
	PETRA_CHK_ERR(ierr);
      }
      ierr = AA->InsertGlobalValues(MyGlobalElements[row], 1, val+row, MyGlobalElements+row);
      if (ierr!=0) {
	cerr << "Error putting  diagonal" << endl;
	PETRA_CHK_ERR(ierr);
      }
    }

    int ierr=AA->TransformToLocal();
    if (ierr!=0) {
      cerr << "Error in Petra_RDP_CRS_Matrix_TransformToLocal" << endl;
      PETRA_CHK_ERR(ierr);
    }
    A = (Petra_RDP_RowMatrix *) AA; // cast CRS pointer to RowMatrix pointer
  }
  else cerr << "Not a supported AZ_MATRIX data type" << endl;

  // Create x vector, note that it is a "long" vector (has ghost entries).
  x = new Petra_RDP_Vector(View, A->ImportMap(),az_x);

  b = new Petra_RDP_Vector (View, *map, az_b);

  if (!Amat->has_global_indices) {
    delete global_bindx;
    delete update;
  }
  return 0;
}
