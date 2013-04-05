//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include <EpetraExt_ConfigDefs.h>
#include <EpetraExt_MatrixMatrix.h>

#include <EpetraExt_MMHelpers.h>

#include <EpetraExt_Transpose_RowMatrix.h>

#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Util.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Directory.h>
#include <Epetra_HashTable.h>
#include <Epetra_Distributor.h>


#include <Teuchos_TimeMonitor.hpp>

#ifdef HAVE_VECTOR
#include <vector>
#endif



namespace EpetraExt {
//=========================================================================
//
//Method for internal use... sparsedot forms a dot-product between two
//sparsely-populated 'vectors'.
//Important assumption: assumes the indices in u_ind and v_ind are sorted.
//
double sparsedot(double* u, int* u_ind, int u_len,
		 double* v, int* v_ind, int v_len)
{
  double result = 0.0;

  int v_idx = 0;
  int u_idx = 0;

  while(v_idx < v_len && u_idx < u_len) {
    int ui = u_ind[u_idx];
    int vi = v_ind[v_idx];

    if (ui < vi) {
      ++u_idx;
    }
    else if (ui > vi) {
      ++v_idx;
    }
    else {
      result += u[u_idx++]*v[v_idx++];
    }
  }

  return(result);
}

//=========================================================================
//kernel method for computing the local portion of C = A*B^T
int mult_A_Btrans(CrsMatrixStruct& Aview,
		  CrsMatrixStruct& Bview,
		  CrsWrapper& C)
{
  int i, j, k;
  int returnValue = 0;

  int maxlen = 0;
  for(i=0; i<Aview.numRows; ++i) {
    if (Aview.numEntriesPerRow[i] > maxlen) maxlen = Aview.numEntriesPerRow[i];
  }
  for(i=0; i<Bview.numRows; ++i) {
    if (Bview.numEntriesPerRow[i] > maxlen) maxlen = Bview.numEntriesPerRow[i];
  }

  //cout << "Aview: " << endl;
  //dumpCrsMatrixStruct(Aview);

  //cout << "Bview: " << endl;
  //dumpCrsMatrixStruct(Bview);

  int numBcols = Bview.colMap->NumMyElements();
  int numBrows = Bview.numRows;

  int iworklen = maxlen*2 + numBcols;
  int* iwork = new int[iworklen];

  int* bcols = iwork+maxlen*2;
  int* bgids = Bview.colMap->MyGlobalElements();
  double* bvals = new double[maxlen*2];
  double* avals = bvals+maxlen;

  int max_all_b = Bview.colMap->MaxAllGID();
  int min_all_b = Bview.colMap->MinAllGID();

  //bcols will hold the GIDs from B's column-map for fast access
  //during the computations below
  for(i=0; i<numBcols; ++i) {
    int blid = Bview.colMap->LID(bgids[i]);
    bcols[blid] = bgids[i];
  }

  //next create arrays indicating the first and last column-index in
  //each row of B, so that we can know when to skip certain rows below.
  //This will provide a large performance gain for banded matrices, and
  //a somewhat smaller gain for *most* other matrices.
  int* b_firstcol = new int[2*numBrows];
  int* b_lastcol = b_firstcol+numBrows;
  int temp;
  for(i=0; i<numBrows; ++i) {
    b_firstcol[i] = max_all_b;
    b_lastcol[i] = min_all_b;

    int Blen_i = Bview.numEntriesPerRow[i];
    if (Blen_i < 1) continue;
    int* Bindices_i = Bview.indices[i];

    if (Bview.remote[i]) {
      for(k=0; k<Blen_i; ++k) {
        temp = Bview.importColMap->GID(Bindices_i[k]);
        if (temp < b_firstcol[i]) b_firstcol[i] = temp;
        if (temp > b_lastcol[i]) b_lastcol[i] = temp;
      }
    }
    else {
      for(k=0; k<Blen_i; ++k) {
        temp = bcols[Bindices_i[k]];
        if (temp < b_firstcol[i]) b_firstcol[i] = temp;
        if (temp > b_lastcol[i]) b_lastcol[i] = temp;
      }
    }
  }

  Epetra_Util util;

  int* Aind = iwork;
  int* Bind = iwork+maxlen;

  bool C_filled = C.Filled();

  //To form C = A*B^T, we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(i,k)*B(j,k) )
  //
  //This is the easiest case of all to code (easier than A*B, A^T*B, A^T*B^T).
  //But it requires the use of a 'sparsedot' function (we're simply forming
  //dot-products with row A_i and row B_j for all i and j).

  //loop over the rows of A.
  for(i=0; i<Aview.numRows; ++i) {
    if (Aview.remote[i]) {
      continue;
    }

    int* Aindices_i = Aview.indices[i];
    double* Aval_i  = Aview.values[i];
    int A_len_i = Aview.numEntriesPerRow[i];
    if (A_len_i < 1) {
      continue;
    }

    for(k=0; k<A_len_i; ++k) {
      Aind[k] = Aview.colMap->GID(Aindices_i[k]);
      avals[k] = Aval_i[k];
    }

    util.Sort(true, A_len_i, Aind, 1, &avals, 0, NULL);

    int mina = Aind[0];
    int maxa = Aind[A_len_i-1];

    if (mina > max_all_b || maxa < min_all_b) {
      continue;
    }

    int global_row = Aview.rowMap->GID(i);

    //loop over the rows of B and form results C_ij = dot(A(i,:),B(j,:))
    for(j=0; j<Bview.numRows; ++j) {
      if (b_firstcol[j] > maxa || b_lastcol[j] < mina) {
        continue;
      }

      int* Bindices_j = Bview.indices[j];
      int B_len_j = Bview.numEntriesPerRow[j];
      if (B_len_j < 1) {
        continue;
      }

      int tmp, Blen = 0;

      if (Bview.remote[j]) {
        for(k=0; k<B_len_j; ++k) {
	  tmp = Bview.importColMap->GID(Bindices_j[k]);
          if (tmp < mina || tmp > maxa) {
            continue;
          }

          bvals[Blen] = Bview.values[j][k];
          Bind[Blen++] = tmp;
	}
      }
      else {
        for(k=0; k<B_len_j; ++k) {
	  tmp = bcols[Bindices_j[k]];
          if (tmp < mina || tmp > maxa) {
            continue;
          }

          bvals[Blen] = Bview.values[j][k];
          Bind[Blen++] = tmp;
	}
      }

      if (Blen < 1) {
        continue;
      }

      util.Sort(true, Blen, Bind, 1, &bvals, 0, NULL);

      double C_ij = sparsedot(avals, Aind, A_len_i,
			      bvals, Bind, Blen);

      if (C_ij == 0.0) {
	continue;
      }
      int global_col = Bview.rowMap->GID(j);

      int err = C_filled ?
        C.SumIntoGlobalValues(global_row, 1, &C_ij, &global_col)
        :
        C.InsertGlobalValues(global_row, 1, &C_ij, &global_col);

      if (err < 0) {
	return(err);
      }
      if (err > 0) {
	if (C_filled) {
	  //C.Filled()==true, and C doesn't have all the necessary nonzero
          //locations, or global_row or global_col is out of range (less
          //than 0 or non local).
          std::cerr << "EpetraExt::MatrixMatrix::Multiply Warning: failed "
              << "to insert value in result matrix at position "<<global_row
             <<"," <<global_col<<", possibly because result matrix has a "
             << "column-map that doesn't include column "<<global_col
             <<" on this proc." <<std::endl;
	  return(err);
	}
      }
    }
  }

  delete [] iwork;
  delete [] bvals;
  delete [] b_firstcol;

  return(returnValue);
}

//=========================================================================
//kernel method for computing the local portion of C = A^T*B
int mult_Atrans_B(CrsMatrixStruct& Aview,
		  CrsMatrixStruct& Bview,
		  CrsWrapper& C)
{
  int C_firstCol = Bview.colMap->MinLID();
  int C_lastCol = Bview.colMap->MaxLID();

  int C_firstCol_import = 0;
  int C_lastCol_import = -1;

  if (Bview.importColMap != NULL) {
    C_firstCol_import = Bview.importColMap->MinLID();
    C_lastCol_import = Bview.importColMap->MaxLID();
  }

  int C_numCols = C_lastCol - C_firstCol + 1;
  int C_numCols_import = C_lastCol_import - C_firstCol_import + 1;

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  double* C_row_i = new double[C_numCols];
  int* C_colInds = new int[C_numCols];

  int i, j, k;

  for(j=0; j<C_numCols; ++j) {
    C_row_i[j] = 0.0;
    C_colInds[j] = 0;
  }

  //To form C = A^T*B, compute a series of outer-product updates.
  //
  // for (ith column of A^T) { 
  //   C_i = outer product of A^T(:,i) and B(i,:)
  // Where C_i is the ith matrix update,
  //       A^T(:,i) is the ith column of A^T, and
  //       B(i,:) is the ith row of B.
  //

  //dumpCrsMatrixStruct(Aview);
  //dumpCrsMatrixStruct(Bview);
  int localProc = Bview.colMap->Comm().MyPID();

  int* Arows = Aview.rowMap->MyGlobalElements();

  bool C_filled = C.Filled();

  //loop over the rows of A (which are the columns of A^T).
  for(i=0; i<Aview.numRows; ++i) {

    int* Aindices_i = Aview.indices[i];
    double* Aval_i  = Aview.values[i];

    //we'll need to get the row of B corresponding to Arows[i],
    //where Arows[i] is the GID of A's ith row.
    int Bi = Bview.rowMap->LID(Arows[i]);
    if (Bi<0) {
      cout << "mult_Atrans_B ERROR, proc "<<localProc<<" needs row "
	   <<Arows[i]<<" of matrix B, but doesn't have it."<<endl;
      return(-1);
    }

    int* Bcol_inds = Bview.indices[Bi];
    double* Bvals_i = Bview.values[Bi];

    //for each column-index Aj in the i-th row of A, we'll update
    //global-row GID(Aj) of the result matrix C. In that row of C,
    //we'll update column-indices given by the column-indices in the
    //ith row of B that we're now holding (Bcol_inds).

    //First create a list of GIDs for the column-indices
    //that we'll be updating.

    int Blen = Bview.numEntriesPerRow[Bi];
    if (Bview.remote[Bi]) {
      for(j=0; j<Blen; ++j) {
        C_colInds[j] = Bview.importColMap->GID(Bcol_inds[j]);
      }
    }
    else {
      for(j=0; j<Blen; ++j) {
        C_colInds[j] = Bview.colMap->GID(Bcol_inds[j]);
      }
    }

    //loop across the i-th row of A (column of A^T)
    for(j=0; j<Aview.numEntriesPerRow[i]; ++j) {

      int Aj = Aindices_i[j];
      double Aval = Aval_i[j];

      int global_row;
      if (Aview.remote[i]) {
	global_row = Aview.importColMap->GID(Aj);
      }
      else {
	global_row = Aview.colMap->GID(Aj);
      }

      if (!C.RowMap().MyGID(global_row)) {
        continue;
      }

      for(k=0; k<Blen; ++k) {
        C_row_i[k] = Aval*Bvals_i[k];
      }

      //
      //Now add this row-update to C.
      //

      int err = C_filled ?
        C.SumIntoGlobalValues(global_row, Blen, C_row_i, C_colInds)
        :
        C.InsertGlobalValues(global_row, Blen, C_row_i, C_colInds);

      if (err < 0) {
        return(err);
      }
      if (err > 0) {
        if (C_filled) {
          //C is Filled, and doesn't have all the necessary nonzero locations.
          return(err);
        }
      }
    }
  }

  delete [] C_row_i;
  delete [] C_colInds;

  return(0);
}

//kernel method for computing the local portion of C = A^T*B^T
int mult_Atrans_Btrans(CrsMatrixStruct& Aview,
		       CrsMatrixStruct& Bview,
		       CrsWrapper& C)
{
  int C_firstCol = Aview.rowMap->MinLID();
  int C_lastCol = Aview.rowMap->MaxLID();

  int C_firstCol_import = 0;
  int C_lastCol_import = -1;

  if (Aview.importColMap != NULL) {
    C_firstCol_import = Aview.importColMap->MinLID();
    C_lastCol_import = Aview.importColMap->MaxLID();
  }

  int C_numCols = C_lastCol - C_firstCol + 1;
  int C_numCols_import = C_lastCol_import - C_firstCol_import + 1;

  if (C_numCols_import > C_numCols) C_numCols = C_numCols_import;

  double* dwork = new double[C_numCols];
  int* iwork = new int[C_numCols];

  double* C_col_j = dwork;
  int* C_inds = iwork;

  //cout << "Aview: " << endl;
  //dumpCrsMatrixStruct(Aview);

  //cout << "Bview: " << endl;
  //dumpCrsMatrixStruct(Bview);


  int i, j, k;

  for(j=0; j<C_numCols; ++j) {
    C_col_j[j] = 0.0;
    C_inds[j] = -1;
  }

  int* A_col_inds = Aview.colMap->MyGlobalElements();
  int* A_col_inds_import = Aview.importColMap ?
    Aview.importColMap->MyGlobalElements() : 0;

  const Epetra_Map* Crowmap = &(C.RowMap());

  //To form C = A^T*B^T, we're going to execute this expression:
  //
  // C(i,j) = sum_k( A(k,i)*B(j,k) )
  //
  //Our goal, of course, is to navigate the data in A and B once, without
  //performing searches for column-indices, etc. In other words, we avoid
  //column-wise operations like the plague...

  int* Brows = Bview.rowMap->MyGlobalElements();

  //loop over the rows of B
  for(j=0; j<Bview.numRows; ++j) {
    int* Bindices_j = Bview.indices[j];
    double* Bvals_j = Bview.values[j];

    int global_col = Brows[j];

    //loop across columns in the j-th row of B and for each corresponding
    //row in A, loop across columns and accumulate product
    //A(k,i)*B(j,k) into our partial sum quantities in C_col_j. In other
    //words, as we stride across B(j,:), we use selected rows in A to
    //calculate updates for column j of the result matrix C.

    for(k=0; k<Bview.numEntriesPerRow[j]; ++k) {

      int bk = Bindices_j[k];
      double Bval = Bvals_j[k];

      int global_k;
      if (Bview.remote[j]) {
	global_k = Bview.importColMap->GID(bk);
      }
      else {
	global_k = Bview.colMap->GID(bk);
      }

      //get the corresponding row in A
      int ak = Aview.rowMap->LID(global_k);
      if (ak<0) {
	continue;
      }

      int* Aindices_k = Aview.indices[ak];
      double* Avals_k = Aview.values[ak];

      int C_len = 0;

      if (Aview.remote[ak]) {
	for(i=0; i<Aview.numEntriesPerRow[ak]; ++i) {
	  C_col_j[C_len] = Avals_k[i]*Bval;
          C_inds[C_len++] = A_col_inds_import[Aindices_k[i]];
	}
      }
      else {
	for(i=0; i<Aview.numEntriesPerRow[ak]; ++i) {
	  C_col_j[C_len] = Avals_k[i]*Bval;
          C_inds[C_len++] = A_col_inds[Aindices_k[i]];
	}
      }

      //Now loop across the C_col_j values and put non-zeros into C.

      for(i=0; i < C_len ; ++i) {
	if (C_col_j[i] == 0.0) continue;

	int global_row = C_inds[i];
	if (!Crowmap->MyGID(global_row)) {
	  continue;
	}

	int err = C.SumIntoGlobalValues(global_row, 1, &(C_col_j[i]), &global_col);

	if (err < 0) {
	  return(err);
	}
	else {
          if (err > 0) {
	    err = C.InsertGlobalValues(global_row, 1, &(C_col_j[i]), &global_col);
	    if (err < 0) {
              return(err);
            }
	  }
	}
      }

    }
  }

  delete [] dwork;
  delete [] iwork;

  return(0);
}

// ==============================================================
int import_and_extract_views(const Epetra_CrsMatrix& M,
			     const Epetra_Map& targetMap,
			     CrsMatrixStruct& Mview,
			     const Epetra_Import * prototypeImporter=0)
{
  //The goal of this method is to populate the 'Mview' struct with views of the
  //rows of M, including all rows that correspond to elements in 'targetMap'.
  //
  //If targetMap includes local elements that correspond to remotely-owned rows
  //of M, then those remotely-owned rows will be imported into
  //'Mview.importMatrix', and views of them will be included in 'Mview'.

  // The prototype importer, if used, has to know who owns all of the PIDs in M's rowmap.  
  // The typical use of this is to provide A's column map when I&XV is called for B, for
  // a C = A * B multiply.  

#ifdef ENABLE_MMM_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor MM(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;
  mtime=MM.getNewTimer("All I&X Alloc");
  mtime->start();
#endif
  Mview.deleteContents();

  Mview.origMatrix          = &M;
  const Epetra_Map& Mrowmap = M.RowMap();
  int numProcs              = Mrowmap.Comm().NumProc();
  Mview.numRows             = targetMap.NumMyElements();
  int* Mrows                = targetMap.MyGlobalElements();

  if (Mview.numRows > 0) {
    Mview.numEntriesPerRow = new int[Mview.numRows];
    Mview.indices = new int*[Mview.numRows];
    Mview.values  = new double*[Mview.numRows];
    Mview.remote  = new bool[Mview.numRows];
  }

  Mview.origRowMap   = &(M.RowMap());
  Mview.rowMap       = &targetMap;
  Mview.colMap       = &(M.ColMap());
  Mview.domainMap    = &(M.DomainMap());
  Mview.importColMap = NULL;
  Mview.numRemote    = 0;


#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("All I&X Extract");
  mtime->start();
#endif

  int i;
  int *rowptr=0, *colind=0;
  double *vals=0;

  EPETRA_CHK_ERR( M.ExtractCrsDataPointers(rowptr,colind,vals) );

  if(Mrowmap.SameAs(targetMap)) {
    // Short Circuit: The Row and Target Maps are the Same
    for(i=0; i<Mview.numRows; ++i) {
      Mview.numEntriesPerRow[i] = rowptr[i+1]-rowptr[i];
      Mview.indices[i]          = colind + rowptr[i];
      Mview.values[i]           = vals + rowptr[i];
      Mview.remote[i]           = false;
    }
#ifdef ENABLE_MMM_TIMINGS
    mtime->stop();
#endif
    return 0;
  }
  else if(prototypeImporter && prototypeImporter->SourceMap().SameAs(M.RowMap()) && prototypeImporter->TargetMap().SameAs(targetMap)){
    // The prototypeImporter can tell me what is local and what is remote
    int * PermuteToLIDs   = prototypeImporter->PermuteToLIDs();
    int * PermuteFromLIDs = prototypeImporter->PermuteFromLIDs();
    int * RemoteLIDs      = prototypeImporter->RemoteLIDs();
    for(int i=0; i<prototypeImporter->NumSameIDs();i++){    
      Mview.numEntriesPerRow[i] = rowptr[i+1]-rowptr[i];
      Mview.indices[i]          = colind + rowptr[i];
      Mview.values[i]           = vals + rowptr[i];
      Mview.remote[i]           = false;
    }
    for(int i=0; i<prototypeImporter->NumPermuteIDs();i++){
      int to                     = PermuteToLIDs[i];
      int from                   = PermuteFromLIDs[i];
      Mview.numEntriesPerRow[to] = rowptr[from+1]-rowptr[from];
      Mview.indices[to]          = colind + rowptr[from];
      Mview.values[to]           = vals + rowptr[from];
      Mview.remote[to]           = false;
    }
    for(int i=0; i<prototypeImporter->NumRemoteIDs();i++){
      Mview.remote[RemoteLIDs[i]] = true;
      ++Mview.numRemote;
    }
  }
  else {
    // Only LID can tell me who is local and who is remote
    for(i=0; i<Mview.numRows; ++i) {
      int mlid = Mrowmap.LID(Mrows[i]);
      if (mlid < 0) {
	Mview.remote[i] = true;
	++Mview.numRemote;
      }
      else {
	Mview.numEntriesPerRow[i] = rowptr[mlid+1]-rowptr[mlid];
	Mview.indices[i]          = colind + rowptr[mlid];
	Mview.values[i]           = vals + rowptr[mlid];
	Mview.remote[i]           = false;
      }
    }
  }

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
 mtime=MM.getNewTimer("All I&X Collective-0");
  mtime->start();
#endif


  if (numProcs < 2) {
    if (Mview.numRemote > 0) {
      cerr << "EpetraExt::MatrixMatrix::Multiply ERROR, numProcs < 2 but "
	   << "attempting to import remote matrix rows."<<endl;
      return(-1);
    }

    //If only one processor we don't need to import any remote rows, so return.
    return(0);
  }

  //
  //Now we will import the needed remote rows of M, if the global maximum
  //value of numRemote is greater than 0.
  //

  int globalMaxNumRemote = 0;
  Mrowmap.Comm().MaxAll(&Mview.numRemote, &globalMaxNumRemote, 1);

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
#endif

  if (globalMaxNumRemote > 0) {
#ifdef ENABLE_MMM_TIMINGS
    mtime=MM.getNewTimer("All I&X Import-1");
    mtime->start();
#endif
    //Create a map that describes the remote rows of M that we need.

    int* MremoteRows = Mview.numRemote>0 ? new int[Mview.numRemote] : NULL;
    int offset = 0;
    for(i=0; i<Mview.numRows; ++i) {
      if (Mview.remote[i]) {
	MremoteRows[offset++] = Mrows[i];
      }
    }

  LightweightMap MremoteRowMap(-1, Mview.numRemote, MremoteRows,Mrowmap.IndexBase());

#ifdef ENABLE_MMM_TIMINGS
    mtime->stop();
    mtime=MM.getNewTimer("All I&X Import-2");
    mtime->start();
#endif
    //Create an importer with target-map MremoteRowMap and 
    //source-map Mrowmap.
    Epetra_Import    * importer=0;
    RemoteOnlyImport * Rimporter=0;
    if(prototypeImporter && prototypeImporter->SourceMap().SameAs(M.RowMap()) && prototypeImporter->TargetMap().SameAs(targetMap)) {
      Rimporter = new RemoteOnlyImport(*prototypeImporter,MremoteRowMap);
    }
    else if(!prototypeImporter) {
      Epetra_Map MremoteRowMap2(-1, Mview.numRemote, MremoteRows,Mrowmap.IndexBase(), Mrowmap.Comm());
      importer=new Epetra_Import(MremoteRowMap2, Mrowmap);
    }
    else
      throw std::runtime_error("prototypeImporter->SourceMap() does not match M.RowMap()!");


#ifdef ENABLE_MMM_TIMINGS
    mtime->stop();
    mtime=MM.getNewTimer("All I&X Import-3");
    mtime->start();
#endif

    if(Mview.importMatrix) delete Mview.importMatrix;
    if(Rimporter) Mview.importMatrix = new LightweightCrsMatrix(M,*Rimporter);
    else Mview.importMatrix = new LightweightCrsMatrix(M,*importer);

#ifdef ENABLE_MMM_TIMINGS
    mtime->stop();
    mtime=MM.getNewTimer("All I&X Import-4");
    mtime->start();
#endif

    //Finally, use the freshly imported data to fill in the gaps in our views
    int N;
    if(Mview.importMatrix->use_lw) N = Mview.importMatrix->RowMapLW_->NumMyElements();
    else N = Mview.importMatrix->RowMapEP_->NumMyElements();

    if(N > 0) {
      rowptr = &Mview.importMatrix->rowptr_[0];
      colind = &Mview.importMatrix->colind_[0];
      vals   = &Mview.importMatrix->vals_[0];
    }


    for(i=0; i<Mview.numRows; ++i) {
      if (Mview.remote[i]) {
	int importLID = MremoteRowMap.LID(Mrows[i]);
	Mview.numEntriesPerRow[i] = rowptr[importLID+1]-rowptr[importLID];
	Mview.indices[i]          = colind + rowptr[importLID];
	Mview.values[i]           = vals + rowptr[importLID];
      }
    }


    int * MyColGIDs = (Mview.importMatrix->ColMap_.NumMyElements()>0)?Mview.importMatrix->ColMap_.MyGlobalElements():0;
    Mview.importColMap = new Epetra_Map(-1,Mview.importMatrix->ColMap_.NumMyElements(),MyColGIDs,Mview.importMatrix->ColMap_.IndexBase(),M.Comm());
    delete [] MremoteRows;
#ifdef ENABLE_MMM_TIMINGS
    mtime->stop();
#endif   

    // Cleanup
    delete Rimporter;
    delete importer;

  }
  return(0);
}




// ==============================================================
int import_only(const Epetra_CrsMatrix& M,
		const Epetra_Map& targetMap,
		CrsMatrixStruct& Mview,
		const Epetra_Import * prototypeImporter=0)
{
  // The goal of this method to populare the Mview object with ONLY the rows of M
  // that correspond to elements in 'targetMap.'  There will be no population of the
  // numEntriesPerRow, indices, values, remote or numRemote arrays.


  // The prototype importer, if used, has to know who owns all of the PIDs in M's rowmap.  
  // The typical use of this is to provide A's column map when I&XV is called for B, for
  // a C = A * B multiply.  

#ifdef ENABLE_MMM_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor MM(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;
  mtime=MM.getNewTimer("Ionly Setup");
  mtime->start();
#endif

  Mview.deleteContents();

  Mview.origMatrix          = &M;
  const Epetra_Map& Mrowmap = M.RowMap();
  int numProcs              = Mrowmap.Comm().NumProc();
  Mview.numRows             = targetMap.NumMyElements();

  Mview.origRowMap   = &(M.RowMap());
  Mview.rowMap       = &targetMap;
  Mview.colMap       = &(M.ColMap());
  Mview.domainMap    = &(M.DomainMap());
  Mview.importColMap = NULL;

  int i;
  int numRemote =0;
  int N = Mview.rowMap->NumMyElements();

  if(Mrowmap.SameAs(targetMap)) {
    numRemote = 0;
    Mview.targetMapToOrigRow.resize(N); 
    for(i=0;i<N; i++) Mview.targetMapToOrigRow[i]=i;

#ifdef ENABLE_MMM_TIMINGS      
  mtime->stop();   
#endif
    return 0;
  }
  else if(prototypeImporter && prototypeImporter->SourceMap().SameAs(M.RowMap()) && prototypeImporter->TargetMap().SameAs(targetMap)){
    numRemote = prototypeImporter->NumRemoteIDs();

    Mview.targetMapToOrigRow.resize(N);    Mview.targetMapToOrigRow.assign(N,-1);
    Mview.targetMapToImportRow.resize(N);  Mview.targetMapToImportRow.assign(N,-1);

    const int * PermuteToLIDs   = prototypeImporter->PermuteToLIDs();
    const int * PermuteFromLIDs = prototypeImporter->PermuteFromLIDs();
    const int * RemoteLIDs      = prototypeImporter->RemoteLIDs();

    for(i=0; i<prototypeImporter->NumSameIDs();i++)
      Mview.targetMapToOrigRow[i] = i;

    // NOTE: These are reversed on purpose since we're doing a revere map.
    for(i=0; i<prototypeImporter->NumPermuteIDs();i++)
      Mview.targetMapToOrigRow[PermuteToLIDs[i]] = PermuteFromLIDs[i];
    
    for(i=0; i<prototypeImporter->NumRemoteIDs();i++)
      Mview.targetMapToImportRow[RemoteLIDs[i]] = i;

  }
  else
    throw "import_only: This routine only works if you either have the right map or no prototypeImporter";

  if (numProcs < 2) {
    if (Mview.numRemote > 0) {
      cerr << "EpetraExt::MatrixMatrix::Multiply ERROR, numProcs < 2 but "
	   << "attempting to import remote matrix rows."<<endl;
      return(-1);
    }
#ifdef ENABLE_MMM_TIMINGS      
    mtime->stop();  
#endif

    //If only one processor we don't need to import any remote rows, so return.
    return(0);
  }

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("Ionly Import-1");
  mtime->start();
#endif
  const int * RemoteLIDs = prototypeImporter->RemoteLIDs();
  
    //Create a map that describes the remote rows of M that we need.
  int* MremoteRows = numRemote>0 ? new int[prototypeImporter->NumRemoteIDs()] : 0;
  for(i=0; i<prototypeImporter->NumRemoteIDs(); i++)
    MremoteRows[i] = targetMap.GID(RemoteLIDs[i]);
  
  LightweightMap MremoteRowMap(-1, numRemote, MremoteRows,Mrowmap.IndexBase());

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("Ionly Import-2");
  mtime->start();
#endif
  //Create an importer with target-map MremoteRowMap and 
  //source-map Mrowmap.
  RemoteOnlyImport * Rimporter=0;
  Rimporter = new RemoteOnlyImport(*prototypeImporter,MremoteRowMap);
  
#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("Ionly Import-3");
  mtime->start();
#endif
  
  Mview.importMatrix = new LightweightCrsMatrix(M,*Rimporter);
  
#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=MM.getNewTimer("Ionly Import-4");
  mtime->start();
#endif

  // Cleanup
  delete Rimporter;
  delete [] MremoteRows;
#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
#endif      
  

  return(0);
}




//=========================================================================
int form_map_union(const Epetra_Map* map1,
		   const Epetra_Map* map2,
		   const Epetra_Map*& mapunion)
{
  //form the union of two maps

  if (map1 == NULL) {
    mapunion = new Epetra_Map(*map2);
    return(0);
  }

  if (map2 == NULL) {
    mapunion = new Epetra_Map(*map1);
    return(0);
  }

  int map1_len       = map1->NumMyElements();
  int* map1_elements = map1->MyGlobalElements();
  int map2_len       = map2->NumMyElements();
  int* map2_elements = map2->MyGlobalElements();

  int* union_elements = new int[map1_len+map2_len];

  int map1_offset = 0, map2_offset = 0, union_offset = 0;

  while(map1_offset < map1_len && map2_offset < map2_len) {
    int map1_elem = map1_elements[map1_offset];
    int map2_elem = map2_elements[map2_offset];

    if (map1_elem < map2_elem) {
      union_elements[union_offset++] = map1_elem;
      ++map1_offset;
    }
    else if (map1_elem > map2_elem) {
      union_elements[union_offset++] = map2_elem;
      ++map2_offset;
    }
    else {
      union_elements[union_offset++] = map1_elem;
      ++map1_offset;
      ++map2_offset;
    }
  }

  int i;
  for(i=map1_offset; i<map1_len; ++i) {
    union_elements[union_offset++] = map1_elements[i];
  }

  for(i=map2_offset; i<map2_len; ++i) {
    union_elements[union_offset++] = map2_elements[i];
  }

  mapunion = new Epetra_Map(-1, union_offset, union_elements,
			    map1->IndexBase(), map1->Comm());

  delete [] union_elements;

  return(0);
}

//=========================================================================
Epetra_Map* find_rows_containing_cols(const Epetra_CrsMatrix& M,
                                      const Epetra_Map& column_map)
{
  //The goal of this function is to find all rows in the matrix M that contain
  //column-indices which are in 'column_map'. A map containing those rows is
  //returned. (The returned map will then be used as the source row-map for
  //importing those rows into an overlapping distribution.)

  std::pair<int,int> my_col_range = get_col_range(column_map);

  const Epetra_Comm& Comm = M.RowMap().Comm();
  int num_procs = Comm.NumProc();
  int my_proc = Comm.MyPID();
  std::vector<int> send_procs;
  send_procs.reserve(num_procs-1);
  std::vector<int> col_ranges;
  col_ranges.reserve(2*(num_procs-1));
  for(int p=0; p<num_procs; ++p) {
    if (p == my_proc) continue;
    send_procs.push_back(p);
    col_ranges.push_back(my_col_range.first);
    col_ranges.push_back(my_col_range.second);
  }

  Epetra_Distributor* distor = Comm.CreateDistributor();

  int num_recv_procs = 0;
  int num_send_procs = send_procs.size();
  bool deterministic = true;
  int* send_procs_ptr = send_procs.size()>0 ? &send_procs[0] : NULL;
  distor->CreateFromSends(num_send_procs, send_procs_ptr, deterministic, num_recv_procs);

  int len_import_chars = 0;
  char* import_chars = NULL;

  char* export_chars = col_ranges.size()>0 ? reinterpret_cast<char*>(&col_ranges[0]) : NULL;
  int err = distor->Do(export_chars, 2*sizeof(int), len_import_chars, import_chars);
  if (err != 0) {
    std::cout << "ERROR returned from Distributor::Do."<<std::endl;
  }
 
  int* IntImports = reinterpret_cast<int*>(import_chars);
  int num_import_pairs = len_import_chars/(2*sizeof(int));
  int offset = 0;
  std::vector<int> send_procs2;
  send_procs2.reserve(num_procs);
  std::vector<int> proc_col_ranges;
  std::pair<int,int> M_col_range = get_col_range(M);
  for(int i=0; i<num_import_pairs; ++i) {
    int first_col = IntImports[offset++];
    int last_col =  IntImports[offset++];
    if (first_col <= M_col_range.second && last_col >= M_col_range.first) {
      send_procs2.push_back(send_procs[i]);
      proc_col_ranges.push_back(first_col);
      proc_col_ranges.push_back(last_col);
    }
  }

  std::vector<int> send_rows;
  std::vector<int> rows_per_send_proc;
  pack_outgoing_rows(M, proc_col_ranges, send_rows, rows_per_send_proc);

  Epetra_Distributor* distor2 = Comm.CreateDistributor();

  int* send_procs2_ptr = send_procs2.size()>0 ? &send_procs2[0] : NULL;
  num_recv_procs = 0;
  err = distor2->CreateFromSends(send_procs2.size(), send_procs2_ptr, deterministic, num_recv_procs);

  export_chars = send_rows.size()>0 ? reinterpret_cast<char*>(&send_rows[0]) : NULL;
  int* rows_per_send_proc_ptr = rows_per_send_proc.size()>0 ? &rows_per_send_proc[0] : NULL;
  len_import_chars = 0;
  err = distor2->Do(export_chars, (int)sizeof(int), rows_per_send_proc_ptr, len_import_chars, import_chars);

  int* recvd_row_ints = reinterpret_cast<int*>(import_chars);
  int num_recvd_rows = len_import_chars/sizeof(int);

  Epetra_Map recvd_rows(-1, num_recvd_rows, recvd_row_ints, column_map.IndexBase(), Comm);

  delete distor;
  delete distor2;
  delete [] import_chars;

  Epetra_Map* result_map = NULL;

  err = form_map_union(&(M.RowMap()), &recvd_rows, (const Epetra_Map*&)result_map);
  if (err != 0) {
    return(NULL);
  }

  return(result_map);
}

//=========================================================================
int MatrixMatrix::Multiply(const Epetra_CrsMatrix& A,
			   bool transposeA,
			   const Epetra_CrsMatrix& B,
			   bool transposeB,
			   Epetra_CrsMatrix& C,
                           bool call_FillComplete_on_result)
{

  // DEBUG
  //  bool NewFlag=!C.IndicesAreLocal() && !C.IndicesAreGlobal();
#ifdef ENABLE_MMM_TIMINGS
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor M(myTime);
  Teuchos::RCP<Teuchos::Time> mtime;  
  mtime=M.getNewTimer("All Setup");
  mtime->start();
#endif
  // end DEBUG

  //
  //This method forms the matrix-matrix product C = op(A) * op(B), where
  //op(A) == A   if transposeA is false,
  //op(A) == A^T if transposeA is true,
  //and similarly for op(B).
  //

  //A and B should already be Filled.
  //(Should we go ahead and call FillComplete() on them if necessary?
  // or error out? For now, we choose to error out.)
  if (!A.Filled() || !B.Filled()) {
    EPETRA_CHK_ERR(-1);
  }

  //We're going to refer to the different combinations of op(A) and op(B)
  //as scenario 1 through 4.

  int scenario = 1;//A*B
  if (transposeB && !transposeA) scenario = 2;//A*B^T
  if (transposeA && !transposeB) scenario = 3;//A^T*B
  if (transposeA && transposeB)  scenario = 4;//A^T*B^T

  //now check size compatibility
  int Aouter = transposeA ? A.NumGlobalCols() : A.NumGlobalRows();
  int Bouter = transposeB ? B.NumGlobalRows() : B.NumGlobalCols();
  int Ainner = transposeA ? A.NumGlobalRows() : A.NumGlobalCols();
  int Binner = transposeB ? B.NumGlobalCols() : B.NumGlobalRows();
  if (Ainner != Binner) {
    cerr << "MatrixMatrix::Multiply: ERROR, inner dimensions of op(A) and op(B) "
         << "must match for matrix-matrix product. op(A) is "
         <<Aouter<<"x"<<Ainner << ", op(B) is "<<Binner<<"x"<<Bouter<<endl;
    return(-1);
  }

  //The result matrix C must at least have a row-map that reflects the
  //correct row-size. Don't check the number of columns because rectangular
  //matrices which were constructed with only one map can still end up
  //having the correct capacity and dimensions when filled.
  if (Aouter > C.NumGlobalRows()) {
    cerr << "MatrixMatrix::Multiply: ERROR, dimensions of result C must "
         << "match dimensions of op(A) * op(B). C has "<<C.NumGlobalRows()
         << " rows, should have at least "<<Aouter << endl;
    return(-1);
  }

  //It doesn't matter whether C is already Filled or not. If it is already
  //Filled, it must have space allocated for the positions that will be
  //referenced in forming C = op(A)*op(B). If it doesn't have enough space,
  //we'll error out later when trying to store result values.

  //We're going to need to import remotely-owned sections of A and/or B
  //if more than 1 processor is performing this run, depending on the scenario.
  int numProcs = A.Comm().NumProc();

  //If we are to use the transpose of A and/or B, we'll need to be able to 
  //access, on the local processor, all rows that contain column-indices in
  //the domain-map.
  const Epetra_Map* domainMap_A = &(A.DomainMap());
  const Epetra_Map* domainMap_B = &(B.DomainMap());

  const Epetra_Map* rowmap_A = &(A.RowMap());
  const Epetra_Map* rowmap_B = &(B.RowMap());

  //Declare some 'work-space' maps which may be created depending on
  //the scenario, and which will be deleted before exiting this function.
  const Epetra_Map* workmap1 = NULL;
  const Epetra_Map* workmap2 = NULL;
  const Epetra_Map* mapunion1 = NULL;

  //Declare a couple of structs that will be used to hold views of the data
  //of A and B, to be used for fast access during the matrix-multiplication.
  CrsMatrixStruct Aview;
  CrsMatrixStruct Bview;

  const Epetra_Map* targetMap_A = rowmap_A;
  const Epetra_Map* targetMap_B = rowmap_B;

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=M.getNewTimer("All I&X");
  mtime->start();
#endif
  if (numProcs > 1) {
    //If op(A) = A^T, find all rows of A that contain column-indices in the
    //local portion of the domain-map. (We'll import any remote rows
    //that fit this criteria onto the local processor.)
    if (transposeA) {
      workmap1 = find_rows_containing_cols(A, *domainMap_A);
      targetMap_A = workmap1;
    }
  }

  //Now import any needed remote rows and populate the Aview struct.
  if(scenario==1 && call_FillComplete_on_result) {
    EPETRA_CHK_ERR(import_only(A,*targetMap_A,Aview));
  }
  else  {
    EPETRA_CHK_ERR( import_and_extract_views(A, *targetMap_A, Aview));
  }


  // NOTE:  Next up is to switch to import_only for B as well, and then modify the THREE SerialCores
  // to add a Acol2Brow and Acol2Bimportrow array for in-algorithm lookups.  
  

  // Make sure B's views are consistent with A even in serial.
  const Epetra_Map* colmap_op_A = NULL;
  if(scenario==1 || numProcs > 1){
    if (transposeA) {
      colmap_op_A = targetMap_A;
    }
    else {
      colmap_op_A = &(A.ColMap());
    }
    targetMap_B = colmap_op_A;
  }

  if (numProcs > 1) {
    //If op(B) = B^T, find all rows of B that contain column-indices in the
    //local-portion of the domain-map, or in the column-map of op(A).
    //We'll import any remote rows that fit this criteria onto the
    //local processor.
    if (transposeB) {
      EPETRA_CHK_ERR( form_map_union(colmap_op_A, domainMap_B, mapunion1) );
      workmap2 = find_rows_containing_cols(B, *mapunion1);
      targetMap_B = workmap2;
    }
  }

  //Now import any needed remote rows and populate the Bview struct.  
  if(scenario==1 && call_FillComplete_on_result) {
    EPETRA_CHK_ERR(import_only(B,*targetMap_B,Bview,A.Importer()));
  }
  else {
    EPETRA_CHK_ERR( import_and_extract_views(B, *targetMap_B, Bview) );
  }

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
  mtime=M.getNewTimer("All Multiply");
  mtime->start();
#endif

  // Zero if filled
  if(C.Filled()) C.PutScalar(0.0);

  //Now call the appropriate method to perform the actual multiplication.
  CrsWrapper_Epetra_CrsMatrix ecrsmat(C);

  switch(scenario) {
  case 1:    EPETRA_CHK_ERR( mult_A_B(A,Aview,B,Bview,C,call_FillComplete_on_result) );
    break;
  case 2:    EPETRA_CHK_ERR( mult_A_Btrans(Aview, Bview, ecrsmat) );
    break;
  case 3:    EPETRA_CHK_ERR( mult_Atrans_B(Aview, Bview, ecrsmat) );
    break;
  case 4:    EPETRA_CHK_ERR( mult_Atrans_Btrans(Aview, Bview, ecrsmat) );
    break;
  }


  if (scenario != 1 && call_FillComplete_on_result) {
    //We'll call FillComplete on the C matrix before we exit, and give
    //it a domain-map and a range-map.
    //The domain-map will be the domain-map of B, unless
    //op(B)==transpose(B), in which case the range-map of B will be used.
    //The range-map will be the range-map of A, unless
    //op(A)==transpose(A), in which case the domain-map of A will be used.

    const Epetra_Map* domainmap =
      transposeB ? &(B.RangeMap()) : &(B.DomainMap());

    const Epetra_Map* rangemap =
      transposeA ? &(A.DomainMap()) : &(A.RangeMap());

    if (!C.Filled()) {
      EPETRA_CHK_ERR( C.FillComplete(*domainmap, *rangemap) );
    }
  }

  //Finally, delete the objects that were potentially created
  //during the course of importing remote sections of A and B.

  delete mapunion1; mapunion1 = NULL;
  delete workmap1; workmap1 = NULL;
  delete workmap2; workmap2 = NULL;

#ifdef ENABLE_MMM_TIMINGS
  mtime->stop();
#endif


  return(0);
}

//=========================================================================
int MatrixMatrix::Add(const Epetra_CrsMatrix& A,
                      bool transposeA,
                      double scalarA,
                      Epetra_CrsMatrix& B,
                      double scalarB )
{
  //
  //This method forms the matrix-matrix sum B = scalarA * op(A) + scalarB * B, where

  //A should already be Filled. It doesn't matter whether B is
  //already Filled, but if it is, then its graph must already contain
  //all nonzero locations that will be referenced in forming the
  //sum.

  if (!A.Filled() ) {
    std::cerr << "EpetraExt::MatrixMatrix::Add ERROR, input matrix A.Filled() is false, it is required to be true. (Result matrix B is not required to be Filled())."<<std::endl;
    EPETRA_CHK_ERR(-1);
  }

  //explicit tranpose A formed as necessary
  Epetra_CrsMatrix * Aprime = 0;
  EpetraExt::RowMatrix_Transpose * Atrans = 0;
  if( transposeA )
  {
    Atrans = new EpetraExt::RowMatrix_Transpose();
    Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
    Aprime = const_cast<Epetra_CrsMatrix*>(&A);

  int MaxNumEntries = EPETRA_MAX( A.MaxNumEntries(), B.MaxNumEntries() );
  int A_NumEntries, B_NumEntries;
  int * A_Indices = new int[MaxNumEntries];
  double * A_Values = new double[MaxNumEntries];
  int* B_Indices;
  double* B_Values;

  int NumMyRows = B.NumMyRows();
  int Row, err;
  int ierr = 0;

  if( scalarA )
  {
    //Loop over B's rows and sum into
    for( int i = 0; i < NumMyRows; ++i )
    {
      Row = B.GRID(i);
      EPETRA_CHK_ERR( Aprime->ExtractGlobalRowCopy( Row, MaxNumEntries, A_NumEntries, A_Values, A_Indices ) );

      if (scalarB != 1.0) {
        if (!B.Filled()) {
          EPETRA_CHK_ERR( B.ExtractGlobalRowView( Row, B_NumEntries,
                                                  B_Values, B_Indices));
        }
        else {
          EPETRA_CHK_ERR( B.ExtractMyRowView( i, B_NumEntries,
                                              B_Values, B_Indices));
        }

        for(int jj=0; jj<B_NumEntries; ++jj) {
          B_Values[jj] = scalarB*B_Values[jj];
        }
      }

      if( scalarA != 1.0 ) {
        for( int j = 0; j < A_NumEntries; ++j ) A_Values[j] *= scalarA;
      }

      if( B.Filled() ) {//Sum In Values
        err = B.SumIntoGlobalValues( Row, A_NumEntries, A_Values, A_Indices );
        assert( err >= 0 );
        if (err < 0) ierr = err;
      }
      else {
        err = B.InsertGlobalValues( Row, A_NumEntries, A_Values, A_Indices );
        assert( err == 0 || err == 1 || err == 3 );
        if (err < 0) ierr = err;
      }
    }
  }
  else {
    EPETRA_CHK_ERR( B.Scale(scalarB) );
  }

  delete [] A_Indices;
  delete [] A_Values;

  if( Atrans ) delete Atrans;

  return(ierr);
}

int MatrixMatrix::Add(const Epetra_CrsMatrix& A,
                      bool transposeA,
                      double scalarA,
                      const Epetra_CrsMatrix & B,
                      bool transposeB,
                      double scalarB,
                      Epetra_CrsMatrix * & C)
{
  //
  //This method forms the matrix-matrix sum C = scalarA * op(A) + scalarB * op(B), where

  //A and B should already be Filled. C should be an empty pointer.

  if (!A.Filled() || !B.Filled() ) {
     std::cerr << "EpetraExt::MatrixMatrix::Add ERROR, input matrix A.Filled() or B.Filled() is false,"
               << "they are required to be true. (Result matrix C should be an empty pointer)" << std::endl;
     EPETRA_CHK_ERR(-1);
  }

  Epetra_CrsMatrix * Aprime = 0, * Bprime=0;
  EpetraExt::RowMatrix_Transpose * Atrans = 0,* Btrans = 0;

  //explicit tranpose A formed as necessary
  if( transposeA ) {
     Atrans = new EpetraExt::RowMatrix_Transpose();
     Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
     Aprime = const_cast<Epetra_CrsMatrix*>(&A);

  //explicit tranpose B formed as necessary
  if( transposeB ) {
     Btrans = new EpetraExt::RowMatrix_Transpose();
     Bprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Btrans)(const_cast<Epetra_CrsMatrix&>(B)))));
  }
  else
     Bprime = const_cast<Epetra_CrsMatrix*>(&B);

  // allocate or zero the new matrix
  if(C!=0)
     C->PutScalar(0.0);
  else
     C = new Epetra_CrsMatrix(Copy,Aprime->RowMap(),0);

  // build arrays  for easy resuse
  int ierr = 0;
  Epetra_CrsMatrix * Mat[] = { Aprime,Bprime};
  double scalar[] = { scalarA, scalarB};

  // do a loop over each matrix to add: A reordering might be more efficient
  for(int k=0;k<2;k++) {
     int MaxNumEntries = Mat[k]->MaxNumEntries();
     int NumEntries;
     int * Indices = new int[MaxNumEntries];
     double * Values = new double[MaxNumEntries];
   
     int NumMyRows = Mat[k]->NumMyRows();
     int Row, err;
     int ierr = 0;
   
     //Loop over rows and sum into C
     for( int i = 0; i < NumMyRows; ++i ) {
        Row = Mat[k]->GRID(i);
        EPETRA_CHK_ERR( Mat[k]->ExtractGlobalRowCopy( Row, MaxNumEntries, NumEntries, Values, Indices));
   
        if( scalar[k] != 1.0 )
           for( int j = 0; j < NumEntries; ++j ) Values[j] *= scalar[k];
   
        if(C->Filled()) { // Sum in values
           err = C->SumIntoGlobalValues( Row, NumEntries, Values, Indices );
           if (err < 0) ierr = err;
        } else { // just add it to the unfilled CRS Matrix
           err = C->InsertGlobalValues( Row, NumEntries, Values, Indices );
           if (err < 0) ierr = err;
        }
     }

     delete [] Indices;
     delete [] Values;
  }

  if( Atrans ) delete Atrans;
  if( Btrans ) delete Btrans;

  return(ierr);
}


} // namespace EpetraExt

