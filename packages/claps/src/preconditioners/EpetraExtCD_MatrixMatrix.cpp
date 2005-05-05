//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
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
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#include <Epetra_ConfigDefs.h>
#include "EpetraExtCD_MatrixMatrix.hpp"
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Util.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <algorithm>
using namespace std;

namespace EpetraExtCD {

int MatrixMatrix::Multiply(const Epetra_CrsMatrix& A,
			   const Epetra_CrsMatrix& B,
			   Epetra_CrsMatrix* & C)
{
  //
  // This function forms the matrix-matrix product C = A * B. This is
  // a modification of EpetraExt MatrixMatrix::Multiply that Alan Williams
  // told me about.
  //

  // A and B should already be Filled. It is assumed that C has not been
  // constructed yet.

  if (!A.Filled() || !B.Filled()) {
    EPETRA_CHK_ERR(-1);
  }

  int i, j, k, err;

  int AnumRows = A.NumMyRows();

  //
  // We  need views of all rows of B corresponding to the column-map of A.
  //

  const Epetra_Map& Acolmap = A.ColMap();
  const Epetra_Map& Browmap = B.RowMap();

  int BnumRows = Acolmap.NumMyElements();
  int* Brows = Acolmap.MyGlobalElements();

  int *BnumEntriesPerRow = NULL;
  int** Bindices = NULL;
  double** Bvalues = NULL;

  //
  // As we extract views of the rows of B, keep track of which elements of
  // Acolmap aren't local rows of B because we'll have to import those next.
  //
  bool* remote = NULL;
  
  if (BnumRows > 0) {
    BnumEntriesPerRow = new int[BnumRows];
    Bindices = new int*[BnumRows];
    Bvalues = new double*[BnumRows];
    remote = new bool[BnumRows];
  }

  const Epetra_Map* Bcolmap = &(B.ColMap());

  int numRemote = 0;

  for(i=0; i<BnumRows; ++i) {
    int blid = Browmap.LID(Brows[i]);
    if (blid < 0) {
      remote[i] = true;
      ++numRemote;
    }
    else {
      EPETRA_CHK_ERR(B.ExtractMyRowView(blid, BnumEntriesPerRow[i],
					Bvalues[i], Bindices[i]) );
      remote[i] = false;
    }
  }

  //
  // Now we will import the needed remote rows of B, if the global maximum
  // value of numRemote is greater than 0.
  //

  int globalMaxNumRemote = 0;
  Acolmap.Comm().MaxAll(&numRemote, &globalMaxNumRemote, 1);

  Epetra_CrsMatrix* importB = NULL;
  const Epetra_Map* importBcolmap = NULL;

  if (globalMaxNumRemote > 0) {
    // Create a map that describes the remote rows of B that we need.

    int* BremoteRows = numRemote>0 ? new int[numRemote] : NULL;
    int offset = 0;
    for(i=0; i<BnumRows; ++i) {
      if (remote[i]) {
	BremoteRows[offset++] = Brows[i];
      }
    }

    Epetra_Map BremoteRowMap(-1, numRemote, BremoteRows,
			     Acolmap.IndexBase(), Acolmap.Comm());

    // Create an importer with target-map BremoteRowMap and
    // source-map Browmap
    Epetra_Import importer(BremoteRowMap, Browmap);

    // Now create a new matrix into which we can import the remote rows of B
    // that we need.
    importB = new Epetra_CrsMatrix(Copy, BremoteRowMap, 1);

    EPETRA_CHK_ERR( importB->Import(B, importer, Insert) );

    EPETRA_CHK_ERR( importB->FillComplete(B.DomainMap(), B.RangeMap()) );

    importBcolmap = &(importB->ColMap());

    // Finally, use the freshly imported data to fill in the gaps in our views
    // of rows of B.

    for(i=0; i<BnumRows; ++i) {
      if (remote[i]) {
	int importLID = BremoteRowMap.LID(Brows[i]);
	EPETRA_CHK_ERR( importB->ExtractMyRowView(importLID, 
			      BnumEntriesPerRow[i], Bvalues[i], Bindices[i]) );
      }
    }
    delete [] BremoteRows;
  }

  //
  // construct maps from the columns in Bcolmap and importBcolmap to local
  // columns of C
  //

  int ncolB = Bcolmap->NumMyPoints();
  int B_firstCol = Bcolmap->MinLID();
  int *Bglobal = NULL; int *Bmap = NULL;
  if (ncolB > 0) {
    Bglobal = new int[ncolB];
    Bcolmap->MyGlobalElements(Bglobal);
    Bmap = new int[ncolB];
  }
  int ncolB_import = 0;
  int B_firstCol_import = 0;
  int *Bglobal_import = NULL; int *Bmap_import = NULL;
  if (importBcolmap != NULL) {
    ncolB_import = importBcolmap->NumMyPoints();
    B_firstCol_import = importBcolmap->MinLID();
    Bglobal_import = new int[ncolB_import];
    importBcolmap->MyGlobalElements(Bglobal_import);
    Bmap_import = new int[ncolB_import];
  }
  int n_all = ncolB + ncolB_import;
  int *gcol_all = NULL;
  int ncolC = 0;
  if (n_all > 0) {
    gcol_all = new int[n_all];
    for (i=0; i<ncolB; i++) gcol_all[i] = Bglobal[i];
    for (i=0; i<ncolB_import; i++) gcol_all[ncolB+i] = Bglobal_import[i];
    sort(&gcol_all[0], &gcol_all[n_all]);
    int previd = gcol_all[0];
    ncolC = 1;
    for (i=1; i<n_all; i++) {
      if (gcol_all[i] != previd) {
	previd = gcol_all[i];
	gcol_all[ncolC] = gcol_all[i];
	ncolC++;
      }
    }
    for (i=0; i<ncolB; i++) {
      int *loc = lower_bound(&gcol_all[0], &gcol_all[ncolC], Bglobal[i]);
      Bmap[i] = (int) (loc - &gcol_all[0]);
    }
    for (i=0; i<ncolB_import; i++) {
      int *loc = lower_bound(&gcol_all[0], &gcol_all[ncolC], 
			     Bglobal_import[i]);
      Bmap_import[i] = (int) (loc - &gcol_all[0]);
    }
  }
  //  int MyPID =   Acolmap.Comm().MyPID();
  delete [] Bglobal;
  delete [] Bglobal_import;
  //
  // construct C matrix
  //
  Epetra_Map ColMap(-1, ncolC, gcol_all, 0, Acolmap.Comm());
  C = new Epetra_CrsMatrix(Copy, A.RowMap(), ColMap, 0);
  delete [] gcol_all;
  //
  // fill in values in C matrix
  //
  double *rowvalues = NULL; if (ncolC > 0) rowvalues = new double[ncolC];
  int *rowflag = NULL; if (ncolC > 0) rowflag = new int[ncolC];
  int *activecols = NULL; if (ncolC > 0) activecols = new int[ncolC];
  for (i=0; i<ncolC; i++) rowflag[i] = -1;
  for (i=0; i<AnumRows; ++i) {
    int ncol_active = 0;
    int AnumEntriesPerRow, *Aindices;
    double *Avalues;
    EPETRA_CHK_ERR(A.ExtractMyRowView(i, AnumEntriesPerRow, Avalues, Aindices));
		  
    // loop across the i-th row of A and for each corresponding row
    // in B, loop across columns and accumulate the product A(i,k)*B(k,j).
    
    for (k=0; k<AnumEntriesPerRow; ++k) {
      int Ak = Aindices[k];
      double Aval = Avalues[k];

      int* Bcol_inds = Bindices[Ak];
      double* Bvals_k = Bvalues[Ak];

      if (remote[Ak]) {
	for (j=0; j<BnumEntriesPerRow[Ak]; ++j) {
	  int loc = Bcol_inds[j] - B_firstCol_import;
	  loc = Bmap_import[loc];
	  if (rowflag[loc] == -1) {
	    rowflag[loc] = ncol_active;
	    activecols[ncol_active] = loc;
	    rowvalues[ncol_active] = Aval*Bvals_k[j];
	    ncol_active++;
	  }
	  else rowvalues[rowflag[loc]] += Aval*Bvals_k[j];
	}
      }
      else {
	for(j=0; j<BnumEntriesPerRow[Ak]; ++j) {
	  int loc = Bcol_inds[j] - B_firstCol;
	  loc = Bmap[loc];
	  if (rowflag[loc] == -1) {
	    rowflag[loc] = ncol_active;
	    activecols[ncol_active] = loc;
	    rowvalues[ncol_active] = Aval*Bvals_k[j];
	    ncol_active++;
	  }
	  else rowvalues[rowflag[loc]] += Aval*Bvals_k[j];
	}
      }
    }
    err = C->InsertMyValues(i, ncol_active, rowvalues, activecols);
    if (err < 0) return(err);
    for (j=0; j<ncol_active; j++) rowflag[activecols[j]] = -1;
  } 
  C->FillComplete(B.DomainMap(), A.RangeMap());

  delete [] BnumEntriesPerRow;
  delete [] Bindices;
  delete [] Bvalues;
  delete [] remote;
  delete importB;
  delete [] Bmap;
  delete [] Bmap_import;
  delete [] rowvalues;
  delete [] rowflag;
  delete [] activecols;

  return(0);
}
} // namespace EpetraExtCD

