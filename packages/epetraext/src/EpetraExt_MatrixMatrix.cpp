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

#include <EpetraExt_MatrixMatrix.h>

#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Util.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>

namespace EpetraExt {

int MatrixMatrix::Multiply(const Epetra_CrsMatrix& A,
			   const Epetra_CrsMatrix& B,
			   Epetra_CrsMatrix& C)
{
  //
  //This function forms the matrix-matrix product C = A * B.
  //

  //A and B should already be Filled. As we'll see below, it doesn't matter
  //whether C is already Filled, but if it is, then its graph must already
  //contain all nonzero locations that will be referenced in forming the
  //product.

  if (!A.Filled() || !B.Filled()) {
    EPETRA_CHK_ERR(-1);
  }

  int i, j, k;

  //First, extract views of the contents of A for fast access later.
  int AnumRows = A.NumMyRows();

  int *AnumEntriesPerRow = NULL;
  int** Aindices = NULL;
  double** Avalues = NULL;

  if (AnumRows > 0) {
    AnumEntriesPerRow = new int[AnumRows];
    Aindices = new int*[AnumRows];
    Avalues = new double*[AnumRows];
  }

  for(i=0; i<AnumRows; ++i) {
    EPETRA_CHK_ERR( A.ExtractMyRowView(i, AnumEntriesPerRow[i],
				       Avalues[i], Aindices[i]) );
  }

  //
  //We will also need views of all rows of B corresponding to the column-map of A.
  //

  const Epetra_Map& Acolmap = A.ColMap();
  const Epetra_Map& Browmap = B.RowMap();

  int BnumRows = Acolmap.NumMyElements();
  int* Brows = Acolmap.MyGlobalElements();

  int *BnumEntriesPerRow = BnumRows>0 ? new int[BnumRows] : NULL;
  int** Bindices = NULL;
  double** Bvalues = NULL;

  //
  //As we extract views of the rows of B, keep track of which elements of
  //Acolmap aren't local rows of B because we'll have to import those next.
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
  //Now we will import the needed remote rows of B, if the global maximum
  //value of numRemote is greater than 0.
  //

  int globalMaxNumRemote = 0;
  Acolmap.Comm().MaxAll(&numRemote, &globalMaxNumRemote, 1);

  Epetra_CrsMatrix* importB = NULL;
  const Epetra_Map* importBcolmap = NULL;

  if (globalMaxNumRemote > 0) {
    //Create a map that describes the remote rows of B that we need.

    int* BremoteRows = numRemote>0 ? new int[numRemote] : NULL;
    int offset = 0;
    for(i=0; i<BnumRows; ++i) {
      if (remote[i]) {
	BremoteRows[offset++] = Brows[i];
      }
    }

    Epetra_Map BremoteRowMap(-1, numRemote, BremoteRows,
			   Acolmap.IndexBase(), Acolmap.Comm());

    //Create an importer with target-map BremoteRowMap and
    //source-map Browmap
    Epetra_Import importer(BremoteRowMap, Browmap);

    //Now create a new matrix into which we can import the remote rows of B
    //that we need.
    importB = new Epetra_CrsMatrix(Copy, BremoteRowMap, 1);

    EPETRA_CHK_ERR( importB->Import(B, importer, Insert) );

    EPETRA_CHK_ERR( importB->FillComplete(B.DomainMap(), B.RangeMap()) );

    importBcolmap = &(importB->ColMap());

    //Finally, use the freshly imported data to fill in the gaps in our views
    //of rows of B.
    for(i=0; i<BnumRows; ++i) {
      if (remote[i]) {
	int importLID = BremoteRowMap.LID(Brows[i]);
	EPETRA_CHK_ERR( importB->ExtractMyRowView(importLID, BnumEntriesPerRow[i],
						  Bvalues[i], Bindices[i]) );
      }
    }

    delete [] BremoteRows;
  }

  //zero the result matrix before we start the calculations.
  EPETRA_CHK_ERR( C.PutScalar(0.0) );

  //
  //A simple algorithm for forming a matrix-matrix product C = A*B is:
  //
  //for(i in rows-of-A) {
  //  for(j in columns-of-B) {
  //    for(k in column-indices of A's i-th row) {
  //      C[i,j] += A[i,k]*B[k,j]
  //    }
  //  }
  //}
  //
  //The following code is based on this. However, for sparse row-oriented
  //matrices, the expression in the innermost loop requires searching for 'j'
  //in the non-dense list of column-indices for B's k-th row. Profiling this
  //code on moderately sized matrices quickly revealed that approach to be
  //prohibitively slow.
  //So instead, the "for(j..." loop has been removed and the innermost expression
  //is replaced with code that loops over B's k-th row and updates a list of
  //values for C's i-th row and no searching is required.
  //

  //First, set up the list of values that will represent temporary storage for the
  //i-th row of C. We'll need two separate lists, because we have two different
  //column-spaces for B to work with. We are working with local indices, and the
  //local numbering of B is different than the local numbering of importB.

  int C_firstCol = Bcolmap->MinLID();
  int C_lastCol = Bcolmap->MaxLID();

  int C_firstCol_import = 0;
  int C_lastCol_import = -1;

  if (importBcolmap != NULL) {
    C_firstCol_import = importBcolmap->MinLID();
    C_lastCol_import = importBcolmap->MaxLID();
  }

  int C_numCols = C_lastCol - C_firstCol + 1;
  double* C_row_i = new double[C_numCols];

  int C_numCols_import = C_lastCol_import - C_firstCol_import + 1;

  double* C_row_i_import = C_numCols_import > 0 ?
                           new double[C_numCols_import] : NULL;

  for(j=0; j<C_numCols; ++j) {
    C_row_i[j] = 0.0;
  }

  for(j=0; j<C_numCols_import; ++j) {
    C_row_i_import[j] = 0.0;
  }

  //
  //Now we're finally ready to begin the calculations.
  //

  //loop over the rows of A.
  for(i=0; i<AnumRows; ++i) {
    int* Aindices_i = Aindices[i];
    double* Aval_i  = Avalues[i];

    //loop across the i-th row of A and for each corresponding row
    //in B, loop across colums and accumulate product
    //A(i,k)*B(k,j) into our partial sum quantities C_row_i.

    for(k=0; k<AnumEntriesPerRow[i]; ++k) {
      int Ak = Aindices_i[k];
      double Aval = Aval_i[k];

      int* Bcol_inds = Bindices[Ak];
      double* Bvals_k = Bvalues[Ak];

      if (remote[Ak]) {
	for(j=0; j<BnumEntriesPerRow[Ak]; ++j) {
	  int loc = Bcol_inds[j] - C_firstCol_import;
	  C_row_i_import[loc] += Aval*Bvals_k[j];
	}
      }
      else {
	for(j=0; j<BnumEntriesPerRow[Ak]; ++j) {
	  int loc = Bcol_inds[j] - C_firstCol;
	  C_row_i[loc] += Aval*Bvals_k[j];
	}
      }
    }

    //
    //Now loop across the C_row_i values and put the non-zeros into C.
    //

    int global_row = A.GRID(i);

    for(j=0; j<C_numCols; ++j) {
      //If this value is zero, don't put it into the C matrix.
      if (C_row_i[j] == 0.0) continue;

      int global_col = Bcolmap->GID(C_firstCol + j);

      //Now put the C_ij quantity into the result C matrix.
      //
      //Try SumInto first, and if that returns a positive error code (meaning
      //that the location doesn't exist) then use Insert.

      int err = C.SumIntoGlobalValues(global_row, 1, &(C_row_i[j]), &global_col);
      if (err < 0) {
	return(err);
      }
      if (err > 0) {
	err = C.InsertGlobalValues(global_row, 1, &(C_row_i[j]), &global_col);
	if (err < 0) {
	  //If we jump out here, it probably means C is Filled, and doesn't
	  //have all the necessary nonzero locations.
	  return(err);
	}
      }

      C_row_i[j] = 0.0;
    }

    //Now loop across the C_row_i_import values and put the non-zeros into C.

    for(j=0; j<C_numCols_import; ++j) {
      //If this value is zero, don't put it into the C matrix.
      if (C_row_i_import[j] == 0.0) continue;

      int global_col = importBcolmap->GID(C_firstCol_import + j);

      //Now put the C_ij quantity into the result C matrix.
      //
      //Try SumInto first, and if that returns a positive error code (meaning
      //that the location doesn't exist) then use Insert.

      int err = C.SumIntoGlobalValues(global_row, 1,
				      &(C_row_i_import[j]), &global_col);
      if (err < 0) {
	return(err);
      }
      if (err > 0) {
	err = C.InsertGlobalValues(global_row, 1,
				   &(C_row_i_import[j]), &global_col);
	if (err < 0) {
	  //If we jump out here, it probably means C is Filled, and doesn't
	  //have all the necessary nonzero locations.
	  return(err);
	}
      }

      C_row_i_import[j] = 0.0;
    }
  }

  if (!C.Filled()) {
    EPETRA_CHK_ERR( C.FillComplete() );
  }

  if (AnumRows > 0) {
    delete [] AnumEntriesPerRow;
    delete [] Aindices;
    delete [] Avalues;
  }

  if (BnumRows > 0) {
    delete [] BnumEntriesPerRow;
    delete [] Bindices;
    delete [] Bvalues;
    delete [] remote;
  }

  if (importB != NULL) {
    delete importB;
  }

  delete [] C_row_i;
  delete [] C_row_i_import;

  return(0);
}

} // namespace EpetraExt

