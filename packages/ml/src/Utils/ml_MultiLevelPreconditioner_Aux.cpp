/*!
 *  \file ml_MultiLevelPreconditioner.cpp
 *
 *  \brief ML black-box preconditioner for Epetra_RowMatrix derived classes.
 *
 *  \author Marzio Sala, SNL, 9214
 *
 *  \date Last update do Doxygen: 22-Jul-04
 *
 */

#include "ml_common.h"
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "ml_memory.h"
#include "ml_DD_prec.h"
#include <iostream>
#include <iomanip>

#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//#include <cstring>
#include "ml_amesos_wrap.h"
#include "ml_agg_METIS.h"
#include "ml_epetra_utils.h"

#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_ParMETIS.h"

#include "ml_anasazi.h"

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::
CreateAuxiliaryMatrixCrs(Epetra_FECrsMatrix* &FakeMatrix)
{

  int NumMyRows = RowMatrix_->NumMyRows();

  const Epetra_Map& RowMap = RowMatrix_->RowMatrixRowMap();
  const Epetra_Map& ColMap = RowMatrix_->RowMatrixColMap();
  FakeMatrix = new Epetra_FECrsMatrix(Copy,RowMap,
				      2*RowMatrix_->MaxNumEntries());

  if (FakeMatrix == 0)
    ML_CHK_ERR(-1); // something went wrong

  int NumDimensions = 0;

  double* x_coord = List_.get("aggregation: x-coordinates", (double *)0);
  if (x_coord != 0) ++NumDimensions;

  // at least x-coordinates must be not null
  if( NumDimensions == 0 ) {
    cerr << ErrorMsg_ << "Option `aggregation: use auxiliary matrix' == true" << endl
         << ErrorMsg_ << "requires x-, y-, or z-coordinates." << endl
         << ErrorMsg_ << "You must specify them using options" << endl
         << ErrorMsg_ << "`aggregation: x-coordinates' (and equivalently for" << endl
         << ErrorMsg_ << "y- and z-." << endl;
    ML_CHK_ERR(-2); // wrong parameters
  }

  double* y_coord = List_.get("aggregation: y-coordinates", (double *)0);
  if (y_coord != 0) ++NumDimensions;

  double* z_coord = List_.get("aggregation: z-coordinates", (double *)0);
  if (z_coord != 0) ++NumDimensions;

  // small check to avoid strange behavior
  if( z_coord != 0 && y_coord == 0 ) {
    cerr << ErrorMsg_ << "Something wrong: `aggregation: y-coordinates'" << endl
         << ErrorMsg_ << "is null, while `aggregation: z-coordinates' is null" << endl;
    ML_CHK_ERR(-3); // something went wrong
  }

  double theta = List_.get("aggregation: theta",0.0);

  bool SymmetricPattern = List_.get("aggregation: use symmetric pattern",false);

  // usual crap to clutter the output
  if( verbose_ ) {
    cout << endl;
    cout << PrintMsg_ << "*** Using auxiliary matrix to create the aggregates" << endl
         << PrintMsg_ << "*** Number of dimensions = " << NumDimensions << endl
         << PrintMsg_ << "*** theta = " << theta;
    if( SymmetricPattern ) cout << ", using symmetric pattern" << endl;
    else                   cout << ", using original pattern" << endl;
    cout << endl;
  }

  // create vectors containing coordinates, replicated for all unknonws
  // FIXME: I don't really need Z in all cases
  // for west
  // I am over-allocating, for large number of equations per node
  // this is not optimal. However, it is a only-once importing
  // of some more data. It should harm too much...
  //
  // The following will work for constant number of equations per
  // node only. The fix should be easy, though.
  Epetra_Vector RowX(RowMap); RowX.PutScalar(0.0);
  Epetra_Vector RowY(RowMap); RowY.PutScalar(0.0);
  Epetra_Vector RowZ(RowMap); RowZ.PutScalar(0.0);

  for (int i = 0 ; i < NumMyRows ; i += NumPDEEqns_) {
    RowX[i] = x_coord[i / NumPDEEqns_];
    if (NumDimensions > 1) RowY[i] = y_coord[i / NumPDEEqns_];
    if (NumDimensions > 2) RowZ[i] = z_coord[i / NumPDEEqns_];
  }

  // create vectors containing coordinates for columns 
  // (this is useful only if MIS/ParMETIS are used)
  Epetra_Vector ColX(ColMap); ColX.PutScalar(0.0);
  Epetra_Vector ColY(ColMap); ColY.PutScalar(0.0);
  Epetra_Vector ColZ(ColMap); ColZ.PutScalar(0.0);
  
  // get coordinates for non-local nodes (in column map)
  Epetra_Import Importer(ColMap,RowMap);
  ColX.Import(RowX,Importer,Insert);
  if (NumDimensions > 1) ColY.Import(RowY,Importer,Insert);
  if (NumDimensions > 2) ColZ.Import(RowZ,Importer,Insert);
  
  // global row and column numbering
  int* MyGlobalRowElements = RowMap.MyGlobalElements();
  int* MyGlobalColElements = ColMap.MyGlobalElements();

  // room for getrow()
  int MaxNnz = RowMatrix_->MaxNumEntries();
  vector<int> colInd(MaxNnz);
  vector<double> colVal(MaxNnz);
  vector<double> coord_i(3);
  vector<double> coord_j(3);

  // =================== //
  // cycle over all rows //
  // =================== //

  for (int i = 0; i < NumMyRows ; i += NumPDEEqns_) {

    int GlobalRow = MyGlobalRowElements[i];

    if( i%NumPDEEqns_ == 0 ) { // do it just once for each block row
      switch( NumDimensions ) {
      case 3:
	coord_i[2] = RowZ[i];
      case 2:
	coord_i[1] = RowY[i];
      case 1:
	coord_i[0] = RowX[i];

      }

      int NumEntries;
      ML_CHK_ERR(RowMatrix_->ExtractMyRowCopy(i,MaxNnz,NumEntries,
                                              &colVal[0],&colInd[0]));

      // NOTE: for VBR matrices, the "real" value that will be used in
      // the subsequent part of the code is only the one for the first
      // equations. For each block, I replace values with the sum of
      // the abs of each block entry.
      for (int j = 0 ; j < NumEntries ; j += NumPDEEqns_) {
	colVal[j] = fabs(colVal[j]);
	for (int k = 1 ; k < NumPDEEqns_ ; ++k) {
	  colVal[j] += fabs(colVal[j+k]);
	}
      }

      // work only on the first equations. Theta will blend the
      // coordinate part with the sub of abs of row elements.
      int GlobalCol;
      double total = 0.0;

      for (int j = 0 ; j < NumEntries ; j += NumPDEEqns_) {

        if (colInd[j] >= NumMyRows)
          continue;

	if (colInd[j]%NumPDEEqns_ == 0) { 

	  // insert diagonal later
	  if (colInd[j] != i) {

	    // get coordinates of this node
	    switch (NumDimensions) {
	    case 3:
	      coord_j[2] = ColZ[colInd[j]];
	    case 2:
	      coord_j[1] = ColY[colInd[j]];
	    case 1:
	      coord_j[0] = ColX[colInd[j]];
	    }

	    // d2 is the square of the distance between node `i' and
	    // node `j'
	    double d2 = (coord_i[0] - coord_j[0]) * (coord_i[0] - coord_j[0]) +
	      (coord_i[1] - coord_j[1]) * (coord_i[1] - coord_j[1]) +		     
	      (coord_i[2] - coord_j[2]) * (coord_i[2] - coord_j[2]);

	    if (d2 == 0.0) {
	      cerr << endl;
	      cerr << ErrorMsg_ << "distance between node " << i/NumPDEEqns_ << " and node " 
                   << colInd[j]/NumPDEEqns_ << endl
                   << ErrorMsg_ << "is zero. Coordinates of these nodes are" << endl
	           << ErrorMsg_ << "x_i = " << coord_i[0] << ", x_j = " << coord_j[0] << endl  
		   << ErrorMsg_ << "y_i = " << coord_i[1] << ", y_j = " << coord_j[1] << endl  
		   << ErrorMsg_ << "z_i = " << coord_i[2] << ", z_j = " << coord_j[2] << endl  
		   << ErrorMsg_ << "Now proceeding with distance = 1.0" << endl;
	      cerr << endl;
	      d2 = 1.0;
	    }

	    // blend d2 with the actual values of the matrix
	    // FIXME: am I useful?
	    double val = -(1.0 - theta) * (1.0 / d2) + theta * (colVal[j]);

	    GlobalCol = MyGlobalColElements[colInd[j]];

	    // insert this value on all rows
	    for (int k = 0 ; k < NumPDEEqns_ ; ++k) {
	      int row = GlobalRow + k;
	      int col = GlobalCol + k;

	      if (FakeMatrix->SumIntoGlobalValues(1,&row,1,&col,&val) != 0) {
		ML_CHK_ERR(FakeMatrix->InsertGlobalValues(1,&row,1,&col,&val));
	      }

	    }

	    total -= val;

	    // put (j,i) element as well, only for in-process stuff.
	    // I have some problems with off-processor elements.
            // It is here that I need the FE matrix.
	    if (SymmetricPattern == true && colInd[j] < NumMyRows ) {

	      for( int k=0 ; k<NumPDEEqns_ ; ++k ) {
		int row = GlobalCol+k;
		int col = GlobalRow+k;

		if( FakeMatrix->SumIntoGlobalValues(1,&row,1,&col,&val) != 0 ) { 
		  ML_CHK_ERR(FakeMatrix->InsertGlobalValues(1,&row,1,&col,&val));
		}

	      }
	      total -= val;
	    }
	  } 
	}
      }

      // create lines with zero-row sum
      for (int k = 0 ; k < NumPDEEqns_ ; ++k) {
	int row = GlobalRow + k;
	if (FakeMatrix->SumIntoGlobalValues(1,&row,1,&row,&total) != 0) {
	  if (FakeMatrix->InsertGlobalValues(1,&row,1,&row,&total) != 0)
	    ML_CHK_ERR(-9); // something went wrong
	}

      }
    }
  }

  if (FakeMatrix->FillComplete())
    ML_CHK_ERR(-5); // something went wrong

  // stick pointer in Amat for level 0 (finest level)
  ml_->Amat[LevelID_[0]].data = (void *)FakeMatrix;

  // tell ML to keep the tentative prolongator
  ML_Aggregate_Set_Reuse(agg_);

  // pray that no bugs will tease us
  return(0);
}

// ============================================================================
int ML_Epetra::MultiLevelPreconditioner::
CreateAuxiliaryMatrixVbr(Epetra_VbrMatrix* &FakeMatrix)
{

  // FakeMatrix has already been created before
  if (FakeMatrix == 0)
    ML_CHK_ERR(-1); // something went wrong

  int NumDimensions = 0;

  double* x_coord = List_.get("aggregation: x-coordinates", (double *)0);
  if (x_coord != 0) ++NumDimensions;

  // at least x-coordinates must be not null
  if( NumDimensions == 0 ) {
    cerr << ErrorMsg_ << "Option `aggregation: use auxiliary matrix' == true" << endl
         << ErrorMsg_ << "requires x-, y-, or z-coordinates." << endl
         << ErrorMsg_ << "You must specify them using options" << endl
         << ErrorMsg_ << "`aggregation: x-coordinates' (and equivalently for" << endl
         << ErrorMsg_ << "y- and z-)." << endl;
    ML_CHK_ERR(-2); // wrong parameters
  }

  double* y_coord = List_.get("aggregation: y-coordinates", (double *)0);
  if (y_coord != 0) ++NumDimensions;

  double* z_coord = List_.get("aggregation: z-coordinates", (double *)0);
  if (z_coord != 0) ++NumDimensions;

  // small check to avoid strange behavior
  if( z_coord != 0 && y_coord == 0 ) {
    cerr << ErrorMsg_ << "Something wrong: `aggregation: y-coordinates'" << endl
         << ErrorMsg_ << "is null, while `aggregation: z-coordinates' is not null" << endl;
    ML_CHK_ERR(-3); // something went wrong
  }

  // usual crap to clutter the output
  if( verbose_ ) {
    cout << endl;
    cout << PrintMsg_ << "*** Using auxiliary matrix to create the aggregates" << endl
         << PrintMsg_ << "*** Number of dimensions = " << NumDimensions << endl
         << PrintMsg_ << "*** (the version for Epetra_VbrMatrix is currently used)" << endl;
    cout << endl;
  }

  const Epetra_BlockMap& RowMap = FakeMatrix->RowMap();
  const Epetra_BlockMap& ColMap = FakeMatrix->ColMap();
  int NumMyRowElements = RowMap.NumMyElements();
  int NumMyColElements = ColMap.NumMyElements();
  int* MyGlobalRowElements = RowMap.MyGlobalElements();
  int* MyGlobalColElements = ColMap.MyGlobalElements();

  // use point map to exchange coordinates
  Epetra_Map PointRowMap(-1,NumMyRowElements,MyGlobalRowElements,0,Comm());
  Epetra_Map PointColMap(-1,NumMyColElements,MyGlobalColElements,0,Comm());
  
  Epetra_Vector RowX(PointRowMap);
  Epetra_Vector RowY(PointRowMap);
  Epetra_Vector RowZ(PointRowMap);

  for (int i = 0 ; i < NumMyRowElements ; ++i) {
    RowX[i] = x_coord[i];
    if (NumDimensions > 1) RowY[i] = y_coord[i];
    if (NumDimensions > 2) RowZ[i] = z_coord[i];
  }

  // create vectors containing coordinates for columns 
  // (this is useful only if MIS/ParMETIS are used)
  Epetra_Vector ColX(PointColMap);
  Epetra_Vector ColY(PointColMap);
  Epetra_Vector ColZ(PointColMap); 
  
  // get coordinates for non-local nodes (in column map)
  Epetra_Import Importer(PointColMap,PointRowMap);
  ColX.Import(RowX,Importer,Insert);
  if (NumDimensions > 1) ColY.Import(RowY,Importer,Insert);
  if (NumDimensions > 2) ColZ.Import(RowZ,Importer,Insert);

  // room for getrow()
  int MaxNnz = FakeMatrix->MaxNumEntries();
  vector<int> colInd(MaxNnz);
  vector<double> colVal(MaxNnz);
  vector<double> coord_i(3);
  vector<double> coord_j(3);

  // change the entries of FakeMatrix so that it corresponds to a discrete
  // Laplacian. Note: This is not exactly the same as in the Crs case.

  FakeMatrix->PutScalar(0.0);

  for (int LocalRow = 0; LocalRow < NumMyRowElements ; ++LocalRow) {

    int RowDim, NumBlockEntries;
    int* BlockIndices;
    Epetra_SerialDenseMatrix** RowValues;

    switch (NumDimensions) {
    case 3:
      coord_i[2] = RowZ[LocalRow];
    case 2:
      coord_i[1] = RowY[LocalRow];
    case 1:
      coord_i[0] = RowX[LocalRow];
    }

    FakeMatrix->ExtractMyBlockRowView(LocalRow,RowDim,NumBlockEntries,
                                      BlockIndices,RowValues);

    // accumulator for zero row-sum
    double total = 0.0;

    for (int j = 0 ; j < NumBlockEntries ; ++j) {

      int LocalCol = BlockIndices[j];

      // insert diagonal later
      if (LocalCol != LocalRow) {

        // get coordinates of this node
        switch (NumDimensions) {
        case 3:
          coord_j[2] = ColZ[LocalCol];
        case 2:
          coord_j[1] = ColY[LocalCol];
        case 1:
          coord_j[0] = ColX[LocalCol];
        }

        // d2 is the square of the distance between node `i' and
        // node `j'
        double d2 = (coord_i[0] - coord_j[0]) * (coord_i[0] - coord_j[0]) +
          (coord_i[1] - coord_j[1]) * (coord_i[1] - coord_j[1]) +
          (coord_i[2] - coord_j[2]) * (coord_i[2] - coord_j[2]);

        if (d2 == 0.0) {
          cerr << endl;
          cerr << ErrorMsg_ << "distance between node " << LocalRow << " and node " 
            << LocalCol << endl
            << ErrorMsg_ << "is zero. Coordinates of these nodes are" << endl
            << ErrorMsg_ << "x_i = " << coord_i[0] << ", x_j = " << coord_j[0] << endl  
            << ErrorMsg_ << "y_i = " << coord_i[1] << ", y_j = " << coord_j[1] << endl  
            << ErrorMsg_ << "z_i = " << coord_i[2] << ", z_j = " << coord_j[2] << endl  
            << ErrorMsg_ << "Now proceeding with distance = 1.0" << endl;
          cerr << endl;
          d2 = 1.0;
        }

        for (int k = 0 ; k < RowValues[j]->M() ; ++k) {
          for (int h = 0 ; h < RowValues[j]->N() ; ++h) {
            if (k == h)
              (*RowValues[j])(k,h) = - 1.0 / d2;
          }
        }

        total += 1.0 /d2;

      } 
    }

    // check that the diagonal block exists
    bool ok = false;
    int DiagonalBlock = 0;
    for (int j = 0 ; j < NumBlockEntries ; ++j) {
      if (BlockIndices[j] == LocalRow) {
        DiagonalBlock = j;
        ok = true;
        break;
      }
    }
    assert (ok == true);

    for (int k = 0 ; k < RowValues[DiagonalBlock]->N() ; ++k) {
      (*RowValues[DiagonalBlock])(k,k) = total;
    }
  }

  // tell ML to keep the tentative prolongator
  ML_Aggregate_Set_Reuse(agg_);

  // pray that no bugs will tease us
  return(0);
}

#endif /*ifdef HAVE_ML_EPETRA && HAVE_ML_TEUCHOS */
