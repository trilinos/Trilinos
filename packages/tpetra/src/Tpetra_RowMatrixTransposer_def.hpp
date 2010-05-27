
//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Export.hpp"
#ifdef DOXYGEN_USE_ONLY
  // #include "Tpetra_RowMatrixTransposer_decl.hpp"
#endif

namespace Tpetra{
	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RowMatrixTransposer(const Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > OrigMatrix)
  : OrigMatrix_(OrigMatrix),
    TransposeMatrix_(0),
    TransposeExporter_(0),
    TransposeRowMap_(0),
    TransposeCreated_(false),
    MakeDataContiguous_(false),
    NumMyRows_(0),
    NumMyCols_(0),
    MaxNumEntries_(0),
    //Indices_(NULL),
    //Values_(NULL),
    //TransNumNz_(NULL),
    //TransIndices_(NULL),
    //TransValues_(NULL),
    //TransMyGlobalEquations_(NULL),
    OrigMatrixIsCrsMatrix_(false)
{
}
//=============================================================================
/*	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RowMatrixTransposer(const RowMatrixTransposer& Source)
  :OrigMatrix_(Source.OrigMatrix_),
   TransposeMatrix_(0),
   TransposeExporter_(0),
   TransposeRowMap_(0),
   TransposeCreated_(Source.TransposeCreated_),
   MakeDataContiguous_(Source.MakeDataContiguous_),
   NumMyRows_(0),
   NumMyCols_(0),
   MaxNumEntries_(0),
   //Indices_(NULL),
   //Values_(NULL),
   //TransNumNz_(NULL),
   //TransIndices_(NULL),
   //TransValues_(NULL),
   //TransMyGlobalEquations_(NULL),
   OrigMatrixIsCrsMatrix_(false)
{
  TransposeMatrix_ = Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*Source.TransposeMatrix_));
  if (MakeDataContiguous_) TransposeMatrix_->MakeDataContiguous();
  TransposeExporter_ = Teuchos::RCP<Export<LocalOrdinal, GlobalOrdinal, Node> >(new Export<LocalOrdinal, GlobalOrdinal, Node>(*Source.TransposeExporter_));
}*/
//=========================================================================
	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~RowMatrixTransposer(){

  DeleteData();

}

//=========================================================================
	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeleteData (){

//  size_t i;

  //if (TransposeMatrix_!=0) {delete TransposeMatrix_; TransposeMatrix_=0;}
  //if (TransposeExporter_!=0) {delete TransposeExporter_; TransposeExporter_=0;}

  // Delete any intermediate storage

  /*if (!OrigMatrixIsCrsMatrix_) {
    delete [] Indices_;
    delete [] Values_;
  }*/
	
	
  /*for(i=0; i<NumMyCols_; i++) {
    int NumIndices = TransNumNz_[i];
    if (NumIndices>0) {
      delete [] TransIndices_[i];
      delete [] TransValues_[i];
    }
  }*/
  //delete [] TransNumNz_;
  //delete [] TransIndices_;
  //delete [] TransValues_;
  //delete [] TransMyGlobalEquations_;
}

//=========================================================================
	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateTranspose (const bool MakeDataContiguous, 
						 Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > TransposeMatrix){
						 //Tpetra_Map * TransposeRowMap_in) {

  size_t i, j;

  if (TransposeCreated_) DeleteData(); // Get rid of existing data first

  //if (TransposeRowMap_in==0)
    //TransposeRowMap_ = (Map<LocalOrdinal, GlobalOrdinal, Node> *) &(OrigMatrix_->getDomainMap()); // Should be replaced with refcount =
  /*else
    TransposeRowMap_ = TransposeRowMap_in; */
  //TransposeRowMap_ = OrigMatrix_->getDomainMap(); 

  // This routine will work for any RowMatrix object, but will attempt cast the matrix to a CrsMatrix if
  // possible (because we can then use a View of the matrix and graph, which is much cheaper).

  // First get the local indices to count how many nonzeros will be in the 
  // transpose graph on each processor


  //Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > OrigCrsMatrix = Teuchos::rcp_dynamic_cast<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(OrigMatrix_);

  //OrigMatrixIsCrsMatrix_ = (OrigCrsMatrix.get()!=0); // If this pointer is non-zero, the cast to CrsMatrix worked

  NumMyRows_ = OrigMatrix_->getNodeNumRows();
  NumMyCols_ = OrigMatrix_->getNodeNumCols();
  //NumMyRows_ = OrigMatrix_->NumMyRows();
  //TransNumNz_ = new size_t[NumMyCols_];
  TransNumNz_ = Teuchos::ArrayRCP<size_t>(NumMyCols_);
  //TransIndices_ = new LocalOrdinal*[NumMyCols_];
  TransIndices_ = Teuchos::Array<Teuchos::ArrayRCP<LocalOrdinal> >(NumMyCols_);
  //TransValues_ = new Scalar*[NumMyCols_];
  TransValues_ = Teuchos::Array<Teuchos::ArrayRCP<Scalar> >(NumMyCols_);


  //size_t NumIndices;

/*  if (OrigMatrixIsCrsMatrix_) {


    Teuchos::RCP<const CrsGraph<Scalar, LocalOrdinal, GlobalOrdinal, Node> > OrigGraph = OrigCrsMatrix->Graph(); // Get matrix graph

    for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0;
    for (i=0; i<NumMyRows_; i++) {
      OrigGraph->ExtractMyRowView(i, Indices_, NumIndices); // Get view of ith row
      for (j=0; j<NumIndices; j++) ++TransNumNz_[Indices_[j]];
    }

  }*/
  //else { // OrigMatrix is not a CrsMatrix

	/*
    MaxNumEntries_ = 0;
    size_t NumEntries;
    for (i=0; i<NumMyRows_; i++) {
      OrigMatrix_->NumMyRowEntries(i, NumEntries);
      MaxNumEntries_ = EPETRA_MAX(MaxNumEntries_, NumEntries);
    }*/
	//MaxNumEntries_ = getNodeMaxNumRowEntries();
    //Indices_ = new LocalOrdinal[MaxNumEntries_];
    //Values_ = new Scalar[MaxNumEntries_];

    for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0;
    for (i=0; i<NumMyRows_; i++) {
      // Get ith row
      //EPETRA_CHK_ERR(OrigMatrix_->ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_)); 
	  OrigMatrix_->getLocalRowView(i, Indices_, Values_);
      //for (j=0; j<NumIndices; j++) ++TransNumNz_[Indices_[j]];
      for (j=0; j<OrigMatrix_->getNumEntriesInLocalRow(i); j++) ++TransNumNz_[Indices_[j]];
    }
  //}


  // Most of remaining code is common to both cases

  for(i=0; i<NumMyCols_; i++) {
    size_t NumIndices = TransNumNz_[i];
    if (NumIndices>0) {
      TransIndices_[i] = Teuchos::ArrayRCP<LocalOrdinal>(NumIndices);
      TransValues_[i] = Teuchos::ArrayRCP<Scalar>(NumIndices);
    }
  }

  // Now copy values and global indices into newly create transpose storage

  for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0; // Reset transpose NumNz counter
  for (i=0; i<NumMyRows_; i++) {
    /*if (OrigMatrixIsCrsMatrix_) {
      //EPETRA_CHK_ERR(OrigCrsMatrix->ExtractMyRowView(i, NumIndices, Values_, Indices_));
    }
    else {
      //EPETRA_CHK_ERR(OrigMatrix_->ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_));
    }*/
	OrigMatrix_->getLocalRowView(i, Indices_, Values_);

    //int ii = OrigMatrix_->RowMatrixRowMap().GID(i);
    GlobalOrdinal ii = OrigMatrix_->getRowMap()->getGlobalElement(i);
    for (j=0; j<OrigMatrix_->getNumEntriesInLocalRow(i); j++) {
      LocalOrdinal TransRow = Indices_[j];
      size_t loc = TransNumNz_[TransRow];
      TransIndices_[TransRow][loc] = ii;
      TransValues_[TransRow][loc] = Values_[j];
      ++TransNumNz_[TransRow]; // increment counter into current transpose row
    }
  }

  //  Build Transpose matrix with some rows being shared across processors.
  //  We will use a view here since the matrix will not be used for anything else

  const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TransMap = OrigMatrix_->getColMap();

  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TempTransA1(TransMap, TransNumNz_);
  //TransMyGlobalEquations_ = new int[NumMyCols_];
  TransMyGlobalEquations_ = TransMap->getNodeElementList();
  
  /* Add  rows one-at-a-time */

  for (i=0; i<NumMyCols_; i++)
    {
   //   EPETRA_CHK_ERR(TempTransA1.InsertGlobalValues(TransMyGlobalEquations_[i], 
	//					    TransNumNz_[i], TransValues_[i], TransIndices_[i]));
//	TempTransA1.insertGlobalValues(TransMyGlobalEquations_[i], Teuchos::ArrayView<LocalOrdinal>(TransIndices_[i](), TransNumNz_[i]), Teuchos::ArrayView<Scalar>(TransValues_[i](),TransNumNz_[i])); 
	TempTransA1.insertGlobalValues(TransMyGlobalEquations_[i], TransIndices_[i](), TransValues_[i]()); 
    }
 
  // Note: The following call to FillComplete is currently necessary because
  //       some global constants that are needed by the Export () are computed in this routine

  const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &domain_map = OrigMatrix_->getDomainMap();
  const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &range_map = OrigMatrix_->getRangeMap();

  TempTransA1.fillComplete(range_map, domain_map, DoNotOptimizeStorage);

  // Now that transpose matrix with shared rows is entered, create a new matrix that will
  // get the transpose with uniquely owned rows (using the same row distribution as A).

  //TransposeMatrix_ = Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(TransposeRowMap_,0));
  TransposeMatrix_ = Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(OrigMatrix_->getDomainMap(),0));

  // Create an Export object that will move TempTransA around

  //TransposeExporter_ = Teuchos::RCP<Export<LocalOrdinal, GlobalOrdinal, Node> >(new Export(TransMap, *TransposeRowMap_));

  //EPETRA_CHK_ERR(TransposeMatrix_->Export(TempTransA1, *TransposeExporter_, Add));
  //TransposeMatrix_->doExport(TempTransA1, *TransposeExporter_, Add);
  
  //EPETRA_CHK_ERR(TransposeMatrix_->FillComplete(range_map, domain_map));
  //TransposeMatrix_->FillComplete(range_map, domain_map);

  if (MakeDataContiguous) {
    //EPETRA_CHK_ERR(TransposeMatrix_->MakeDataContiguous());
    TransposeMatrix_->fillComplete(DoOptimizeStorage);
  }
  else{
    TransposeMatrix_->fillComplete(DoNotOptimizeStorage);
  }

  TransposeMatrix = TransposeMatrix_;
  TransposeCreated_ = true;

  return(0);
}
//=========================================================================
/*int Tpetra_RowMatrixTransposer::UpdateTransposeValues(Tpetra_RowMatrix * MatrixWithNewValues){

  int i, j, NumIndices;

  if (!TransposeCreated_) EPETRA_CHK_ERR(-1); // Transpose must be already created

  // Sanity check of incoming matrix.  Perform some tests to see if it is compatible with original input matrix
  if (OrigMatrix_!=MatrixWithNewValues) { // Check if pointer of new matrix is same as previous input matrix
    OrigMatrix_ = MatrixWithNewValues; // Reset this pointer if not, then check for other attributes
    if (NumMyRows_ != OrigMatrix_->NumMyRows() ||
	NumMyCols_ != OrigMatrix_->NumMyCols() ||
	NumMyRows_ != OrigMatrix_->NumMyRows()) {
      EPETRA_CHK_ERR(-2); // New matrix not compatible with previous
    }
  }

  Tpetra_CrsMatrix * OrigCrsMatrix = dynamic_cast<Tpetra_CrsMatrix *>(MatrixWithNewValues);

	
  OrigMatrixIsCrsMatrix_ = (OrigCrsMatrix!=0); // If this pointer is non-zero, the cast to CrsMatrix worked


  // Now copy values and global indices into newly create transpose storage

  for (i=0;i<NumMyCols_; i++) TransNumNz_[i] = 0; // Reset transpose NumNz counter
  for (i=0; i<NumMyRows_; i++) {
    if (OrigMatrixIsCrsMatrix_) {
      EPETRA_CHK_ERR(OrigCrsMatrix->ExtractMyRowView(i, NumIndices, Values_, Indices_));
    }
    else {
      EPETRA_CHK_ERR(OrigMatrix_->ExtractMyRowCopy(i, MaxNumEntries_, NumIndices, Values_, Indices_));
    }

    int ii = OrigMatrix_->RowMatrixRowMap().GID(i);
    for (j=0; j<NumIndices; j++) {
      int TransRow = Indices_[j];
      int loc = TransNumNz_[TransRow];
      TransIndices_[TransRow][loc] = ii;
      TransValues_[TransRow][loc] = Values_[j];
      ++TransNumNz_[TransRow]; // increment counter into current transpose row
    }
  }

  //  Build Transpose matrix with some rows being shared across processors.
  //  We will use a view here since the matrix will not be used for anything else

  const Tpetra_Map & TransMap = OrigMatrix_->RowMatrixColMap();

  Tpetra_CrsMatrix TempTransA1(View, TransMap, TransNumNz_);
  TransMap.MyGlobalElements(TransMyGlobalEquations_);*/
  
  /* Add  rows one-at-a-time */

  /*for (i=0; i<NumMyCols_; i++)
    {
      EPETRA_CHK_ERR(TempTransA1.InsertGlobalValues(TransMyGlobalEquations_[i], 
						    TransNumNz_[i], TransValues_[i], TransIndices_[i]));
    }
 
  // Note: The following call to FillComplete is currently necessary because
  //       some global constants that are needed by the Export () are computed in this routine
  const Tpetra_Map& domain_map = OrigMatrix_->OperatorDomainMap();
  const Tpetra_Map& range_map = OrigMatrix_->OperatorRangeMap();

  EPETRA_CHK_ERR(TempTransA1.FillComplete(range_map, domain_map, false));

  // Now that transpose matrix with shared rows is entered, update values of target transpose matrix

  TransposeMatrix_->PutScalar(0.0);  // Zero out all values of the matrix

  EPETRA_CHK_ERR(TransposeMatrix_->Export(TempTransA1, *TransposeExporter_, Add));

  return(0);
}*/

	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator=(const RowMatrixTransposer& src)
{
  (void)src;//prevents unused variable compiler warning

  //not currently supported
  bool throw_error = true;
  if (throw_error) {
    std::cerr << std::endl
	      << "RowMatrixTransposer::operator= not supported."
	      <<std::endl;
    throw -1;
  }

  return(*this);
}
//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIXTRANSPOSE_INSTANT(LO,GO,NODE) \
  \
  template class RowMatrixTransposer< LO , GO , NODE >;


}

