// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER

#ifndef TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP
#define TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP

#include "Tpetra_Export.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef DOXYGEN_USE_ONLY
  // #include "Tpetra_RowMatrixtransposer_decl.hpp"
#endif

namespace Tpetra {

template<class Scalar, 
	 class LocalOrdinal, 
	 class GlobalOrdinal, 
	 class Node, 
	 class SpMatOps>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::
RowMatrixTransposer (const Teuchos::RCP<const crs_matrix_type>& origMatrix)
  : origMatrix_ (origMatrix) {}

template<class Scalar, 
	 class LocalOrdinal, 
	 class GlobalOrdinal, 
	 class Node, 
	 class SpMatOps>
TEUCHOS_DEPRECATED
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::
RowMatrixTransposer (const crs_matrix_type& origMatrix)
  : origMatrix_ (Teuchos::rcpFromRef (origMatrix)) {}

// mfh 03 Feb 2013: In a definition outside the class like this, the
// return value is considered outside the class scope (for things like
// resolving typedefs), but the arguments are considered inside the
// class scope.
template<class Scalar, 
	 class LocalOrdinal,
	 class GlobalOrdinal, 
	 class Node, 
	 class SpMatOps>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> >
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::
createTranspose (const OptimizeOption optimizeTranspose,
		 Teuchos::RCP<const map_type> transposeRowMap)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;

  //
  // This transpose is based upon the approach in EpetraExt.
  //

  // mfh 03 Feb 2013: The domain Map of the input matrix will become
  // the range Map of the transpose, so it's a good default choice for
  // the row Map of the transpose.
  RCP<const map_type> newRowMap = transposeRowMap.is_null () ? 
    origMatrix_->getDomainMap () : transposeRowMap;

  global_size_t numLocalCols = origMatrix_->getNodeNumCols();
  global_size_t numLocalRows = origMatrix_->getNodeNumRows();
  ArrayView<const LO> localIndices;
  ArrayView<const Scalar> localValues;
  // Determine how many nonzeros there are per row in the transpose.
  ArrayRCP<size_t> TransNumNz(numLocalCols,0);
  for (size_t i=0; i<numLocalRows; ++i) {
    const size_t numEntriesInRow = origMatrix_->getNumEntriesInLocalRow(i);
    origMatrix_->getLocalRowView(i, localIndices, localValues);
    for (size_t j=0; j<numEntriesInRow; ++j) {
      ++TransNumNz[ localIndices[j] ];
    }
  }

  //create temporary row-major storage for the transposed matrix
  size_t numLocalNnz  = origMatrix_->getNodeNumEntries();
  ArrayRCP<GO> TransIndices(numLocalNnz);
  ArrayRCP<Scalar> TransValues(numLocalNnz);

  //ptr is an vector of indices into TransIndices and TransValues
  //ptr[i]..ptr[i+1]-1 correspond to nonzero rows/values in the original matrix's local column i 
  std::vector<size_t> ptr(numLocalCols+1);
  ptr[0]=0;
  for (size_t i=1; i<numLocalCols+1; ++i) {
    ptr[i] = ptr[i-1] + TransNumNz[i-1];
  }

  //TransNumNz will now be used to track how many entries have been placed thus far in each
  //row of the row-major storage.
  for (size_t i=0; i<numLocalCols; ++i) TransNumNz[i] = 0;

  //populate the row-major storage so that the data for the transposed matrix is easy to access
  const map_type& origRowMap = *(origMatrix_->getRowMap());
  for (size_t i=0; i<numLocalRows; ++i) {

    const size_t numEntriesInRow = origMatrix_->getNumEntriesInLocalRow(i);
    origMatrix_->getLocalRowView(i, localIndices, localValues);

    const GO transCol = origRowMap.getGlobalElement(i);
    for (size_t j=0; j<numEntriesInRow; ++j) {
      const LO transRow = localIndices[j];
      const LO offset = TransNumNz[transRow];
      TransIndices[ptr[transRow]+offset] = transCol;
      TransValues [ptr[transRow]+offset] = localValues[j];
      ++TransNumNz[transRow];
    }
    
  } //for (size_t i=0; i<numLocalRows; ++i)

  //Allocate and populate temporary matrix with rows not uniquely owned
  RCP<const map_type> transMap = origMatrix_->getColMap();
  RCP<crs_matrix_type> transMatrixWithSharedRows(new crs_matrix_type (transMap, TransNumNz, StaticProfile));

  const map_type& transMapRef = *transMap;
  for (size_t i=0; i<numLocalCols; ++i) {
     const size_t leng=ptr[i+1]-ptr[i];
     transMatrixWithSharedRows->insertGlobalValues(transMapRef.getGlobalElement(i),
                                                   TransIndices.view(ptr[i],leng),
                                                   TransValues.view(ptr[i],leng) );
  }

  //A comment in EpetraExt says
  // "Note: The following call to FillComplete is currently necessary because
  //        some global constants that are needed by the Export () are computed in this routine"
  //TODO determine whether this is necessary here
  RCP<ParameterList> params = parameterList ();
  params->set ("Optimize Storage", false); //TODO should storage be optimized for this temporary matrix? EpetraExt does not.
  transMatrixWithSharedRows->fillComplete(origMatrix_->getRangeMap(), origMatrix_->getDomainMap(), params);

  //exporter that is used in the transfers of nnz per row and rows themselves
  RCP<Tpetra::Export<LocalOrdinal,GlobalOrdinal,Node> > exporter = rcp( new Tpetra::Export<LocalOrdinal,GlobalOrdinal,Node>(transMap,newRowMap) );

  //RCP<Vector<size_t, LO, GO, Node> > partialNnzPerRow = rcp(new Vector<size_t,LO,GO,Node>(transMap,TransNumNz()));
  //RCP<Vector<size_t, LO, GO, Node> > fullNnzPerRow = rcp(new Vector<size_t,LO,GO,Node>(newRowMap,false));
  //fullNnzPerRow->doExport(*partialNnzPerRow,*exporter,Tpetra::ADD);
  //const ArrayRCP<const size_t> nnzPerRow = fullNnzPerRow->getData();
  //RCP<crs_matrix_type> transposeMatrix(new crs_matrix_type (newRowMap, nnzPerRow, StaticProfile));
  
  // The following code is to avoid compilation problems with explicit instantiation.
  // It replaces the 5 line of code abovve.  Here's the issue:
  // The CRS matrix ctor requires size_t for specifying nz per row.  But to avoid required EI of Vector<size_t>,
  // I copy to a Vector<LO>, communicate it, then copy it back to a Vector<size_t>.
  // The overhead is allocating and filling two additional ArrayRCPs, each with global length about the #rows in the matrix.
  ArrayRCP<LO> TransNumNzAsLO(TransNumNz.size());
  for (LO i=0; i<TransNumNz.size(); ++i) TransNumNzAsLO[i] = Teuchos::as<LO>(TransNumNz[i]);
  RCP<Vector<LO, LO, GO, Node> > partialNnzPerRow = rcp(new Vector<LO,LO,GO,Node>(transMap,TransNumNzAsLO()));
  RCP<Vector<LO, LO, GO, Node> > fullNnzPerRow = rcp(new Vector<LO,LO,GO,Node>(newRowMap,false));
  fullNnzPerRow->doExport(*partialNnzPerRow,*exporter,Tpetra::ADD);
  const ArrayRCP<const LO> nnzPerRow = fullNnzPerRow->getData();
  ArrayRCP<size_t> nnzPerRowAsSizeT(nnzPerRow.size());
  for (LO i=0; i<nnzPerRowAsSizeT.size(); ++i) nnzPerRowAsSizeT[i] = Teuchos::as<size_t>(nnzPerRow[i]);

  RCP<crs_matrix_type> transposeMatrix(new crs_matrix_type (newRowMap, nnzPerRowAsSizeT, StaticProfile));

  transposeMatrix->doExport(*transMatrixWithSharedRows,*exporter,Tpetra::ADD);

  const bool optimizeStorage = (optimizeTranspose == DoOptimizeStorage);
  params->set ("Optimize Storage", optimizeStorage);
  transposeMatrix->fillComplete(origMatrix_->getRangeMap(), origMatrix_->getDomainMap(), params);
  return transposeMatrix;


}
//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIXTRANSPOSE_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class RowMatrixTransposer< SCALAR, LO , GO , NODE >;


}

#endif
