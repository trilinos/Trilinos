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
#include "Tpetra_Import.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef DOXYGEN_USE_ONLY
  // #include "Tpetra_RowMatrixtransposer_decl.hpp"
#endif

namespace Tpetra {

template<class Scalar,
     class LocalOrdinal,
     class GlobalOrdinal,
     class Node>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
RowMatrixTransposer (const Teuchos::RCP<const crs_matrix_type>& origMatrix)
  : origMatrix_ (origMatrix) {}

template<class Scalar,
     class LocalOrdinal,
     class GlobalOrdinal,
     class Node>
TEUCHOS_DEPRECATED
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
RowMatrixTransposer (const crs_matrix_type& origMatrix)
  : origMatrix_ (Teuchos::rcpFromRef (origMatrix)) {}

template<class Scalar,
     class LocalOrdinal,
     class GlobalOrdinal,
     class Node>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createTranspose()
{
  using Teuchos::RCP;

  // Do the local transpose
  RCP<crs_matrix_type> transMatrixWithSharedRows = createTransposeLocal ();

  // If transMatrixWithSharedRows has an exporter, that's what we
  // want.  If it doesn't, the rows aren't actually shared, and we're
  // done!
  RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > exporter =
    transMatrixWithSharedRows->getGraph ()->getExporter ();
  if (exporter.is_null ()) {
    return transMatrixWithSharedRows;
  }
  else {
    // Use the Export object to do a fused Export and fillComplete.
    return exportAndFillCompleteCrsMatrix<crs_matrix_type> (transMatrixWithSharedRows, *exporter);
  }
}


// mfh 03 Feb 2013: In a definition outside the class like this, the
// return value is considered outside the class scope (for things like
// resolving typedefs), but the arguments are considered inside the
// class scope.
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
createTransposeLocal ()
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Import<LO, GO, Node> import_type;
  typedef Tpetra::Export<LO, GO, Node> export_type;

  //
  // This transpose is based upon the approach in EpetraExt.
  //
  size_t numLocalCols = origMatrix_->getNodeNumCols();
  size_t numLocalRows = origMatrix_->getNodeNumRows();
  size_t numLocalNnz  = origMatrix_->getNodeNumEntries();

  // Determine how many nonzeros there are per row in the transpose.
  Array<size_t> CurrentStart(numLocalCols,0);
  ArrayView<const LO> localIndices;
  ArrayView<const Scalar> localValues;
  RCP<const crs_matrix_type> crsMatrix =
    rcp_dynamic_cast<const crs_matrix_type> (origMatrix_);
  if (crsMatrix == Teuchos::null) {
    for (size_t i=0; i<numLocalRows; ++i) {
      const size_t numEntriesInRow = origMatrix_->getNumEntriesInLocalRow(i);
      origMatrix_->getLocalRowView(i, localIndices, localValues);
      for (size_t j=0; j<numEntriesInRow; ++j) {
        ++CurrentStart[ localIndices[j] ];
      }
    }
  } else {
    ArrayRCP<const size_t> origRowPtr_rcp;
    ArrayRCP<const LO>     origColInd_rcp;
    ArrayRCP<const Scalar> origValues_rcp;
    crsMatrix->getAllValues(origRowPtr_rcp, origColInd_rcp, origValues_rcp);
    ArrayView<const LO> origColInd = origColInd_rcp();
    for (LO j=0; j<origColInd.size(); ++j) {
      ++CurrentStart[ origColInd[j] ];
    }
  }

  // create temporary row-major storage for the transposed matrix

  ArrayRCP<size_t> rowptr_rcp(numLocalCols+1);
  ArrayRCP<LO>     colind_rcp(numLocalNnz);
  ArrayRCP<Scalar> values_rcp(numLocalNnz);

  // Since ArrayRCP's are slow...
  ArrayView<size_t> TransRowptr = rowptr_rcp();
  ArrayView<LO>     TransColind = colind_rcp();
  ArrayView<Scalar> TransValues = values_rcp();

  // Scansum the TransRowptr; reset CurrentStart
  TransRowptr[0]=0;
  for (size_t i=1; i<numLocalCols+1; ++i) TransRowptr[i]  = CurrentStart[i-1] + TransRowptr[i-1];
  for (size_t i=0; i<numLocalCols;   ++i) CurrentStart[i] = TransRowptr[i];

  // populate the row-major storage so that the data for the transposed
  // matrix is easy to access
  if (crsMatrix == Teuchos::null) {
    for (size_t i=0; i<numLocalRows; ++i) {
      const size_t numEntriesInRow = origMatrix_->getNumEntriesInLocalRow (i);
      origMatrix_->getLocalRowView(i, localIndices, localValues);

      for (size_t j=0; j<numEntriesInRow; ++j) {
        size_t idx = CurrentStart[localIndices[j]];
        TransColind[idx] = Teuchos::as<LO>(i);
        TransValues[idx] = localValues[j];
        ++CurrentStart[localIndices[j]];
      }
    } //for (size_t i=0; i<numLocalRows; ++i)
  } else {
    ArrayRCP<const size_t> origRowPtr_rcp;
    ArrayRCP<const LO>     origColInd_rcp;
    ArrayRCP<const Scalar> origValues_rcp;
    crsMatrix->getAllValues(origRowPtr_rcp, origColInd_rcp, origValues_rcp);
    ArrayView<const size_t>   origRowPtr = origRowPtr_rcp();
    ArrayView<const LO> origColInd = origColInd_rcp();
    ArrayView<const Scalar>   origValues = origValues_rcp();
    size_t k=0;
    for (LO i=0; i<origRowPtr.size()-1; ++i) {
      const LO rowIndex = Teuchos::as<LO>(i);
      for (size_t j=origRowPtr[i]; j<origRowPtr[i+1]; ++j) {
        size_t idx = CurrentStart[origColInd[k]];
        TransColind[idx] = rowIndex;
        TransValues[idx] = origValues[k];
        ++CurrentStart[origColInd[k++]];
      }
    }
  }

  //Allocate and populate temporary matrix with rows not uniquely owned
  RCP<crs_matrix_type> transMatrixWithSharedRows =
    rcp (new crs_matrix_type (origMatrix_->getColMap (),
                              origMatrix_->getRowMap (), 0));
  transMatrixWithSharedRows->setAllValues (rowptr_rcp, colind_rcp, values_rcp);

  // Prebuild the importers and exporters the no-communication way,
  // flipping the importers and exporters around.
  RCP<const import_type> myImport;
  RCP<const export_type> myExport;
  if (! origMatrix_->getGraph ()->getImporter ().is_null ()) {
    myExport = rcp (new export_type (*origMatrix_->getGraph ()->getImporter ()));
  }
  if (! origMatrix_->getGraph ()->getExporter ().is_null ()) {
    myImport = rcp (new import_type (*origMatrix_->getGraph ()->getExporter ()));
  }

  // Call ESFC & return
  transMatrixWithSharedRows->expertStaticFillComplete (origMatrix_->getRangeMap (),
                                                       origMatrix_->getDomainMap (),
                                                       myImport, myExport);
  return transMatrixWithSharedRows;
}
//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIXTRANSPOSER_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class RowMatrixTransposer< SCALAR, LO , GO , NODE >;


}

#endif
