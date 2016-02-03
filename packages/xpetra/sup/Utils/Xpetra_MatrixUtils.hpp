// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_XPETRA_SUP_MATRIX_UTILS_HPP_
#define PACKAGES_XPETRA_SUP_MATRIX_UTILS_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"

namespace Xpetra {

/*!
  @class MatrixUtils
  @brief Xpetra utility class for common matrix-related routines

  The routines should be independent from Epetra/Tpetra and be purely implemented in Xpetra.
  Other matrix-related routines are out-sourced into other helper classes (e.g. MatrixMatrix for
  MM multiplication and addition).

*/
template <class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node>
class MatrixUtils {
#undef XPETRA_MATRIXUTILS_SHORT
#include "Xpetra_UseShortNames.hpp"

public:

  /** Given CrsMatrix objects A, B and C, form the product C = A*B.
    In a parallel setting, A and B need not have matching distributions,
    but C needs to have the same row-map as A (if transposeA is false).
    At this time C=AT*B and C=A*BT are known to not work. However,
    C=A*B and C=AT*BT are known to work, Kurtis Nusbaum 03/24/2011

    @param input Input matrix, must already have had 'FillComplete()' called.
    @param rangeMapExtractor MapExtractor object describing the splitting of rows of the output block matrix
    @param domainMapExtractor MapExtractor object describing the splitting of columns of the output block matrix
  */
  static Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > SplitMatrix(
                       const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& input,
                       Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rangeMapExtractor,
                       Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domainMapExtractor) {
    size_t numRows  = rangeMapExtractor->NumMaps();
    size_t numCols  = domainMapExtractor->NumMaps();

    RCP<const Map> fullRangeMap  = rangeMapExtractor->getFullMap();
    RCP<const Map> fullDomainMap = domainMapExtractor->getFullMap();

    /*TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap.getMaxAllGlobalIndex() != input.getRowMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap.getGlobalNumElements() != input.getRowMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullRangeMap.getNodeNumElements() != input.getRowMap()->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: RangeMapExtractor incompatible to row map of input matrix.")

    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap.getMaxAllGlobalIndex() != input.getColMap()->getMaxAllGlobalIndex(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap.getGlobalNumElements() != input.getColMap()->getGlobalNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")
    TEUCHOS_TEST_FOR_EXCEPTION(fullDomainMap.getNodeNumElements() != input.getColMap()->getNodeNumElements(), Xpetra::Exceptions::Incompatible, "Xpetra::MatrixUtils::Split: DomainMapExtractor incompatible to domain map of input matrix.")*/

    std::vector<Teuchos::RCP<CrsMatrix> > subMatrices(numRows*numCols, Teuchos::null);

    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        // create empty CrsMatrix objects (we're reserving a little bit too much memory, but should be still reasonable)
        subMatrices[r*numCols+c] = Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build (rangeMapExtractor->getMap(r),input.getNodeMaxNumRowEntries());
      }
    }

#if 1

    // 1) loop over all rows of input matrix
    for (size_t rr = 0; rr < input.getRowMap()->getNodeNumElements(); rr++) {

      GlobalOrdinal growid = input.getRowMap()->getGlobalElement(rr); // LID -> GID (column)

      // Find the block id in range map extractor that belongs to same global row id.
      // We assume the global ids to be unique in the map extractor
      // TODO what about thyra mode? There we generate our global GIDs but they are ordered
      // by the way how the submaps are ordered.
      size_t rowBlockId = rangeMapExtractor->getMapIndexForGID(growid);

      // extract matrix entries from input matrix
      // we use global ids since we have to determine the corresponding
      // block column ids using the global ids anyway
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      input.getLocalRowView(rr, indices, vals);

      std::vector<Teuchos::Array<GlobalOrdinal> > blockColIdx (numCols, Teuchos::Array<GlobalOrdinal>());
      std::vector<Teuchos::Array<Scalar> >        blockColVals(numCols, Teuchos::Array<Scalar>());

      for(size_t i=0; i<(size_t)indices.size(); i++) {
        GlobalOrdinal gcolid = input.getColMap()->getGlobalElement(indices[i]);
        size_t colBlockId = domainMapExtractor->getMapIndexForGID(gcolid);
        blockColIdx [colBlockId].push_back(gcolid);
        blockColVals[colBlockId].push_back( vals[i]);
      }

      for (size_t c = 0; c < numCols; c++) {
        subMatrices[rowBlockId*numCols+c]->insertGlobalValues(growid,blockColIdx[c].view(0,blockColIdx[c].size()),blockColVals[c].view(0,blockColVals[c].size()));
      }

    }

    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        subMatrices[r*numCols+c]->fillComplete(domainMapExtractor->getMap(c), rangeMapExtractor->getMap(r));
      }
    }

    Teuchos::RCP<BlockedCrsMatrix> bA = Teuchos::rcp(new BlockedCrsMatrix(
        rangeMapExtractor, domainMapExtractor, 10 /*input.getRowMap()->getNodeMaxNumRowEntries()*/));


    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        bA->setMatrix(r,c,subMatrices[r*numCols+c]);
      }
    }

    return bA;

#else
    // loop over all block rows
    for (size_t r = 0; r < numRows; r++) {
      Teuchos::RCP<const Map> blockRowMap = rangeMapExtractor->getMap(r);


      // loop over all rows within block row
      for (size_t rr = 0; rr < blockRowMap->getNodeNumElements(); rr++) {
        std::vector<Teuchos::Array<GlobalOrdinal> > blockColIdx (numCols, Teuchos::Array<GlobalOrdinal>());
        std::vector<Teuchos::Array<Scalar> >        blockColVals(numCols, Teuchos::Array<Scalar>());

        // we need global ids since we might jump over rows in the input matrix!
        GlobalOrdinal growid = input.getRowMap()->getGlobalElement(rr); // LID -> GID (column)

        size_t nnz = input.getNodeMaxNumRowEntries();

        Teuchos::ArrayView<const GlobalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        input.getGlobalRowView(growid, indices, vals);

        //TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Xpetra::Exceptions::RuntimeError, "Xpetra::MatrixUtils::Split: number of nonzeros not equal to number of indices? Error.");

        // loop through all entries in current row
        for(size_t i=0; i<(size_t)indices.size(); i++) {
          GlobalOrdinal gcolid = indices[i];
          for (size_t c = 0; c < numCols; c++) {
            Teuchos::RCP<const Map> blockDomainMap = domainMapExtractor->getMap(c);
            if(blockDomainMap->isNodeGlobalElement(gcolid)) {
              Scalar gcolval = vals[i];
              blockColIdx[c].push_back(gcolid);
              blockColVals[c].push_back(gcolval);
            }
          }
        }

        for (size_t c = 0; c < numCols; c++) {
          std::cout << "BLOCK " << c << std::endl;
          for(size_t s = 0; s < blockColIdx[c].size(); s++) {
            Teuchos::Array<GlobalOrdinal> bci = blockColIdx[c];
            std::cout << "PROC " << input.getRowMap()->getComm()->getRank() << " colid " << bci[s] << " -> " << c << std::endl;
          }
          subMatrices[r*numCols+c]->insertGlobalValues(growid,blockColIdx[c].view(0,blockColIdx[c].size()),blockColVals[c].view(0,blockColVals[c].size()));
          std::cout << "added matrix entries row " << rr << " block column " << c << std::endl;

          //blockColIdx.clear();
          //blockColVals.clear();
        }
      }

      for (size_t c = 0; c < numCols; c++) {
        subMatrices[r*numCols+c]->fillComplete(domainMapExtractor->getMap(c), blockRowMap);
      }
    }

    Teuchos::RCP<BlockedCrsMatrix> bA = Teuchos::rcp(new BlockedCrsMatrix(
        rangeMapExtractor, domainMapExtractor, 10 /*input.getRowMap()->getNodeMaxNumRowEntries()*/));


    for (size_t r = 0; r < numRows; r++) {
      for (size_t c = 0; c < numCols; c++) {
        bA->setMatrix(r,c,subMatrices[r*numCols+c]);
      }
    }

    return bA;
#endif
  }

};

} // end namespace Xpetra

#define XPETRA_MATRIXUTILS_SHORT

#endif // PACKAGES_XPETRA_SUP_MATRIX_UTILS_HPP_
