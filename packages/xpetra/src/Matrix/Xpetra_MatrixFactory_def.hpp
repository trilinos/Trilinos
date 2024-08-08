// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MATRIXFACTORY_DEF_HPP
#define XPETRA_MATRIXFACTORY_DEF_HPP

#include "Xpetra_MatrixFactory2_decl.hpp"
#include "Xpetra_MatrixFactory_decl.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Map>& rowMap) {
  return rcp(new CrsMatrixWrap(rowMap));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Map>& rowMap, size_t maxNumEntriesPerRow) {
  return rcp(new CrsMatrixWrap(rowMap, maxNumEntriesPerRow));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Map>& rowMap, const RCP<const Map>& colMap, size_t maxNumEntriesPerRow) {
  return rcp(new CrsMatrixWrap(rowMap, colMap, maxNumEntriesPerRow));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Map>& rowMap, const RCP<const Map>& colMap, const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc) {
  return rcp(new CrsMatrixWrap(rowMap, colMap, NumEntriesPerRowToAlloc));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
    const Teuchos::RCP<const Map>& rowMap,
    const Teuchos::RCP<const Map>& colMap,
    const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
    const Teuchos::RCP<Teuchos::ParameterList>& params) {
  XPETRA_MONITOR("MatrixFactory::Build");
  return rcp(new CrsMatrixWrap(rowMap, colMap, lclMatrix, params));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
    const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
    const Teuchos::RCP<const Map>& rowMap,
    const Teuchos::RCP<const Map>& colMap,
    const Teuchos::RCP<const Map>& domainMap,
    const Teuchos::RCP<const Map>& rangeMap,
    const Teuchos::RCP<Teuchos::ParameterList>& params) {
  XPETRA_MONITOR("MatrixFactory::Build");
  return rcp(new CrsMatrixWrap(lclMatrix, rowMap, colMap, domainMap, rangeMap, params));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Map>& rowMap, const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc) {
  return rcp(new CrsMatrixWrap(rowMap, NumEntriesPerRowToAlloc));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const CrsGraph>& graph, const RCP<ParameterList>& paramList) {
  return rcp(new CrsMatrixWrap(graph, paramList));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const CrsGraph>& graph,
                                                                                                                               typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::values_type& values,
                                                                                                                               const RCP<ParameterList>& paramList) {
  return rcp(new CrsMatrixWrap(graph, values, paramList));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Vector>& diagonal) {
  RCP<const BlockedVector> bdiagonal = Teuchos::rcp_dynamic_cast<const BlockedVector>(diagonal);
  Teuchos::RCP<Matrix> mtx           = Teuchos::null;

  if (bdiagonal == Teuchos::null) {
    Teuchos::ArrayRCP<const Scalar> vals                     = diagonal->getData(0);
    LocalOrdinal numMyElements                               = diagonal->getMap()->getLocalNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = diagonal->getMap()->getLocalElementList();

    mtx = Teuchos::rcp(new CrsMatrixWrap(diagonal->getMap(), 1));

    for (LocalOrdinal i = 0; i < numMyElements; ++i) {
      mtx->insertGlobalValues(myGlobalElements[i],
                              Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i]),
                              Teuchos::tuple<Scalar>(vals[i]));
    }
    mtx->fillComplete();
  } else {
    RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(bdiagonal->getBlockedMap(), bdiagonal->getBlockedMap(), 1));

    for (size_t r = 0; r < bdiagonal->getBlockedMap()->getNumMaps(); ++r) {
      if (!bdiagonal->getMultiVector(r).is_null()) {
        const RCP<MultiVector> subvec = bdiagonal->getMultiVector(r);
        bop->setMatrix(r, r, Build(subvec->getVector(0)));
      }
    }
    bop->fillComplete();
    mtx = BuildCopy(bop);
  }

  return mtx;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Matrix>& sourceMatrix, const Import& importer, const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, const Teuchos::RCP<Teuchos::ParameterList>& params) {
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(sourceMatrix);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<CrsMatrix> originalCrs = crsOp->getCrsMatrix();
  RCP<CrsMatrix> newCrs      = CrsMatrixFactory::Build(originalCrs, importer, domainMap, rangeMap, params);
  if (newCrs->hasMatrix())
    return rcp(new CrsMatrixWrap(newCrs));
  else
    return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Matrix>& sourceMatrix, const Export& exporter, const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, const Teuchos::RCP<Teuchos::ParameterList>& params) {
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(sourceMatrix);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<CrsMatrix> originalCrs = crsOp->getCrsMatrix();
  return rcp(new CrsMatrixWrap(CrsMatrixFactory::Build(originalCrs, exporter, domainMap, rangeMap, params)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Matrix>& sourceMatrix, const Import& RowImporter, const Import& DomainImporter, const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, const Teuchos::RCP<Teuchos::ParameterList>& params) {
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(sourceMatrix);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<CrsMatrix> originalCrs = crsOp->getCrsMatrix();
  RCP<CrsMatrix> newCrs      = CrsMatrixFactory::Build(originalCrs, RowImporter, Teuchos::rcpFromRef(DomainImporter), domainMap, rangeMap, params);
  if (newCrs->hasMatrix())
    return rcp(new CrsMatrixWrap(newCrs));
  else
    return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& sourceMatrix, const Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node>& RowExporter, const Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node>& DomainExporter, const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap, const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap, const Teuchos::RCP<Teuchos::ParameterList>& params) {
  RCP<const CrsMatrixWrap> crsOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(sourceMatrix);
  if (crsOp == Teuchos::null)
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

  RCP<CrsMatrix> originalCrs = crsOp->getCrsMatrix();
  RCP<CrsMatrix> newCrs      = CrsMatrixFactory::Build(originalCrs, RowExporter, Teuchos::rcpFromRef(DomainExporter), domainMap, rangeMap, params);
  if (newCrs->hasMatrix())
    return rcp(new CrsMatrixWrap(newCrs));
  else
    return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A, bool setFixedBlockSize) {
  RCP<const BlockedCrsMatrix> input = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(A);
  if (input == Teuchos::null)
    return Xpetra::MatrixFactory2<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(A, setFixedBlockSize);

  // deep copy of MapExtractors (and underlying maps)
  RCP<const MapExtractor> rgMapExt = Teuchos::rcp(new MapExtractor(*(input->getRangeMapExtractor())));
  RCP<const MapExtractor> doMapExt = Teuchos::rcp(new MapExtractor(*(input->getDomainMapExtractor())));

  // create new BlockedCrsMatrix object
  RCP<BlockedCrsMatrix> bop = Teuchos::rcp(new BlockedCrsMatrix(rgMapExt, doMapExt, input->getLocalMaxNumRowEntries()));

  for (size_t r = 0; r < input->Rows(); ++r) {
    for (size_t c = 0; c < input->Cols(); ++c)
      if (input->getMatrix(r, c) != Teuchos::null) {
        // make a deep copy of the matrix
        // This is a recursive call to this function
        RCP<Matrix> mat =
            Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildCopy(input->getMatrix(r, c), setFixedBlockSize);
        bop->setMatrix(r, c, mat);
      }
  }

  if (input->isFillComplete())
    bop->fillComplete();
  return bop;
}

}  // namespace Xpetra

#endif
