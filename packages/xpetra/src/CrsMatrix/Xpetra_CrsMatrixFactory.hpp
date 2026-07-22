// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_CRSMATRIXFACTORY_HPP
#define XPETRA_CRSMATRIXFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_CrsMatrix.hpp"

#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Xpetra_TpetraBlockCrsMatrix.hpp"

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class CrsMatrixFactory {
 private:
  //! Private constructor. This is a static class.
  CrsMatrixFactory() {}

 public:
  //! Constructor for empty matrix (intended use is an import/export target - can't insert entries directly)
  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap) {
    if (rowMap->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, 0));

    XPETRA_FACTORY_END;
  }

  //! Constructor specifying fixed number of entries for each row.
  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
        size_t maxNumEntriesPerRow,
        const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (rowMap->lib() == UseTpetra)
      return Teuchos::rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, maxNumEntriesPerRow, plist));

    XPETRA_FACTORY_END;
  }

  //! Constructor specifying (possibly different) number of entries in each row.
  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
        const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
        const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null) {
    if (rowMap->lib() == UseTpetra)
      return Teuchos::rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, NumEntriesPerRowToAlloc, plist));

    XPETRA_FACTORY_END;
  }

  //! Constructor specifying column Map and fixed number of entries for each row.
  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
        const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& colMap,
        size_t maxNumEntriesPerRow,
        const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (rowMap->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, maxNumEntriesPerRow, plist));

    XPETRA_FACTORY_END;
  }

  //! Constructor specifying column Map and number of entries in each row.
  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap, const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& colMap, const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc, const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (rowMap->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, NumEntriesPerRowToAlloc, plist));

    XPETRA_FACTORY_END;
  }

  //! Constructor specifying a previously constructed graph.
  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(const Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>& graph, const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (graph->getRowMap()->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(graph, plist));

    XPETRA_FACTORY_END;
  }

  //! Constructor specifying a previously constructed graph and values array
  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(const Teuchos::RCP<const CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>& graph, typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::values_type& values, const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (graph->getRowMap()->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(graph, values, plist));

    XPETRA_FACTORY_END;
  }

  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(
      const Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& sourceMatrix,
      const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
      const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap = Teuchos::null,
      const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap  = Teuchos::null,
      const Teuchos::RCP<Teuchos::ParameterList>& params                 = Teuchos::null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (sourceMatrix->getRowMap()->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix, importer, domainMap, rangeMap, params));

    XPETRA_FACTORY_END;
  }

  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(
      const Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& sourceMatrix,
      const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
      const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap = Teuchos::null,
      const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap  = Teuchos::null,
      const Teuchos::RCP<Teuchos::ParameterList>& params                 = Teuchos::null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (sourceMatrix->getRowMap()->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix, exporter, domainMap, rangeMap, params));

    XPETRA_FACTORY_END;
  }

  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(
      const Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& sourceMatrix,
      const Import<LocalOrdinal, GlobalOrdinal, Node>& RowImporter,
      const RCP<const Import<LocalOrdinal, GlobalOrdinal, Node>> DomainImporter,
      const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap,
      const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap,
      const Teuchos::RCP<Teuchos::ParameterList>& params) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (sourceMatrix->getRowMap()->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix, RowImporter, DomainImporter, domainMap, rangeMap, params));

    XPETRA_FACTORY_END;
  }

  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(
      const Teuchos::RCP<const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& sourceMatrix,
      const Export<LocalOrdinal, GlobalOrdinal, Node>& RowExporter,
      const RCP<const Export<LocalOrdinal, GlobalOrdinal, Node>> DomainExporter,
      const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap,
      const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap,
      const Teuchos::RCP<Teuchos::ParameterList>& params) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (sourceMatrix->getRowMap()->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(sourceMatrix, RowExporter, DomainExporter, domainMap, rangeMap, params));

    XPETRA_FACTORY_END;
  }

  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& colMap,
      const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
      const Teuchos::RCP<Teuchos::ParameterList>& params = null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (rowMap->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, lclMatrix, params));

    XPETRA_FACTORY_END;
  }

  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(
      const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& colMap,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap = Teuchos::null,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap  = Teuchos::null,
      const Teuchos::RCP<Teuchos::ParameterList>& params                          = null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (rowMap->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, rowMap, colMap, domainMap, rangeMap, params));

    XPETRA_FACTORY_END;
  }

  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(
      const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& colMap,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap,
      const Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node>>& importer,
      const Teuchos::RCP<const Export<LocalOrdinal, GlobalOrdinal, Node>>& exporter,
      const Teuchos::RCP<Teuchos::ParameterList>& params = null) {
    XPETRA_MONITOR("CrsMatrixFactory::Build");

    if (rowMap->lib() == UseTpetra)
      return rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lclMatrix, rowMap, colMap, domainMap, rangeMap, importer, exporter, params));

    XPETRA_FACTORY_END;
  }

  // Builds a BlockCrsMatrix
  static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> BuildBlock(
      const Teuchos::RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>& blockGraph,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap,
      const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap,
      LocalOrdinal blockSize) {
    XPETRA_MONITOR("CrsMatrixFactory::BuildBlock");

    if (domainMap->lib() == UseTpetra) {
      return rcp(new Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(blockGraph, domainMap, rangeMap, blockSize));
    }
    XPETRA_FACTORY_END;
  }
};

}  // namespace Xpetra

#define XPETRA_CRSMATRIXFACTORY_SHORT
#endif
