// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MATRIXFACTORY_DECL_HPP
#define XPETRA_MATRIXFACTORY_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_MapExtractor_fwd.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_BlockedCrsMatrix_fwd.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_BlockedMap.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_BlockedVector.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class MatrixFactory {
#undef XPETRA_MATRIXFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

 private:
  //! Private constructor. This is a static class.
  MatrixFactory() {}

 public:
  /// Constructor for an empty, DynamicProfile matrix.
  /// Supports Epetra only, as DynamicProfile no longer exists in Tpetra.
  static RCP<Matrix> Build(const RCP<const Map>& rowMap);

  //! Constructor specifying the number of non-zeros for all rows.
  static RCP<Matrix> Build(const RCP<const Map>& rowMap, size_t maxNumEntriesPerRow);

  //! Constructor specifying the max number of non-zeros per row and providing column map
  static RCP<Matrix> Build(const RCP<const Map>& rowMap, const RCP<const Map>& colMap, size_t maxNumEntriesPerRow);

  //! Constructor specifying the (possibly different) number of entries per row and providing column map
  static RCP<Matrix> Build(const RCP<const Map>& rowMap, const RCP<const Map>& colMap, const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc);

  //! Constructor providing a local Kokkos::CrsMatrix together with a row and column map
  static RCP<Matrix> Build(
      const Teuchos::RCP<const Map>& rowMap,
      const Teuchos::RCP<const Map>& colMap,
      const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
      const Teuchos::RCP<Teuchos::ParameterList>& params = null);
  //! Constructor providing a local Kokkos::CrsMatrix together with all maps
  static RCP<Matrix> Build(
      const typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type& lclMatrix,
      const Teuchos::RCP<const Map>& rowMap,
      const Teuchos::RCP<const Map>& colMap,
      const Teuchos::RCP<const Map>& domainMap           = Teuchos::null,
      const Teuchos::RCP<const Map>& rangeMap            = Teuchos::null,
      const Teuchos::RCP<Teuchos::ParameterList>& params = null);

  //! Constructor specifying (possibly different) number of entries in each row.
  static RCP<Matrix> Build(const RCP<const Map>& rowMap, const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc);

  //! Constructor specifying graph
  static RCP<Matrix> Build(const RCP<const CrsGraph>& graph, const RCP<ParameterList>& paramList = Teuchos::null);

  //! Constructor specifying graph and values array
  static RCP<Matrix> Build(const RCP<const CrsGraph>& graph,
                           typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type::values_type& values,
                           const RCP<ParameterList>& paramList = Teuchos::null);

  //! Constructor for creating a diagonal Xpetra::Matrix using the entries of a given vector for the diagonal
  static RCP<Matrix> Build(const RCP<const Vector>& diagonal);

  //! Constructor to create a Matrix using a fusedImport-style construction.  The originalMatrix must be a Xpetra::CrsMatrixWrap under the hood or this will fail.
  static RCP<Matrix> Build(const RCP<const Matrix>& sourceMatrix, const Import& importer, const RCP<const Map>& domainMap = Teuchos::null, const RCP<const Map>& rangeMap = Teuchos::null, const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  //! Constructor to create a Matrix using a fusedExport-style construction.  The originalMatrix must be a Xpetra::CrsMatrixWrap under the hood or this will fail.
  static RCP<Matrix> Build(const RCP<const Matrix>& sourceMatrix, const Export& exporter, const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, const Teuchos::RCP<Teuchos::ParameterList>& params);

  //! Constructor to create a Matrix using a fusedImport-style construction.  The originalMatrix must be a Xpetra::CrsMatrixWrap under the hood or this will fail.
  static RCP<Matrix> Build(const RCP<const Matrix>& sourceMatrix, const Import& RowImporter, const Import& DomainImporter, const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, const Teuchos::RCP<Teuchos::ParameterList>& params);

  //! Constructor to create a Matrix using a fusedExport-style construction.  The originalMatrix must be a Xpetra::CrsMatrixWrap under the hood or this will fail.
  static RCP<Matrix> Build(const RCP<const Matrix>& sourceMatrix, const Export& RowExporter, const Export& DomainExporter, const RCP<const Map>& domainMap = Teuchos::null, const RCP<const Map>& rangeMap = Teuchos::null, const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  //! create an explicit copy of a given matrix
  //! This routine supports blocked and single-block operators
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildCopy(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A, bool setFixedBlockSize = true);
};
#define XPETRA_MATRIXFACTORY_SHORT

}  // namespace Xpetra

#define XPETRA_MATRIXFACTORY_SHORT
#endif  // ifndef XPETRA_MATRIXFACTORY_DECL_HPP
