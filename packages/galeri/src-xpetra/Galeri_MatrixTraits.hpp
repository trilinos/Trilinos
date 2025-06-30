// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/*
  Direct translation of parts of Galeri to use Tpetra or Xpetra rather than Epetra. Epetra also supported.
*/

#ifndef GALERI_XPETRAMATRIXTRAITS_HPP
#define GALERI_XPETRAMATRIXTRAITS_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>

#include "Galeri_config.h"
#ifdef HAVE_GALERI_XPETRA
// needed for the specialized traits:
#include "Xpetra_Map.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#endif

namespace Galeri {

namespace Xpetra {

/* Default traits */
/* These traits work for the following couples of (Map,Matrix):
   - Map = Tpetra::Map<...>,       and Matrix = Tpetra::CrsMatrix<...>
   - Map = Xpetra::TpetraMap<...> and Matrix = Xpetra::TpetraCrsMatrix<...>
   - Map = Xpetra::EpetraMap,     and Matrix = Xpetra::EpetraCrsMatrix
*/
template <class Map, class Matrix>
class MatrixTraits {
  using scalar_type         = typename Matrix::scalar_type;
  using local_ordinal_type  = typename Matrix::local_ordinal_type;
  using global_ordinal_type = typename Matrix::global_ordinal_type;
  using node_type           = typename Matrix::node_type;
  using local_matrix_type   = typename ::Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>::local_matrix_device_type;

 public:
  static Teuchos::RCP<Matrix> Build(const Teuchos::RCP<const Map> &rowMap, size_t maxNumEntriesPerRow) { return Teuchos::rcp(new Matrix(rowMap, maxNumEntriesPerRow)); };
  static Teuchos::RCP<Matrix> Build(local_matrix_type &lclMatrix, const Teuchos::RCP<const Map> &rowMap, const Teuchos::RCP<const Map> &colMap, const Teuchos::RCP<const Map> &domainMap, const Teuchos::RCP<const Map> &rangeMap) { return Teuchos::rcp(new Matrix(lclMatrix, rowMap, colMap, domainMap, rangeMap)); };
};

#ifdef HAVE_GALERI_XPETRA

/* Specialized traits for:
   - Map = Xpetra::Map<...>, Matrix = Xpetra::CrsMatrix<...> */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class MatrixTraits<::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>, ::Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> {
  using local_matrix_type = typename ::Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;

 public:
  static Teuchos::RCP<::Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap, size_t maxNumEntriesPerRow)
  // Use the CrsMatrixFactory to decide what kind of matrix to create (Xpetra::TpetraCrsMatrix or Xpetra::EpetraCrsMatrix).
  { return ::Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, maxNumEntriesPerRow); };
  static Teuchos::RCP<::Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(local_matrix_type &lclMatrix, const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap, const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &colMap, const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap, const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap) {
    auto crs = ::Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(lclMatrix, rowMap, colMap, domainMap, rangeMap);
    return crs;
  }
};

/* Specialized traits for:
   - Map = Xpetra::Map<...>, Matrix = Xpetra::Matrix<...> */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class MatrixTraits<::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>, ::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> {
  using local_matrix_type = typename ::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;

 public:
  static Teuchos::RCP<::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap, size_t maxNumEntriesPerRow)
  // Use the CrsMatrixFactory to decide what kind of matrix to create (Xpetra::TpetraCrsMatrix or Xpetra::EpetraCrsMatrix).
  { return ::Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap, maxNumEntriesPerRow); };
  static Teuchos::RCP<::Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Build(local_matrix_type &lclMatrix, const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &rowMap, const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &colMap, const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &domainMap, const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &rangeMap) { return ::Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(lclMatrix, rowMap, colMap, domainMap, rangeMap); }
};

#endif

}  // namespace Xpetra
}  // namespace Galeri

#endif  // ifndef GALERI_XPETRAMATRIXTRAITS_HPP
