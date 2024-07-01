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
//
// ***********************************************************************
//
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MatrixFactory2 {
#undef XPETRA_MATRIXFACTORY2_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildCopy(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A, bool setFixedBlockSize = true);
};
#define XPETRA_MATRIXFACTORY2_SHORT

// template<>
// class MatrixFactory2<double,int,int,typename Xpetra::Matrix<double, int, int>::node_type> {
template <class Node>
class MatrixFactory2<double, int, int, Node> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  // typedef Matrix<double, int, GlobalOrdinal>::node_type Node;
#undef XPETRA_MATRIXFACTORY2_SHORT
#include "Xpetra_UseShortNames.hpp"
 public:
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildCopy(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A, bool setFixedBlockSize = true) {
    RCP<const CrsMatrixWrap> oldOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A);
    if (oldOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> oldCrsOp = oldOp->getCrsMatrix();

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
    RCP<const EpetraCrsMatrixT<GlobalOrdinal, Node> > oldECrsOp = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GlobalOrdinal, Node> >(oldCrsOp);
    if (oldECrsOp != Teuchos::null) {
      // Underlying matrix is Epetra
      RCP<CrsMatrix> newECrsOp(new EpetraCrsMatrixT<GlobalOrdinal, Node>(*oldECrsOp));
      RCP<CrsMatrixWrap> newOp(new CrsMatrixWrap(newECrsOp));
      if (setFixedBlockSize)
        newOp->SetFixedBlockSize(A->GetFixedBlockSize());
      return newOp;
    }
#endif
#endif

#ifdef HAVE_XPETRA_TPETRA
    // Underlying matrix is Tpetra
    RCP<const TpetraCrsMatrix> oldTCrsOp = Teuchos::rcp_dynamic_cast<const TpetraCrsMatrix>(oldCrsOp);
    if (oldTCrsOp != Teuchos::null) {
      RCP<CrsMatrix> newTCrsOp(new TpetraCrsMatrix(*oldTCrsOp));
      RCP<CrsMatrixWrap> newOp(new CrsMatrixWrap(newTCrsOp));
      if (setFixedBlockSize)
        newOp->SetFixedBlockSize(A->GetFixedBlockSize());
      return newOp;
    }
    return Teuchos::null;
#else
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::EpetraCrsMatrix or Xpetra::TpetraCrsMatrix failed");
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);  // make compiler happy
#endif

  }  // BuildCopy
};

#define XPETRA_MATRIXFACTORY2_SHORT

#ifdef HAVE_XPETRA_INT_LONG_LONG
// template<>
// class MatrixFactory2<double,int,long long,typename Xpetra::Matrix<double, int, long long>::node_type> {
template <class Node>
class MatrixFactory2<double, int, long long, Node> {
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  // typedef Matrix<double, int, GlobalOrdinal>::node_type Node;
#undef XPETRA_MATRIXFACTORY2_SHORT
#include "Xpetra_UseShortNames.hpp"
 public:
  static RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildCopy(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A, bool setFixedBlockSize = true) {
    RCP<const CrsMatrixWrap> oldOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A);
    if (oldOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> oldCrsOp = oldOp->getCrsMatrix();

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
    RCP<const EpetraCrsMatrixT<GlobalOrdinal, Node> > oldECrsOp = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GlobalOrdinal, Node> >(oldCrsOp);
    if (oldECrsOp != Teuchos::null) {
      // Underlying matrix is Epetra
      RCP<CrsMatrix> newECrsOp(new EpetraCrsMatrixT<GlobalOrdinal, Node>(*oldECrsOp));
      RCP<CrsMatrixWrap> newOp(new CrsMatrixWrap(newECrsOp));
      if (setFixedBlockSize)
        newOp->SetFixedBlockSize(A->GetFixedBlockSize());
      return newOp;
    }
#endif
#endif

#ifdef HAVE_XPETRA_TPETRA
    // Underlying matrix is Tpetra
    RCP<const TpetraCrsMatrix> oldTCrsOp = Teuchos::rcp_dynamic_cast<const TpetraCrsMatrix>(oldCrsOp);
    if (oldTCrsOp != Teuchos::null) {
      RCP<CrsMatrix> newTCrsOp(new TpetraCrsMatrix(*oldTCrsOp));
      RCP<CrsMatrixWrap> newOp(new CrsMatrixWrap(newTCrsOp));
      if (setFixedBlockSize)
        newOp->SetFixedBlockSize(A->GetFixedBlockSize());
      return newOp;
    }
#else
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::EpetraCrsMatrix or Xpetra::TpetraCrsMatrix failed");
#endif

    return Teuchos::null;  // make compiler happy
  }
};
#endif  // HAVE_XPETRA_INT_LONG_LONG

#define XPETRA_MATRIXFACTORY2_SHORT

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
#define XPETRA_MATRIXFACTORY2_SHORT
#endif
