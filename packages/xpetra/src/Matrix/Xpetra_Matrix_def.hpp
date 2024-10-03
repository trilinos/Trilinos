// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MATRIX_DEF_HPP
#define XPETRA_MATRIX_DEF_HPP

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"
#include "Xpetra_Matrix_decl.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"
#include "Xpetra_MatrixView.hpp"
#include "Xpetra_Operator.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Matrix() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Matrix() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateView(viewLabel_t viewLabel, const RCP<const Map> &rowMap, const RCP<const Map> &colMap) {
  TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == true, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.CreateView(): a view labeled '" + viewLabel + "' already exist.");
  RCP<MatrixView> view = rcp(new MatrixView(rowMap, colMap));
  operatorViewTable_.put(viewLabel, view);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateView(const viewLabel_t viewLabel, const RCP<const Matrix> &A, bool transposeA, const RCP<const Matrix> &B, bool transposeB) {
  RCP<const Map> domainMap = Teuchos::null;
  RCP<const Map> rangeMap  = Teuchos::null;

  const size_t blkSize = 1;
  std::vector<size_t> stridingInfo(1, blkSize);
  LocalOrdinal stridedBlockId = -1;

  if (A->IsView(viewLabel)) {
    rangeMap  = transposeA ? A->getColMap(viewLabel) : A->getRowMap(viewLabel);
    domainMap = transposeA ? A->getRowMap(viewLabel) : A->getColMap(viewLabel);  // will be overwritten if B != Teuchos::null

  } else {
    rangeMap  = transposeA ? A->getDomainMap() : A->getRangeMap();
    domainMap = transposeA ? A->getRangeMap() : A->getDomainMap();

    if (viewLabel == "stridedMaps") {
      rangeMap  = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(rangeMap, stridingInfo, stridedBlockId);
      domainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(domainMap, stridingInfo, stridedBlockId);
    }
  }

  if (B != Teuchos::null) {
    // B has strided Maps

    if (B->IsView(viewLabel)) {
      domainMap = transposeB ? B->getRowMap(viewLabel) : B->getColMap(viewLabel);

    } else {
      domainMap = transposeB ? B->getRangeMap() : B->getDomainMap();

      if (viewLabel == "stridedMaps")
        domainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(domainMap, stridingInfo, stridedBlockId);
    }
  }

  if (IsView(viewLabel))
    RemoveView(viewLabel);

  CreateView(viewLabel, rangeMap, domainMap);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PrintViews(Teuchos::FancyOStream &out) const {
  int last = out.getOutputToRootOnly();
  Teuchos::OSTab tab(out);
  out.setOutputToRootOnly(0);
  Teuchos::Array<viewLabel_t> viewLabels;
  Teuchos::Array<RCP<MatrixView>> viewList;
  operatorViewTable_.arrayify(viewLabels, viewList);
  out << "views associated with this operator" << std::endl;
  for (int i = 0; i < viewLabels.size(); ++i)
    out << viewLabels[i] << std::endl;
  out.setOutputToRootOnly(last);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RemoveView(const viewLabel_t viewLabel) {
  TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.RemoveView(): view '" + viewLabel + "' does not exist.");
  TEUCHOS_TEST_FOR_EXCEPTION(viewLabel == GetDefaultViewLabel(), Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.RemoveView(): view '" + viewLabel + "' is the default view and cannot be removed.");
  operatorViewTable_.remove(viewLabel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const std::string Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SwitchToView(const viewLabel_t viewLabel) {
  TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.SwitchToView(): view '" + viewLabel + "' does not exist.");
  viewLabel_t oldViewLabel = GetCurrentViewLabel();
  currentViewLabel_        = viewLabel;
  return oldViewLabel;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::IsView(const viewLabel_t viewLabel) const {
  return operatorViewTable_.containsKey(viewLabel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const std::string Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SwitchToDefaultView() { return SwitchToView(GetDefaultViewLabel()); }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const std::string &Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetDefaultViewLabel() const { return defaultViewLabel_; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const std::string &Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetCurrentViewLabel() const { return currentViewLabel_; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRowMap() const { return getRowMap(GetCurrentViewLabel()); }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getRowMap(viewLabel_t viewLabel) const {
  TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.GetRowMap(): view '" + viewLabel + "' does not exist.");
  return operatorViewTable_.get(viewLabel)->GetRowMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getColMap() const { return getColMap(GetCurrentViewLabel()); }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> &Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getColMap(viewLabel_t viewLabel) const {
  TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Matrix.GetColMap(): view '" + viewLabel + "' does not exist.");
  return operatorViewTable_.get(viewLabel)->GetColMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetFixedBlockSize(LocalOrdinal blksize, GlobalOrdinal offset) {
  TEUCHOS_TEST_FOR_EXCEPTION(isFillComplete() == false, Exceptions::RuntimeError, "Xpetra::Matrix::SetFixedBlockSize(): operator is not filled and completed.");  // TODO: do we need this? we just wanna "copy" the domain and range map
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(Teuchos::as<size_t>(blksize));
  LocalOrdinal stridedBlockId = -1;

  RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>> stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      this->getRangeMap(),
      stridingInfo,
      stridedBlockId,
      offset);
  RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
      this->getDomainMap(),
      stridingInfo,
      stridedBlockId,
      offset);

  if (IsFixedBlockSizeSet()) RemoveView("stridedMaps");
  CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetFixedBlockSize() const {
  if (IsFixedBlockSizeSet()) {
    Teuchos::RCP<const StridedMap<LocalOrdinal, GlobalOrdinal, Node>> rangeMap  = Teuchos::rcp_dynamic_cast<const StridedMap<LocalOrdinal, GlobalOrdinal, Node>>(getRowMap("stridedMaps"));
    Teuchos::RCP<const StridedMap<LocalOrdinal, GlobalOrdinal, Node>> domainMap = Teuchos::rcp_dynamic_cast<const StridedMap<LocalOrdinal, GlobalOrdinal, Node>>(getColMap("stridedMaps"));
    TEUCHOS_TEST_FOR_EXCEPTION(rangeMap == Teuchos::null, Exceptions::BadCast, "Xpetra::Matrix::GetFixedBlockSize(): rangeMap is not of type StridedMap");
    TEUCHOS_TEST_FOR_EXCEPTION(domainMap == Teuchos::null, Exceptions::BadCast, "Xpetra::Matrix::GetFixedBlockSize(): domainMap is not of type StridedMap");
    TEUCHOS_TEST_FOR_EXCEPTION(domainMap->getFixedBlockSize() != rangeMap->getFixedBlockSize(), Exceptions::RuntimeError, "Xpetra::Matrix::GetFixedBlockSize(): block size of rangeMap and domainMap are different.");
    return Teuchos::as<LocalOrdinal>(domainMap->getFixedBlockSize());  // TODO: why LocalOrdinal?
  } else
    // TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "Xpetra::Matrix::GetFixedBlockSize(): no strided maps available."); // TODO remove this
    return 1;
}  // TODO: why LocalOrdinal?

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::IsFixedBlockSizeSet() const {
  return IsView("stridedMaps");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetMaxEigenvalueEstimate(Scalar const &sigma) {
  operatorViewTable_.get(GetCurrentViewLabel())->SetMaxEigenvalueEstimate(sigma);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetMaxEigenvalueEstimate() const {
  return operatorViewTable_.get(GetCurrentViewLabel())->GetMaxEigenvalueEstimate();
}

}  // namespace Xpetra

#define XPETRA_MATRIX_SHORT
#endif  // XPETRA_MATRIX_DECL_HPP
