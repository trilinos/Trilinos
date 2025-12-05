// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_THYRAUTILS_HPP
#define XPETRA_THYRAUTILS_HPP

#include "Xpetra_ConfigDefs.hpp"
#ifdef HAVE_XPETRA_THYRA

#include <typeinfo>

#include "Tpetra_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_BlockedMap.hpp"
#include "Xpetra_BlockedMultiVector.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"
#include "Xpetra_MapUtils.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_SpmdVectorSpaceBase.hpp>
#include <Thyra_ProductVectorSpaceBase.hpp>
#include <Thyra_ProductMultiVectorBase.hpp>
#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_DefaultProductVectorSpace.hpp>
#include <Thyra_DefaultBlockedLinearOp.hpp>
#include <Thyra_LinearOpBase.hpp>
#include "Thyra_DiagonalLinearOpBase.hpp"
#include <Thyra_DetachedMultiVectorView.hpp>
#include <Thyra_MultiVectorStdOps.hpp>

#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraVectorSpace.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Xpetra_TpetraMap.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>

namespace Xpetra {

template <class Scalar,
          class LocalOrdinal  = int,
          class GlobalOrdinal = LocalOrdinal,
          class Node          = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class ThyraUtils {
 private:
#undef XPETRA_THYRAUTILS_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  static Teuchos::RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int>>& comm, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId = -1, GlobalOrdinal offset = 0);

  static Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  // const version
  static Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> v, const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  // non-const version
  static Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> v, const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  static bool isTpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op);

  static bool isBlockedOperator(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op);

  static Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op);

  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<Thyra::LinearOpBase<Scalar>>& op);

  static Teuchos::RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetraOperator(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op);

  static Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetraOperator(const Teuchos::RCP<Thyra::LinearOpBase<Scalar>>& op);

  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<Thyra::DiagonalLinearOpBase<Scalar>>& op);

  static Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::DiagonalLinearOpBase<Scalar>>& op);

  static Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
  toThyra(Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map);

  static Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
  toThyraMultiVector(Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec);

  static Teuchos::RCP<const Thyra::VectorBase<Scalar>>
  toThyraVector(Teuchos::RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec);

  // update Thyra multi vector with data from Xpetra multi vector
  // In case of a Thyra::ProductMultiVector the Xpetra::MultiVector is splitted into its subparts using a provided MapExtractor
  static void
  updateThyra(Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> source, Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mapExtractor, const Teuchos::RCP<Thyra::MultiVectorBase<Scalar>>& target);

  static Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>
  toThyra(const Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat);

  static Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
  toThyra(const Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat);

  static Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
  toThyra(const Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat);

};  // end Utils class

}  // end namespace Xpetra

#define XPETRA_THYRAUTILS_SHORT
#endif  // HAVE_XPETRA_THYRA

#endif  // XPETRA_THYRAUTILS_HPP
