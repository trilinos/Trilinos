// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef HAVE_XPETRA_THYRAUTILS_DEF_HPP
#define HAVE_XPETRA_THYRAUTILS_DEF_HPP

#ifdef HAVE_XPETRA_THYRA

#include "Xpetra_BlockedCrsMatrix.hpp"

#include "Xpetra_ThyraUtils_decl.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int>>& comm, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId, GlobalOrdinal offset) {
  Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map = toXpetra(vectorSpace, comm);

  if (stridedBlockId == -1) {
    TEUCHOS_TEST_FOR_EXCEPT(map->getLocalNumElements() % stridingInfo.size() != 0);
  } else {
    TEUCHOS_TEST_FOR_EXCEPT(map->getLocalNumElements() % stridingInfo[stridedBlockId] != 0);
  }

  Teuchos::RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>> ret = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, stridingInfo, stridedBlockId, offset);
  return ret;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using ThyVecSpaceBase     = Thyra::VectorSpaceBase<Scalar>;
  using ThyProdVecSpaceBase = Thyra::ProductVectorSpaceBase<Scalar>;
  using ThyUtils            = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  RCP<Map> resultMap                             = Teuchos::null;
  RCP<const ThyProdVecSpaceBase> prodVectorSpace = rcp_dynamic_cast<const ThyProdVecSpaceBase>(vectorSpace);
  if (prodVectorSpace != Teuchos::null) {
    // SPECIAL CASE: product Vector space
    // collect all submaps to store them in a hierarchical BlockedMap object
    TEUCHOS_TEST_FOR_EXCEPTION(prodVectorSpace->numBlocks() == 0, std::logic_error, "Found a product vector space with zero blocks.");
    std::vector<RCP<const Map>> mapsThyra(prodVectorSpace->numBlocks(), Teuchos::null);
    std::vector<RCP<const Map>> mapsXpetra(prodVectorSpace->numBlocks(), Teuchos::null);
    for (int b = 0; b < prodVectorSpace->numBlocks(); ++b) {
      RCP<const ThyVecSpaceBase> bv = prodVectorSpace->getBlock(b);
      // can be of type Map or BlockedMap (containing Thyra GIDs)
      mapsThyra[b] = ThyUtils::toXpetra(bv, comm);  // recursive call
    }

    // get offsets for submap GIDs
    // we need that for the full map (Xpetra GIDs)
    std::vector<GlobalOrdinal> gidOffsets(prodVectorSpace->numBlocks(), 0);
    for (int i = 1; i < prodVectorSpace->numBlocks(); ++i) {
      gidOffsets[i] = mapsThyra[i - 1]->getMaxAllGlobalIndex() + gidOffsets[i - 1] + 1;
    }

    for (int b = 0; b < prodVectorSpace->numBlocks(); ++b) {
      RCP<const ThyVecSpaceBase> bv = prodVectorSpace->getBlock(b);
      // map can be of type Map or BlockedMap (containing Xpetra style GIDs)
      mapsXpetra[b] = MapUtils::transformThyra2XpetraGIDs(*mapsThyra[b], gidOffsets[b]);
    }

    resultMap = Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(mapsXpetra, mapsThyra));
  } else {
#ifdef HAVE_XPETRA_TPETRA
    // STANDARD CASE: no product map
    // check whether we have a Tpetra based Thyra operator
    Teuchos::RCP<const Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetra_vsc = Teuchos::rcp_dynamic_cast<const Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(vectorSpace);
    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(tpetra_vsc) == true, Xpetra::Exceptions::RuntimeError, "toXpetra failed to convert provided vector space to Thyra::TpetraVectorSpace. This is the general implementation for Tpetra only.");

    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> TpMap;
    typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpVector;
    typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
    typedef Thyra::VectorBase<Scalar> ThyVecBase;
    RCP<ThyVecBase> rgVec = Thyra::createMember<Scalar>(vectorSpace, std::string("label"));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgVec));
    RCP<const TpVector> rgTpetraVec = TOE::getTpetraVector(rgVec);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgTpetraVec));
    RCP<const TpMap> rgTpetraMap = rgTpetraVec->getMap();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgTpetraMap));
    resultMap = Xpetra::toXpetraNonConst(rgTpetraMap);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Cannot transform Thyra::VectorSpace to Xpetra::Map. This is the general implementation for Tpetra only, but Tpetra is disabled.");
#endif
  }  // end standard case (no product map)
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(resultMap));
  return resultMap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetra(Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> v, const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  using ThyMultVecBase     = Thyra::MultiVectorBase<Scalar>;
  using ThyProdMultVecBase = Thyra::ProductMultiVectorBase<Scalar>;
  using ThyUtils           = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  // return value
  RCP<MultiVector> xpMultVec = Teuchos::null;

  // check whether v is a product multi vector
  Teuchos::RCP<const ThyProdMultVecBase> thyProdVec = rcp_dynamic_cast<const ThyProdMultVecBase>(v);
  if (thyProdVec != Teuchos::null) {
    // SPECIAL CASE: create a nested BlockedMultiVector
    // generate nested BlockedMap (containing Thyra and Xpetra GIDs)
    RCP<Map> fullMap = ThyUtils::toXpetra(v->range(), comm);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(Teuchos::rcp_dynamic_cast<BlockedMap>(fullMap)));

    // create new Xpetra::BlockedMultiVector
    xpMultVec = MultiVectorFactory::Build(fullMap, as<size_t>(thyProdVec->domain()->dim()));

    RCP<BlockedMultiVector> xpBlockedMultVec = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(xpMultVec, true);

    // loop over all blocks, transform Thyra MultiVectors to Xpetra MultiVectors recursively
    for (int b = 0; b < thyProdVec->productSpace()->numBlocks(); ++b) {
      RCP<const ThyMultVecBase> thyBlockMV = thyProdVec->getMultiVectorBlock(b);
      // xpBlockMV can be of type MultiVector or BlockedMultiVector
      RCP<const MultiVector> xpBlockMV = ThyUtils::toXpetra(thyBlockMV, comm);  // recursive call
      xpBlockedMultVec->setMultiVector(b, xpBlockMV, true /* Thyra mode */);
    }
  } else {
    // STANDARD CASE: no product vector
#ifdef HAVE_XPETRA_TPETRA
    typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> ConverterT;
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpMultVec;
    typedef Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpTpMultVec;
    typedef Thyra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyTpMultVec;
    typedef Thyra::SpmdMultiVectorBase<Scalar> ThySpmdMultVecBase;

    RCP<const ThySpmdMultVecBase> thyraSpmdMultiVector = rcp_dynamic_cast<const ThySpmdMultVecBase>(v);
    RCP<const ThyTpMultVec> thyraTpetraMultiVector     = rcp_dynamic_cast<const ThyTpMultVec>(thyraSpmdMultiVector);
    TEUCHOS_TEST_FOR_EXCEPTION(thyraTpetraMultiVector == Teuchos::null, Xpetra::Exceptions::RuntimeError, "toXpetra failed to convert provided multi vector to Thyra::TpetraMultiVector. This is the general implementation for Tpetra only.");
    const RCP<const TpMultVec> tpMultVec = ConverterT::getConstTpetraMultiVector(v);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpMultVec));
    RCP<TpMultVec> tpNonConstMultVec = Teuchos::rcp_const_cast<TpMultVec>(tpMultVec);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpNonConstMultVec));
    xpMultVec = rcp(new XpTpMultVec(tpNonConstMultVec));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Cannot transform Thyra::MultiVector to Xpetra::MultiVector. This is the general implementation for Tpetra only, but Teptra is disabled.");
#endif
  }  // end standard case
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpMultVec));
  return xpMultVec;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetra(Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> v, const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> cv =
      Teuchos::rcp_const_cast<const Thyra::MultiVectorBase<Scalar>>(v);
  Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> r =
      toXpetra(cv, comm);
  return Teuchos::rcp_const_cast<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(r);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    isTpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
  TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(op));

  // check whether we have a Tpetra based Thyra operator
  bool bIsTpetra = false;
#ifdef HAVE_XPETRA_TPETRA
  Teuchos::RCP<const Thyra::TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetra_op = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(op);
  bIsTpetra                                                                                      = Teuchos::is_null(tpetra_op) ? false : true;

  // for debugging purposes: find out why dynamic cast failed
  if (!bIsTpetra &&
#ifdef HAVE_XPETRA_EPETRA
      Teuchos::rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op) == Teuchos::null &&
#endif
      Teuchos::rcp_dynamic_cast<const Thyra::DefaultBlockedLinearOp<Scalar>>(op) == Teuchos::null) {
    // If op is not blocked and not an Epetra object, it should be in fact an Tpetra object
    typedef Thyra::TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraLinearOp_t;
    if (Teuchos::rcp_dynamic_cast<const TpetraLinearOp_t>(op) == Teuchos::null) {
      std::cout << "ATTENTION: The dynamic cast to the TpetraLinearOp failed even though it should be a TpetraLinearOp." << std::endl;
      std::cout << "           If you are using Panzer or Stratimikos you might check that the template parameters are " << std::endl;
      std::cout << "           properly set!" << std::endl;
      std::cout << Teuchos::rcp_dynamic_cast<const TpetraLinearOp_t>(op, true) << std::endl;
    }
  }
#endif

#if 0
  // Check whether it is a blocked operator.
  // If yes, grab the (0,0) block and check the underlying linear algebra
  // Note: we expect that the (0,0) block exists!
  if(bIsTpetra == false) {
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> > ThyBlockedOp =
      Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar> >(op);
    if(ThyBlockedOp != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPT(ThyBlockedOp->blockExists(0,0)==false);
      Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > b00 =
        ThyBlockedOp->getBlock(0,0);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(b00));
      bIsTpetra = isTpetra(b00);
    }
  }
#endif

  return bIsTpetra;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    isEpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    isBlockedOperator(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
  // Check whether it is a blocked operator.
  Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar>> ThyBlockedOp =
      Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar>>(op);
  if (ThyBlockedOp != Teuchos::null) {
    return true;
  }
  return false;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
#ifdef HAVE_XPETRA_TPETRA
  if (isTpetra(op)) {
    typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
    Teuchos::RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraOp = TOE::getConstTpetraOperator(op);
    // we should also add support for the const versions!
    // getConstTpetraOperator(const RCP<const LinearOpBase<Scalar> > &op);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraOp));
    Teuchos::RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraRowMat = Teuchos::rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraOp, true);
    Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMat = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraRowMat, true);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraNcnstCrsMat  = Teuchos::rcp_const_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraCrsMat);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraNcnstCrsMat));

    Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpetraCrsMat =
        Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(TpetraNcnstCrsMat));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpetraCrsMat));

    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsMat =
        Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xTpetraCrsMat, true);
    Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsWrap =
        Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpCrsMat));
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpMat =
        Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpCrsWrap, true);
    return xpMat;
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  if (isEpetra(op)) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  }
#endif
  return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetra(const Teuchos::RCP<Thyra::LinearOpBase<Scalar>>& op) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;

#ifdef HAVE_XPETRA_TPETRA
  if (isTpetra(op)) {
    typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
    typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraOperator_t;
    RCP<const TpetraOperator_t> TpetraOp = TOE::getConstTpetraOperator(op);
    TEUCHOS_TEST_FOR_EXCEPT(is_null(TpetraOp));
    RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraRowMat =
        rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraOp, true);
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMat =
        rcp_dynamic_cast<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraRowMat, true);

    RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpetraCrsMat =
        rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
            rcp_const_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraCrsMat)));
    TEUCHOS_TEST_FOR_EXCEPT(is_null(xTpetraCrsMat));
    RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsMat =
        rcp_dynamic_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xTpetraCrsMat, true);
    RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsWrap =
        rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpCrsMat));
    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpMat =
        rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpCrsWrap, true);
    return xpMat;
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  if (isEpetra(op)) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  }
#endif
  return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetraOperator(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
#ifdef HAVE_XPETRA_TPETRA
  if (isTpetra(op)) {
    typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
    Teuchos::RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraOp = TOE::getConstTpetraOperator(op);

    auto nonConstTpetraOp = Teuchos::rcp_const_cast<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraOp);

    Teuchos::RCP<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpetraOp =
        Teuchos::rcp(new Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(nonConstTpetraOp));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpetraOp));

    Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpOp =
        Teuchos::rcp_dynamic_cast<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xTpetraOp, true);
    return xpOp;
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  if (isEpetra(op)) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  }
#endif
  return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetraOperator(const Teuchos::RCP<Thyra::LinearOpBase<Scalar>>& op) {
#ifdef HAVE_XPETRA_TPETRA
  if (isTpetra(op)) {
    typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
    Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraOp = TOE::getTpetraOperator(op);

    Teuchos::RCP<Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpetraOp =
        Teuchos::rcp(new Xpetra::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(TpetraOp));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpetraOp));

    Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpOp =
        Teuchos::rcp_dynamic_cast<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xTpetraOp, true);
    return xpOp;
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  if (isEpetra(op)) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  }
#endif
  return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetra(const Teuchos::RCP<Thyra::DiagonalLinearOpBase<Scalar>>& op) {
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;

  RCP<const Thyra::VectorBase<Scalar>> diag = op->getDiag();

  RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpDiag;
#ifdef HAVE_XPETRA_TPETRA
  using thyTpV = Thyra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using tV     = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  if (!rcp_dynamic_cast<const thyTpV>(diag).is_null()) {
    RCP<const tV> tDiag = Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getConstTpetraVector(diag);
    if (!tDiag.is_null())
      xpDiag = Xpetra::toXpetra(tDiag);
  }
#endif
  TEUCHOS_ASSERT(!xpDiag.is_null());
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> M = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xpDiag);
  return M;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toXpetra(const Teuchos::RCP<const Thyra::DiagonalLinearOpBase<Scalar>>& op) {
  return toXpetra(Teuchos::rcp_const_cast<Thyra::DiagonalLinearOpBase<Scalar>>(op));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toThyra(Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map) {
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> thyraMap = Teuchos::null;

  // check whether map is of type BlockedMap
  RCP<const BlockedMap> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap>(map);
  if (bmap.is_null() == false) {
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>> vecSpaces(bmap->getNumMaps());
    for (size_t i = 0; i < bmap->getNumMaps(); i++) {
      // we need Thyra GIDs for all the submaps
      Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>> vs =
          Xpetra::ThyraUtils<Scalar, LO, GO, Node>::toThyra(bmap->getMap(i, true));
      vecSpaces[i] = vs;
    }

    thyraMap = Thyra::productVectorSpace<Scalar>(vecSpaces());
    return thyraMap;
  }

  // standard case
#ifdef HAVE_XPETRA_TPETRA
  if (map->lib() == Xpetra::UseTpetra) {
    Teuchos::RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>> tpetraMap = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
    if (tpetraMap == Teuchos::null)
      throw Exceptions::BadCast("Xpetra::ThyraUtils::toThyra: Cast from Xpetra::Map to Xpetra::TpetraMap failed");
    RCP<Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>> thyraTpetraMap = Thyra::tpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpetraMap->getTpetra_Map());
    thyraMap                                                                                = thyraTpetraMap;
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  if (map->lib() == Xpetra::UseEpetra) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  }
#endif

  return thyraMap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toThyraMultiVector(Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec) {
  // create Thyra MultiVector
#ifdef HAVE_XPETRA_TPETRA
  if (vec->getMap()->lib() == Xpetra::UseTpetra) {
    auto thyTpMap                                                            = Thyra::tpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Teuchos::rcp_dynamic_cast<const TpetraMap>(vec->getMap())->getTpetra_Map());
    RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpMV = Teuchos::rcp_dynamic_cast<const TpetraMultiVector>(vec)->getTpetra_MultiVector();
    auto thyDomMap                                                           = Thyra::tpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Tpetra::createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(vec->getNumVectors(), vec->getMap()->getComm()));
    auto thyMV                                                               = rcp(new Thyra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
    thyMV->initialize(thyTpMap, thyDomMap, tpMV);
    return thyMV;
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  if (vec->getMap()->lib() == Xpetra::UseEpetra) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  }
#endif

  TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "MultiVector cannot be converted to Thyra.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toThyraVector(Teuchos::RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec) {
  // create Thyra Vector
#ifdef HAVE_XPETRA_TPETRA
  if (vec->getMap()->lib() == Xpetra::UseTpetra) {
    auto thyTpMap                                                        = Thyra::tpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Teuchos::rcp_dynamic_cast<const TpetraMap>(vec->getMap())->getTpetra_Map());
    RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpVec = Teuchos::rcp_dynamic_cast<const TpetraVector>(vec)->getTpetra_Vector();
    auto thyVec                                                          = rcp(new Thyra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>());
    thyVec->initialize(thyTpMap, tpVec);
    return thyVec;
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  if (vec->getMap()->lib() == Xpetra::UseEpetra) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  }
#endif

  TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Vector cannot be converted to Thyra.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    updateThyra(Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> source, Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mapExtractor, const Teuchos::RCP<Thyra::MultiVectorBase<Scalar>>& target) {
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::VectorSpaceBase<Scalar> ThyVecSpaceBase;
  typedef Thyra::SpmdVectorSpaceBase<Scalar> ThySpmdVecSpaceBase;
  typedef Thyra::MultiVectorBase<Scalar> ThyMultVecBase;
  // typedef Thyra::SpmdMultiVectorBase<Scalar> ThySpmdMultVecBase;
  // typedef Thyra::ProductVectorSpaceBase<Scalar> ThyProdVecSpaceBase;
  typedef Thyra::ProductMultiVectorBase<Scalar> ThyProdMultVecBase;

  // copy data from tY_inout to Y_inout
  RCP<ThyProdMultVecBase> prodTarget = rcp_dynamic_cast<ThyProdMultVecBase>(target);
  if (prodTarget != Teuchos::null) {
    RCP<const BlockedMultiVector> bSourceVec = rcp_dynamic_cast<const BlockedMultiVector>(source);
    if (bSourceVec.is_null() == true) {
      // SPECIAL CASE: target vector is product vector:
      // update Thyra product multi vector with data from a merged Xpetra multi vector

      TEUCHOS_TEST_FOR_EXCEPTION(mapExtractor == Teuchos::null, std::logic_error, "Found a Thyra product vector, but user did not provide an Xpetra::MapExtractor.");
      TEUCHOS_TEST_FOR_EXCEPTION(prodTarget->productSpace()->numBlocks() != as<int>(mapExtractor->NumMaps()), std::logic_error, "Inconsistent numbers of sub maps in Thyra::ProductVectorSpace and Xpetra::MapExtractor.");

      for (int bbb = 0; bbb < prodTarget->productSpace()->numBlocks(); ++bbb) {
        // access Xpetra data
        RCP<MultiVector> xpSubBlock = mapExtractor->ExtractVector(source, bbb, false);  // use Xpetra ordering (doesn't really matter)

        // access Thyra data
        Teuchos::RCP<ThyMultVecBase> thySubBlock = prodTarget->getNonconstMultiVectorBlock(bbb);
        RCP<const ThyVecSpaceBase> vs            = thySubBlock->range();
        RCP<const ThySpmdVecSpaceBase> mpi_vs    = rcp_dynamic_cast<const ThySpmdVecSpaceBase>(vs);
        const LocalOrdinal localOffset           = (mpi_vs != Teuchos::null ? mpi_vs->localOffset() : 0);
        const LocalOrdinal localSubDim           = (mpi_vs != Teuchos::null ? mpi_vs->localSubDim() : vs->dim());
        RCP<Thyra::DetachedMultiVectorView<Scalar>> thyData =
            Teuchos::rcp(new Thyra::DetachedMultiVectorView<Scalar>(*thySubBlock, Teuchos::Range1D(localOffset, localOffset + localSubDim - 1)));

        // loop over all vectors in multivector
        for (size_t j = 0; j < xpSubBlock->getNumVectors(); ++j) {
          Teuchos::ArrayRCP<const Scalar> xpData = xpSubBlock->getData(j);  // access const data from Xpetra object

          // loop over all local rows
          for (LocalOrdinal i = 0; i < localSubDim; ++i) {
            (*thyData)(i, j) = xpData[i];
          }
        }
      }
    } else {
      // source vector is a blocked multivector
      // TODO test me
      TEUCHOS_TEST_FOR_EXCEPTION(prodTarget->productSpace()->numBlocks() != as<int>(bSourceVec->getBlockedMap()->getNumMaps()), std::logic_error, "Inconsistent numbers of sub maps in Thyra::ProductVectorSpace and Xpetra::BlockedMultiVector.");

      for (int bbb = 0; bbb < prodTarget->productSpace()->numBlocks(); ++bbb) {
        // access Thyra data
        RCP<MultiVector> xpSubBlock = bSourceVec->getMultiVector(bbb, true);  // use Thyra ordering

        Teuchos::RCP<const ThyMultVecBase> thyXpSubBlock = toThyraMultiVector(xpSubBlock);

        // access Thyra data
        Teuchos::RCP<ThyMultVecBase> thySubBlock = prodTarget->getNonconstMultiVectorBlock(bbb);
        Thyra::assign(thySubBlock.ptr(), *thyXpSubBlock);
      }
    }
  } else {
    // STANDARD case:
    // update Thyra::MultiVector with data from an Xpetra::MultiVector

    // access Thyra data
    RCP<const ThySpmdVecSpaceBase> mpi_vs = rcp_dynamic_cast<const ThySpmdVecSpaceBase>(target->range());
    TEUCHOS_TEST_FOR_EXCEPTION(mpi_vs == Teuchos::null, std::logic_error, "Failed to cast Thyra::VectorSpaceBase to Thyra::SpmdVectorSpaceBase.");
    const LocalOrdinal localOffset = (mpi_vs != Teuchos::null ? mpi_vs->localOffset() : 0);
    const LocalOrdinal localSubDim = (mpi_vs != Teuchos::null ? mpi_vs->localSubDim() : target->range()->dim());
    RCP<Thyra::DetachedMultiVectorView<Scalar>> thyData =
        Teuchos::rcp(new Thyra::DetachedMultiVectorView<Scalar>(*target, Teuchos::Range1D(localOffset, localOffset + localSubDim - 1)));

    // loop over all vectors in multivector
    for (size_t j = 0; j < source->getNumVectors(); ++j) {
      Teuchos::ArrayRCP<const Scalar> xpData = source->getData(j);  // access const data from Xpetra object
      // loop over all local rows
      for (LocalOrdinal i = 0; i < localSubDim; ++i) {
        (*thyData)(i, j) = xpData[i];
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toThyra(const Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat) {
  // create a Thyra operator from Xpetra::CrsMatrix
  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> thyraOp = Teuchos::null;

  // bool bIsTpetra = false;

#ifdef HAVE_XPETRA_TPETRA
  Teuchos::RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetraMat = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat);
  if (tpetraMat != Teuchos::null) {
    // bIsTpetra = true;
    Teuchos::RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpCrsMat = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat, true);
    Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpCrsMat        = xTpCrsMat->getTpetra_CrsMatrix();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));

    Teuchos::RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpRowMat  = Teuchos::rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpCrsMat, true);
    Teuchos::RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpOperator = Teuchos::rcp_dynamic_cast<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpRowMat, true);

    thyraOp = Thyra::createConstLinearOp(tpOperator);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraOp));
  } else {
#ifdef HAVE_XPETRA_EPETRA
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Cast to Tpetra::CrsMatrix failed. Assume matrix should be Epetra then. Epetra needs SC=double, LO=int, and GO=int or GO=long long");
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Cast to Tpetra::CrsMatrix failed. Assume matrix should be Epetra then. No Epetra available");
#endif
  }
  return thyraOp;
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
#endif
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toThyra(const Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat) {
  // create a Thyra operator from Xpetra::CrsMatrix
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>> thyraOp = Teuchos::null;

  // bool bIsTpetra = false;

#ifdef HAVE_XPETRA_TPETRA
  Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat);
  if (tpetraMat != Teuchos::null) {
    // bIsTpetra = true;
    Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpCrsMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat, true);
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpCrsMat        = xTpCrsMat->getTpetra_CrsMatrixNonConst();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));

    Teuchos::RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpRowMat  = Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpCrsMat, true);
    Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpOperator = Teuchos::rcp_dynamic_cast<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpRowMat, true);

    thyraOp = Thyra::createLinearOp(tpOperator);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraOp));
  } else {
    // cast to TpetraCrsMatrix failed
#ifdef HAVE_XPETRA_EPETRA
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Cast to TpetraCrsMatrix failed. Assuming matrix supposed to be Epetra. Epetra needs SC=double, LO=int, and GO=int or GO=long long");
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Cast to TpetraCrsMatrix failed. Guess, matrix should be Epetra then, but no Epetra available.");
#endif
  }
  return thyraOp;
#else
  TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
#endif
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    toThyra(const Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat) {
  int nRows = mat->Rows();
  int nCols = mat->Cols();

  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Ablock             = mat->getInnermostCrsMatrix();
  Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> Ablock_wrap = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Ablock);
  TEUCHOS_TEST_FOR_EXCEPT(Ablock_wrap.is_null() == true);

#ifdef HAVE_XPETRA_TPETRA
  Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Ablock_wrap->getCrsMatrix());
  if (tpetraMat != Teuchos::null) {
    // create new Thyra blocked operator
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar>> blockMat =
        Thyra::defaultBlockedLinearOp<Scalar>();

    blockMat->beginBlockFill(nRows, nCols);

    for (int r = 0; r < nRows; ++r) {
      for (int c = 0; c < nCols; ++c) {
        Teuchos::RCP<Matrix> xpmat = mat->getMatrix(r, c);

        if (xpmat == Teuchos::null) continue;  // shortcut for empty blocks

        Teuchos::RCP<Thyra::LinearOpBase<Scalar>> thBlock = Teuchos::null;

        // check whether the subblock is again a blocked operator
        Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpblock =
            Teuchos::rcp_dynamic_cast<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpmat);
        if (xpblock != Teuchos::null) {
          if (xpblock->Rows() == 1 && xpblock->Cols() == 1) {
            // If it is a single block operator, unwrap it
            Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpblock->getCrsMatrix());
            TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
            thBlock = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpwrap->getCrsMatrix());
          } else {
            // recursive call for general blocked operators
            thBlock = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpblock);
          }
        } else {
          // check whether it is a CRSMatrix object
          Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpmat);
          TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
          thBlock = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpwrap->getCrsMatrix());
        }

        blockMat->setBlock(r, c, thBlock);
      }
    }

    blockMat->endBlockFill();

    return blockMat;
  } else {
    // tpetraMat == Teuchos::null
#ifdef HAVE_XPETRA_EPETRA
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Cast to TpetraCrsMatrix failed. Assuming matrix supposed to be Epetra. Epetra needs SC=double, LO=int, and GO=int or GO=long long");
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Cast to TpetraCrsMatrix failed. Guess, matrix should be Epetra then, but no Epetra available.");
#endif
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }
#endif  // endif HAVE_XPETRA_TPETRA

  // 4-Aug-2017 JJH Added 2nd condition to avoid "warning: dynamic initialization in unreachable code"
  //                If HAVE_XPETRA_TPETRA is defined, then this method will always return or throw in the if-then-else above.
#if defined(HAVE_XPETRA_EPETRA) && !defined(HAVE_XPETRA_TPETRA)
  TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
#endif  // endif HAVE_XPETRA_EPETRA
}

#ifdef HAVE_XPETRA_EPETRA

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
// implementation of "toThyra" for full specialization on SC=double, LO=GO=int and NO=EpetraNode
// We need the specialization in the cpp file due to a circle dependency in the .hpp files for BlockedCrsMatrix
Teuchos::RCP<Thyra::LinearOpBase<double>>
ThyraUtils<double, int, int, EpetraNode>::toThyra(const Teuchos::RCP<Xpetra::BlockedCrsMatrix<double, int, int, EpetraNode>>& mat) {
  int nRows = mat->Rows();
  int nCols = mat->Cols();

  Teuchos::RCP<Xpetra::Matrix<double, int, int, EpetraNode>> Ablock             = mat->getInnermostCrsMatrix();
  Teuchos::RCP<Xpetra::CrsMatrixWrap<double, int, int, EpetraNode>> Ablock_wrap = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<double, int, int, EpetraNode>>(Ablock);
  TEUCHOS_TEST_FOR_EXCEPT(Ablock_wrap.is_null() == true);

  bool bTpetra = false;
  bool bEpetra = false;
#ifdef HAVE_XPETRA_TPETRA
  // Note: Epetra is enabled
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))
  Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Ablock_wrap->getCrsMatrix());
  if (tpetraMat != Teuchos::null) bTpetra = true;
#else
  bTpetra = false;
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
  Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> epetraMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(Ablock_wrap->getCrsMatrix());
  if (epetraMat != Teuchos::null) bEpetra = true;
#endif

  TEUCHOS_TEST_FOR_EXCEPT(bTpetra == bEpetra);  // we only allow Epetra OR Tpetra

  // create new Thyra blocked operator
  Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar>> blockMat =
      Thyra::defaultBlockedLinearOp<Scalar>();

  blockMat->beginBlockFill(nRows, nCols);

  for (int r = 0; r < nRows; ++r) {
    for (int c = 0; c < nCols; ++c) {
      Teuchos::RCP<Matrix> xpmat = mat->getMatrix(r, c);

      if (xpmat == Teuchos::null) continue;  // shortcut for empty blocks

      Teuchos::RCP<Thyra::LinearOpBase<Scalar>> thBlock = Teuchos::null;

      // check whether the subblock is again a blocked operator
      Teuchos::RCP<BlockedCrsMatrix> xpblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(xpmat);
      if (xpblock != Teuchos::null) {
        if (xpblock->Rows() == 1 && xpblock->Cols() == 1) {
          // If it is a single block operator, unwrap it
          Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpblock->getCrsMatrix());
          TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
          thBlock = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpwrap->getCrsMatrix());
        } else {
          // recursive call for general blocked operators
          thBlock = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpblock);
        }
      } else {
        // check whether it is a CRSMatrix object
        Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpmat);
        TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
        thBlock = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpwrap->getCrsMatrix());
      }

      blockMat->setBlock(r, c, thBlock);
    }
  }

  blockMat->endBlockFill();

  return blockMat;
}
#endif  // #ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
// implementation of "toThyra" for full specialization on SC=double, LO=int, GO=long long and NO=EpetraNode
// We need the specialization in the cpp file due to a circle dependency in the .hpp files for BlockedCrsMatrix
Teuchos::RCP<Thyra::LinearOpBase<double>>
ThyraUtils<double, int, long long, EpetraNode>::toThyra(const Teuchos::RCP<Xpetra::BlockedCrsMatrix<double, int, long long, EpetraNode>>& mat) {
  int nRows = mat->Rows();
  int nCols = mat->Cols();

  Teuchos::RCP<Xpetra::Matrix<double, int, long long, EpetraNode>> Ablock             = mat->getInnermostCrsMatrix();
  Teuchos::RCP<Xpetra::CrsMatrixWrap<double, int, long long, EpetraNode>> Ablock_wrap = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<double, int, long long, EpetraNode>>(Ablock);
  TEUCHOS_TEST_FOR_EXCEPT(Ablock_wrap.is_null() == true);

  bool bTpetra = false;
  bool bEpetra = false;
#ifdef HAVE_XPETRA_TPETRA
  // Note: Epetra is enabled
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))
  Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(Ablock_wrap->getCrsMatrix());
  if (tpetraMat != Teuchos::null) bTpetra = true;
#else
  bTpetra = false;
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
  Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> epetraMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(Ablock_wrap->getCrsMatrix());
  if (epetraMat != Teuchos::null) bEpetra = true;
#endif

  TEUCHOS_TEST_FOR_EXCEPT(bTpetra == bEpetra);  // we only allow Epetra OR Tpetra

  // create new Thyra blocked operator
  Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar>> blockMat =
      Thyra::defaultBlockedLinearOp<Scalar>();

  blockMat->beginBlockFill(nRows, nCols);

  for (int r = 0; r < nRows; ++r) {
    for (int c = 0; c < nCols; ++c) {
      Teuchos::RCP<Matrix> xpmat = mat->getMatrix(r, c);

      if (xpmat == Teuchos::null) continue;  // shortcut for empty blocks

      Teuchos::RCP<Thyra::LinearOpBase<Scalar>> thBlock = Teuchos::null;

      // check whether the subblock is again a blocked operator
      Teuchos::RCP<BlockedCrsMatrix> xpblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(xpmat);
      if (xpblock != Teuchos::null) {
        if (xpblock->Rows() == 1 && xpblock->Cols() == 1) {
          // If it is a single block operator, unwrap it
          Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpblock->getCrsMatrix());
          TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
          thBlock = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpwrap->getCrsMatrix());
        } else {
          // recursive call for general blocked operators
          thBlock = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpblock);
        }
      } else {
        // check whether it is a CRSMatrix object
        Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpmat);
        TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
        thBlock = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpwrap->getCrsMatrix());
      }

      blockMat->setBlock(r, c, thBlock);
    }
  }

  blockMat->endBlockFill();

  return blockMat;
}
#endif  // #ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES

#endif

}  // namespace Xpetra

#endif

#endif
