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

#ifdef HAVE_XPETRA_TPETRA
#include "Tpetra_ConfigDefs.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Epetra_config.h"
#include "Epetra_CombineMode.h"
#endif

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

#ifdef HAVE_XPETRA_TPETRA
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
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_SpmdVectorBase.hpp>
#include <Thyra_get_Epetra_Operator.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif

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

  static bool isEpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op);

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

// full specialization for Epetra support
// Note, that Thyra only has support for Epetra (not for Epetra64)
#ifdef HAVE_XPETRA_EPETRA

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template <>
class ThyraUtils<double, int, int, EpetraNode> {
 public:
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

 private:
#undef XPETRA_THYRAUTILS_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  static Teuchos::RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int>>& comm, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId = -1, GlobalOrdinal offset = 0) {
    Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map = ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toXpetra(vectorSpace, comm);

    if (stridedBlockId == -1) {
      TEUCHOS_TEST_FOR_EXCEPT(map->getLocalNumElements() % stridingInfo.size() != 0);
    } else {
      TEUCHOS_TEST_FOR_EXCEPT(map->getLocalNumElements() % stridingInfo[stridedBlockId] != 0);
    }

    Teuchos::RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>> ret = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, stridingInfo, stridedBlockId, offset);
    return ret;
  }

  static Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    typedef Thyra::VectorSpaceBase<Scalar> ThyVecSpaceBase;
    typedef Thyra::ProductVectorSpaceBase<Scalar> ThyProdVecSpaceBase;
    typedef Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyUtils;

    RCP<Map> resultMap = Teuchos::null;

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
      // STANDARD CASE: no product map
      // Epetra/Tpetra specific code to access the underlying map data

      // check whether we have a Tpetra based Thyra operator
      bool bIsTpetra = false;
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))
      Teuchos::RCP<const Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetra_vsc = Teuchos::rcp_dynamic_cast<const Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(vectorSpace);
      bIsTpetra                                                                                          = Teuchos::is_null(tpetra_vsc) ? false : true;
#endif
#endif

      // check whether we have an Epetra based Thyra operator
      bool bIsEpetra = !bIsTpetra;  // note: this is a little bit fragile!

#ifdef HAVE_XPETRA_TPETRA
      if (bIsTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))
        typedef Thyra::VectorBase<Scalar> ThyVecBase;
        typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> TpMap;
        typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpVector;
        typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
        RCP<ThyVecBase> rgVec = Thyra::createMember<Scalar>(vectorSpace, std::string("label"));
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgVec));
        RCP<const TpVector> rgTpetraVec = TOE::getTpetraVector(rgVec);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgTpetraVec));
        RCP<const TpMap> rgTpetraMap = rgTpetraVec->getMap();
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgTpetraMap));

        resultMap = Xpetra::toXpetraNonConst(rgTpetraMap);
#else
        throw Xpetra::Exceptions::RuntimeError("Problem AAA. Add TPETRA_INST_INT_INT:BOOL=ON in your configuration.");
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (bIsEpetra) {
        // RCP<const Epetra_Map> epMap = Teuchos::null;
        RCP<const Epetra_Map>
            epetra_map = Teuchos::get_extra_data<RCP<const Epetra_Map>>(vectorSpace, "epetra_map");
        if (!Teuchos::is_null(epetra_map)) {
          resultMap = Teuchos::rcp(new Xpetra::EpetraMapT<GlobalOrdinal, Node>(epetra_map));
          TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(resultMap));
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No Epetra_Map data found in Thyra::VectorSpace.");
        }
      }
#endif
    }  // end standard case (no product map)
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(resultMap));
    return resultMap;
  }

  // const version
  static Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> v, const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    using ThyProdMultVecBase = Thyra::ProductMultiVectorBase<Scalar>;
    using ThyMultVecBase     = Thyra::MultiVectorBase<Scalar>;
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

      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpMultVec));
      return xpMultVec;
    } else {
      // STANDARD CASE: no product vector
      // Epetra/Tpetra specific code to access the underlying map data
      bool bIsTpetra = false;
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))

      // typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> TpMap;
      // typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpVector;
      typedef Thyra::SpmdMultiVectorBase<Scalar> ThySpmdMultVecBase;
      typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> ConverterT;
      typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpMultVec;
      typedef Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpTpMultVec;
      typedef Thyra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyTpMultVec;

      RCP<const ThySpmdMultVecBase> thyraSpmdMultiVector = rcp_dynamic_cast<const ThySpmdMultVecBase>(v);
      RCP<const ThyTpMultVec> thyraTpetraMultiVector     = rcp_dynamic_cast<const ThyTpMultVec>(thyraSpmdMultiVector);
      if (thyraTpetraMultiVector != Teuchos::null) {
        bIsTpetra                            = true;
        const RCP<const TpMultVec> tpMultVec = ConverterT::getConstTpetraMultiVector(v);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpMultVec));
        RCP<TpMultVec> tpNonConstMultVec = Teuchos::rcp_const_cast<TpMultVec>(tpMultVec);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpNonConstMultVec));
        xpMultVec = rcp(new XpTpMultVec(tpNonConstMultVec));
      }
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (bIsTpetra == false) {
        // no product vector
        Teuchos::RCP<Map> map = ThyUtils::toXpetra(v->range(), comm);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(map));
        RCP<Xpetra::EpetraMapT<GlobalOrdinal, Node>> xeMap = rcp_dynamic_cast<Xpetra::EpetraMapT<GlobalOrdinal, Node>>(map, true);
        RCP<const Epetra_Map> eMap                         = xeMap->getEpetra_MapRCP();
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(eMap));
        Teuchos::RCP<const Epetra_MultiVector> epMultVec = Thyra::get_Epetra_MultiVector(*eMap, v);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epMultVec));
        RCP<Epetra_MultiVector> epNonConstMultVec = Teuchos::rcp_const_cast<Epetra_MultiVector>(epMultVec);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epNonConstMultVec));
        xpMultVec = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>(epNonConstMultVec));
      }
#endif
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpMultVec));
      return xpMultVec;
    }  // end standard case
  }

  // non-const version
  static Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> v, const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> cv =
        Teuchos::rcp_const_cast<const Thyra::MultiVectorBase<Scalar>>(v);
    Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> r =
        toXpetra(cv, comm);
    return Teuchos::rcp_const_cast<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(r);
  }

  static bool isTpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
    // check whether we have a Tpetra based Thyra operator
    bool bIsTpetra = false;
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))

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

  static bool isEpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
    // check whether we have an Epetra based Thyra operator
    bool bIsEpetra = false;

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<const Thyra::EpetraLinearOp> epetra_op = Teuchos::rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op, false);
    bIsEpetra                                           = Teuchos::is_null(epetra_op) ? false : true;
#endif

#if 0
    // Check whether it is a blocked operator.
    // If yes, grab the (0,0) block and check the underlying linear algebra
    // Note: we expect that the (0,0) block exists!
    if(bIsEpetra == false) {
      Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> > ThyBlockedOp =
          Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar> >(op,false);
      if(ThyBlockedOp != Teuchos::null) {
        TEUCHOS_TEST_FOR_EXCEPT(ThyBlockedOp->blockExists(0,0)==false);
        Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > b00 =
          ThyBlockedOp->getBlock(0,0);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(b00));
        bIsEpetra = isEpetra(b00);
      }
    }
#endif

    return bIsEpetra;
  }

  static bool isBlockedOperator(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
    // Check whether it is a blocked operator.
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar>> ThyBlockedOp =
        Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar>>(op);
    if (ThyBlockedOp != Teuchos::null) {
      return true;
    }
    return false;
  }

  static Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
#ifdef HAVE_XPETRA_TPETRA
    if (isTpetra(op)) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))

      typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
      Teuchos::RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraOp = TOE::getConstTpetraOperator(op);
      // we should also add support for the const versions!
      // getConstTpetraOperator(const RCP<const LinearOpBase<Scalar> > &op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraOp));
      Teuchos::RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraRowMat = Teuchos::rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraOp);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraRowMat));
      Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMat = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraRowMat, true);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraNcnstCrsMat  = Teuchos::rcp_const_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraCrsMat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraNcnstCrsMat));

      Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpetraCrsMat =
          Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(TpetraNcnstCrsMat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpetraCrsMat));
      Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ret =
          Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xTpetraCrsMat);
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsMat =
          Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xTpetraCrsMat, true);
      Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsWrap =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpCrsMat));
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpMat =
          Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpCrsWrap, true);
      return xpMat;
#else
      throw Xpetra::Exceptions::RuntimeError("Problem BBB. Add TPETRA_INST_INT_INT:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    if (isEpetra(op)) {
      Teuchos::RCP<const Epetra_Operator> epetra_op = Thyra::get_Epetra_Operator(*op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_op));
      Teuchos::RCP<const Epetra_RowMatrix> epetra_rowmat = Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(epetra_op, true);
      Teuchos::RCP<const Epetra_CrsMatrix> epetra_crsmat = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(epetra_rowmat, true);
      Teuchos::RCP<Epetra_CrsMatrix> epetra_ncnstcrsmat  = Teuchos::rcp_const_cast<Epetra_CrsMatrix>(epetra_crsmat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_ncnstcrsmat));

      Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> xEpetraCrsMat =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>(epetra_ncnstcrsmat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xEpetraCrsMat));

      Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsMat =
          Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xEpetraCrsMat, true);
      Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsWrap =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpCrsMat));
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpMat =
          Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpCrsWrap, true);
      return xpMat;
    }
#endif
    return Teuchos::null;
  }

  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<Thyra::LinearOpBase<Scalar>>& op) {
#ifdef HAVE_XPETRA_TPETRA
    if (isTpetra(op)) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))

      typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
      Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraOp = TOE::getTpetraOperator(op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraOp));
      Teuchos::RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraRowMat = Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraOp, true);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMat = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraRowMat, true);

      Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpetraCrsMat =
          Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(TpetraCrsMat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpetraCrsMat));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsMat =
          Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xTpetraCrsMat, true);
      Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsWrap =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpCrsMat));
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpMat =
          Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpCrsWrap, true);
      return xpMat;
#else
      throw Xpetra::Exceptions::RuntimeError("Problem CCC. Add TPETRA_INST_INT_INT:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    if (isEpetra(op)) {
      Teuchos::RCP<Epetra_Operator> epetra_op = Thyra::get_Epetra_Operator(*op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_op));
      Teuchos::RCP<Epetra_RowMatrix> epetra_rowmat = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(epetra_op, true);
      Teuchos::RCP<Epetra_CrsMatrix> epetra_crsmat = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(epetra_rowmat, true);

      Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> xEpetraCrsMat =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>(epetra_crsmat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xEpetraCrsMat));
      Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsMat =
          Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xEpetraCrsMat, true);
      Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsWrap =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpCrsMat));
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpMat =
          Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpCrsWrap, true);
      return xpMat;
    }
#endif
    return Teuchos::null;
  }

  static Teuchos::RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetraOperator(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
    return toXpetraOperator(Teuchos::rcp_const_cast<Thyra::LinearOpBase<Scalar>>(op));

    // #ifdef HAVE_XPETRA_TPETRA
    //     if(isTpetra(op)) {
    //       typedef Thyra::TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node> TOE;
    //       Teuchos::RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > TpetraOp = TOE::getConstTpetraOperator(op);

    //       Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > nonConstTpetraOp =
    //         Teuchos::rcp_const_cast<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(TpetraOp);

    //       Teuchos::RCP<Xpetra::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xTpetraOp =
    //         Teuchos::rcp(new Xpetra::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>(nonConstTpetraOp));
    //       TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpetraOp));

    //       Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xpOp =
    //         Teuchos::rcp_dynamic_cast<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(xTpetraOp, true);
    //       return xpOp;
    //     }
    // #endif

    // #ifdef HAVE_XPETRA_EPETRA
    //     if(isEpetra(op)) {
    //       TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Epetra needs SC=double, LO=int, and GO=int or GO=long long");
    //     }
    // #endif
    //     return Teuchos::null;
  }

  static Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
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

  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
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
#ifdef HAVE_XPETRA_EPETRA
    using ThyVSBase = Thyra::SpmdVectorSpaceBase<Scalar>;
    if (xpDiag.is_null()) {
      RCP<const Epetra_Comm> comm = Thyra::get_Epetra_Comm(*rcp_dynamic_cast<const ThyVSBase>(op->range())->getComm());
      RCP<const Epetra_Map> map   = Thyra::get_Epetra_Map(*(op->range()), comm);
      if (!map.is_null()) {
        RCP<const Epetra_Vector> eDiag                 = Thyra::get_Epetra_Vector(*map, diag);
        RCP<Epetra_Vector> nceDiag                     = rcp_const_cast<Epetra_Vector>(eDiag);
        RCP<Xpetra::EpetraVectorT<int, Node>> xpEpDiag = rcp(new Xpetra::EpetraVectorT<int, Node>(nceDiag));
        xpDiag                                         = rcp_dynamic_cast<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpEpDiag, true);
      }
    }
#endif
    TEUCHOS_ASSERT(!xpDiag.is_null());
    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> M = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xpDiag);
    return M;
  }

  static Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::DiagonalLinearOpBase<Scalar>>& op) {
    return toXpetra(Teuchos::rcp_const_cast<Thyra::DiagonalLinearOpBase<Scalar>>(op));
  }

  static Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
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
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))
      Teuchos::RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>> tpetraMap = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
      if (tpetraMap == Teuchos::null)
        throw Exceptions::BadCast("Xpetra::ThyraUtils::toThyra: Cast from Xpetra::Map to Xpetra::TpetraMap failed");
      RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> tpMap                         = tpetraMap->getTpetra_Map();
      RCP<Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>> thyraTpetraMap = Thyra::tpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpMap);
      thyraMap                                                                                = thyraTpetraMap;
#else
      throw Xpetra::Exceptions::RuntimeError("Problem DDD. Add TPETRA_INST_INT_INT:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    if (map->lib() == Xpetra::UseEpetra) {
      Teuchos::RCP<const Xpetra::EpetraMapT<GlobalOrdinal, Node>> epetraMap = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GlobalOrdinal, Node>>(map);
      if (epetraMap == Teuchos::null)
        throw Exceptions::BadCast("Xpetra::ThyraUtils::toThyra: Cast from Xpetra::Map to Xpetra::EpetraMap failed");
      RCP<const Thyra::VectorSpaceBase<Scalar>> thyraEpetraMap = Thyra::create_VectorSpace(epetraMap->getEpetra_MapRCP());
      thyraMap                                                 = thyraEpetraMap;
    }
#endif

    return thyraMap;
  }

  static Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
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
      auto thyEpMap = Thyra::create_VectorSpace(Teuchos::rcp_dynamic_cast<const EpetraMapT<GlobalOrdinal, Node>>(vec->getMap())->getEpetra_MapRCP());
      auto epMV     = Teuchos::rcp_dynamic_cast<const EpetraMultiVectorT<GlobalOrdinal, Node>>(vec)->getEpetra_MultiVector();
      auto thyMV    = Thyra::create_MultiVector(epMV, thyEpMap);
      return thyMV;
    }
#endif

    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "MultiVector cannot be converted to Thyra.");
  }

  static Teuchos::RCP<const Thyra::VectorBase<Scalar>>
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
      auto thyEpMap = Thyra::create_VectorSpace(Teuchos::rcp_dynamic_cast<const EpetraMapT<GlobalOrdinal, Node>>(vec->getMap())->getEpetra_MapRCP());
      auto epVec    = rcp(Teuchos::rcp_dynamic_cast<const EpetraVectorT<GlobalOrdinal, Node>>(vec)->getEpetra_Vector(), false);
      auto thyVec   = Thyra::create_Vector(epVec, thyEpMap);
      return thyVec;
    }
#endif

    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Vector cannot be converted to Thyra.");
  }

  static void updateThyra(Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> source, Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mapExtractor, const Teuchos::RCP<Thyra::MultiVectorBase<Scalar>>& target) {
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

  static Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>
  toThyra(const Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat) {
    // create a Thyra operator from Xpetra::CrsMatrix
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> thyraOp = Teuchos::null;

#ifdef HAVE_XPETRA_TPETRA
    Teuchos::RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetraMat = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat);
    if (tpetraMat != Teuchos::null) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))

      Teuchos::RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpCrsMat = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat, true);
      Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpCrsMat        = xTpCrsMat->getTpetra_CrsMatrix();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));

      Teuchos::RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpRowMat  = Teuchos::rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpCrsMat, true);
      Teuchos::RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpOperator = Teuchos::rcp_dynamic_cast<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpRowMat, true);

      thyraOp = Thyra::createConstLinearOp(tpOperator);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraOp));
#else
      throw Xpetra::Exceptions::RuntimeError("Add TPETRA_INST_INT_INT:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> epetraMat = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(mat);
    if (epetraMat != Teuchos::null) {
      Teuchos::RCP<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> xEpCrsMat = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(mat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xEpCrsMat));
      Teuchos::RCP<const Epetra_CrsMatrix> epCrsMat = xEpCrsMat->getEpetra_CrsMatrix();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epCrsMat));

      Teuchos::RCP<const Thyra::EpetraLinearOp> thyraEpOp = Thyra::epetraLinearOp(epCrsMat, "op");
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraEpOp));
      thyraOp = thyraEpOp;
    }
#endif
    return thyraOp;
  }

  static Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
  toThyra(const Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat) {
    // create a Thyra operator from Xpetra::CrsMatrix
    Teuchos::RCP<Thyra::LinearOpBase<Scalar>> thyraOp = Teuchos::null;

#ifdef HAVE_XPETRA_TPETRA
    Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat);
    if (tpetraMat != Teuchos::null) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))

      Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpCrsMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat, true);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpCrsMat        = xTpCrsMat->getTpetra_CrsMatrixNonConst();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));

      Teuchos::RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpRowMat  = Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpCrsMat, true);
      Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpOperator = Teuchos::rcp_dynamic_cast<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpRowMat, true);

      thyraOp = Thyra::createLinearOp(tpOperator);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraOp));
#else
      throw Xpetra::Exceptions::RuntimeError("Add TPETRA_INST_INT_INT:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> epetraMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(mat);
    if (epetraMat != Teuchos::null) {
      Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> xEpCrsMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(mat, true);
      Teuchos::RCP<Epetra_CrsMatrix> epCrsMat                               = xEpCrsMat->getEpetra_CrsMatrixNonConst();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epCrsMat));

      Teuchos::RCP<Thyra::EpetraLinearOp> thyraEpOp = Thyra::nonconstEpetraLinearOp(epCrsMat, "op");
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraEpOp));
      thyraOp = thyraEpOp;
    }
#endif
    return thyraOp;
  }

  static Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
  toThyra(const Teuchos::RCP<Xpetra::BlockedCrsMatrix<double, int, int, EpetraNode>>& mat);

};      // specialization on SC=double, LO=GO=int and NO=EpetraNode
#endif  // #ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template <>
class ThyraUtils<double, int, long long, EpetraNode> {
 public:
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

 private:
#undef XPETRA_THYRAUTILS_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  static Teuchos::RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int>>& comm, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId = -1, GlobalOrdinal offset = 0) {
    Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map = ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toXpetra(vectorSpace, comm);

    if (stridedBlockId == -1) {
      TEUCHOS_TEST_FOR_EXCEPT(map->getLocalNumElements() % stridingInfo.size() != 0);
    } else {
      TEUCHOS_TEST_FOR_EXCEPT(map->getLocalNumElements() % stridingInfo[stridedBlockId] != 0);
    }

    Teuchos::RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>> ret = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(map, stridingInfo, stridedBlockId, offset);
    return ret;
  }

  static Teuchos::RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    typedef Thyra::VectorSpaceBase<Scalar> ThyVecSpaceBase;
    typedef Thyra::ProductVectorSpaceBase<Scalar> ThyProdVecSpaceBase;
    typedef Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyUtils;

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

      Teuchos::RCP<Map> resultMap = Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(mapsXpetra, mapsThyra));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(resultMap));
      return resultMap;
    } else {
      // STANDARD CASE: no product map
      // Epetra/Tpetra specific code to access the underlying map data

      // check whether we have a Tpetra based Thyra operator
      bool bIsTpetra = false;
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))
      Teuchos::RCP<const Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetra_vsc = Teuchos::rcp_dynamic_cast<const Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(vectorSpace);
      bIsTpetra                                                                                          = Teuchos::is_null(tpetra_vsc) ? false : true;
#endif
#endif

      // check whether we have an Epetra based Thyra operator
      bool bIsEpetra = !bIsTpetra;  // note: this is a little bit fragile!

#ifdef HAVE_XPETRA_TPETRA
      if (bIsTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))
        typedef Thyra::VectorBase<Scalar> ThyVecBase;
        typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> TpMap;
        typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpVector;
        typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
        RCP<ThyVecBase> rgVec = Thyra::createMember<Scalar>(vectorSpace, std::string("label"));
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgVec));
        RCP<const TpVector> rgTpetraVec = TOE::getTpetraVector(rgVec);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgTpetraVec));
        RCP<const TpMap> rgTpetraMap = rgTpetraVec->getMap();
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgTpetraMap));

        RCP<Map> rgXpetraMap = Xpetra::toXpetraNonConst(rgTpetraMap);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgXpetraMap));
        return rgXpetraMap;
#else
        throw Xpetra::Exceptions::RuntimeError("Add TPETRA_INST_INT_LONG_LONG:BOOL=ON in your configuration.");
#endif
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (bIsEpetra) {
        // RCP<const Epetra_Map> epMap = Teuchos::null;
        RCP<const Epetra_Map>
            epetra_map = Teuchos::get_extra_data<RCP<const Epetra_Map>>(vectorSpace, "epetra_map");
        if (!Teuchos::is_null(epetra_map)) {
          Teuchos::RCP<Map> rgXpetraMap = Teuchos::rcp(new Xpetra::EpetraMapT<GlobalOrdinal, Node>(epetra_map));
          TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgXpetraMap));
          return rgXpetraMap;
        } else {
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No Epetra_Map data found in Thyra::VectorSpace.");
        }
      }
#endif
    }  // end standard case (no product map)
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Cannot transform Thyra::VectorSpace to Xpetra::Map.");
    // return Teuchos::null; // unreachable
  }

  // const version
  static Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> v, const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;
    typedef Thyra::ProductMultiVectorBase<Scalar> ThyProdMultVecBase;
    typedef Thyra::MultiVectorBase<Scalar> ThyMultVecBase;
    typedef Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyUtils;

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

      RCP<BlockedMultiVector> xpBlockedMultVec = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(xpMultVec);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpBlockedMultVec));

      // loop over all blocks, transform Thyra MultiVectors to Xpetra MultiVectors recursively
      for (int b = 0; b < thyProdVec->productSpace()->numBlocks(); ++b) {
        RCP<const ThyMultVecBase> thyBlockMV = thyProdVec->getMultiVectorBlock(b);
        // xpBlockMV can be of type MultiVector or BlockedMultiVector
        RCP<const MultiVector> xpBlockMV = ThyUtils::toXpetra(thyBlockMV, comm);  // recursive call
        xpBlockedMultVec->setMultiVector(b, xpBlockMV, true /* Thyra mode */);
      }

      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpMultVec));
      return xpMultVec;
    } else {
      // STANDARD CASE: no product vector
      // Epetra/Tpetra specific code to access the underlying map data
      bool bIsTpetra = false;
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))

      // typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> TpMap;
      // typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpVector;
      typedef Thyra::SpmdMultiVectorBase<Scalar> ThySpmdMultVecBase;
      typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> ConverterT;
      typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpMultVec;
      typedef Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> XpTpMultVec;
      typedef Thyra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyTpMultVec;

      RCP<const ThySpmdMultVecBase> thyraSpmdMultiVector = rcp_dynamic_cast<const ThySpmdMultVecBase>(v);
      RCP<const ThyTpMultVec> thyraTpetraMultiVector     = rcp_dynamic_cast<const ThyTpMultVec>(thyraSpmdMultiVector);
      if (thyraTpetraMultiVector != Teuchos::null) {
        bIsTpetra                            = true;
        const RCP<const TpMultVec> tpMultVec = ConverterT::getConstTpetraMultiVector(v);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpMultVec));
        RCP<TpMultVec> tpNonConstMultVec = Teuchos::rcp_const_cast<TpMultVec>(tpMultVec);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpNonConstMultVec));
        xpMultVec = rcp(new XpTpMultVec(tpNonConstMultVec));
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpMultVec));
        return xpMultVec;
      }
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (bIsTpetra == false) {
        // no product vector
        Teuchos::RCP<Map> map = ThyUtils::toXpetra(v->range(), comm);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(map));
        RCP<Xpetra::EpetraMapT<GlobalOrdinal, Node>> xeMap = rcp_dynamic_cast<Xpetra::EpetraMapT<GlobalOrdinal, Node>>(map, true);
        RCP<const Epetra_Map> eMap                         = xeMap->getEpetra_MapRCP();
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(eMap));
        Teuchos::RCP<const Epetra_MultiVector> epMultVec = Thyra::get_Epetra_MultiVector(*eMap, v);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epMultVec));
        RCP<Epetra_MultiVector> epNonConstMultVec = Teuchos::rcp_const_cast<Epetra_MultiVector>(epMultVec);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epNonConstMultVec));
        xpMultVec = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>(epNonConstMultVec));
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpMultVec));
        return xpMultVec;
      }
#endif
    }  // end standard case
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Cannot transform Thyra::MultiVector to Xpetra::MultiVector.");
    // return Teuchos::null; // unreachable
  }

  // non-const version
  static Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> v, const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> cv =
        Teuchos::rcp_const_cast<const Thyra::MultiVectorBase<Scalar>>(v);
    Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> r =
        toXpetra(cv, comm);
    return Teuchos::rcp_const_cast<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(r);
  }

  static bool isTpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
    // check whether we have a Tpetra based Thyra operator
    bool bIsTpetra = false;
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))

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

  static bool isEpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
    // check whether we have an Epetra based Thyra operator
    bool bIsEpetra = false;

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<const Thyra::EpetraLinearOp> epetra_op = Teuchos::rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op, false);
    bIsEpetra                                           = Teuchos::is_null(epetra_op) ? false : true;
#endif

#if 0
    // Check whether it is a blocked operator.
    // If yes, grab the (0,0) block and check the underlying linear algebra
    // Note: we expect that the (0,0) block exists!
    if(bIsEpetra == false) {
      Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> > ThyBlockedOp =
          Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar> >(op,false);
      if(ThyBlockedOp != Teuchos::null) {
        TEUCHOS_TEST_FOR_EXCEPT(ThyBlockedOp->blockExists(0,0)==false);
        Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > b00 =
          ThyBlockedOp->getBlock(0,0);
        TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(b00));
        bIsEpetra = isEpetra(b00);
      }
    }
#endif

    return bIsEpetra;
  }

  static bool isBlockedOperator(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
    // Check whether it is a blocked operator.
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar>> ThyBlockedOp =
        Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar>>(op);
    if (ThyBlockedOp != Teuchos::null) {
      return true;
    }
    return false;
  }

  static Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
#ifdef HAVE_XPETRA_TPETRA
    if (isTpetra(op)) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))

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
#else
      throw Xpetra::Exceptions::RuntimeError("Add TPETRA_INST_INT_LONG_LONG:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    if (isEpetra(op)) {
      Teuchos::RCP<const Epetra_Operator> epetra_op = Thyra::get_Epetra_Operator(*op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_op));
      Teuchos::RCP<const Epetra_RowMatrix> epetra_rowmat = Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(epetra_op, true);
      Teuchos::RCP<const Epetra_CrsMatrix> epetra_crsmat = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(epetra_rowmat, true);
      Teuchos::RCP<Epetra_CrsMatrix> epetra_ncnstcrsmat  = Teuchos::rcp_const_cast<Epetra_CrsMatrix>(epetra_crsmat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_ncnstcrsmat));

      Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> xEpetraCrsMat =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>(epetra_ncnstcrsmat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xEpetraCrsMat));

      Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsMat =
          Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xEpetraCrsMat, true);
      Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsWrap =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpCrsMat));
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpMat =
          Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpCrsWrap, true);
      return xpMat;
    }
#endif
    return Teuchos::null;
  }

  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<Thyra::LinearOpBase<Scalar>>& op) {
#ifdef HAVE_XPETRA_TPETRA
    if (isTpetra(op)) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))

      typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
      Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraOp = TOE::getTpetraOperator(op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraOp));
      Teuchos::RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraRowMat = Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraOp, true);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraCrsMat = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(TpetraRowMat, true);

      Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpetraCrsMat =
          Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(TpetraCrsMat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpetraCrsMat));

      Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsMat =
          Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xTpetraCrsMat, true);
      Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsWrap =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpCrsMat));
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpMat =
          Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpCrsWrap, true);
      return xpMat;
#else
      throw Xpetra::Exceptions::RuntimeError("Add TPETRA_INST_INT_LONG_LONG:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    if (isEpetra(op)) {
      Teuchos::RCP<Epetra_Operator> epetra_op = Thyra::get_Epetra_Operator(*op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_op));
      Teuchos::RCP<Epetra_RowMatrix> epetra_rowmat = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(epetra_op, true);
      Teuchos::RCP<Epetra_CrsMatrix> epetra_crsmat = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(epetra_rowmat, true);

      Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> xEpetraCrsMat =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>(epetra_crsmat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xEpetraCrsMat));

      Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsMat =
          Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xEpetraCrsMat, true);
      Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpCrsWrap =
          Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xpCrsMat));
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpMat =
          Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpCrsWrap, true);
      return xpMat;
    }
#endif
    return Teuchos::null;
  }

  static Teuchos::RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetraOperator(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>& op) {
#ifdef HAVE_XPETRA_TPETRA
    if (isTpetra(op)) {
      typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node> TOE;
      Teuchos::RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TpetraOp = TOE::getConstTpetraOperator(op);

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

  static Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
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

  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<Thyra::DiagonalLinearOpBase<Scalar>>& op) {
    using Teuchos::rcp_const_cast;
    using Teuchos::rcp_dynamic_cast;
    using ThyVSBase = Thyra::SpmdVectorSpaceBase<Scalar>;
    using thyTpV    = Thyra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using tV        = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    RCP<const Thyra::VectorBase<Scalar>> diag = op->getDiag();

    RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xpDiag;
#ifdef HAVE_XPETRA_TPETRA
    if (!rcp_dynamic_cast<const thyTpV>(diag).is_null()) {
      RCP<const tV> tDiag = Thyra::TpetraOperatorVectorExtraction<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getConstTpetraVector(diag);
      if (!tDiag.is_null())
        xpDiag = Xpetra::toXpetra(tDiag);
    }
#endif
#ifdef HAVE_XPETRA_EPETRA
    if (xpDiag.is_null()) {
      RCP<const Epetra_Comm> comm = Thyra::get_Epetra_Comm(*rcp_dynamic_cast<const ThyVSBase>(op->range())->getComm());
      RCP<const Epetra_Map> map   = Thyra::get_Epetra_Map(*(op->range()), comm);
      if (!map.is_null()) {
        RCP<const Epetra_Vector> eDiag                 = Thyra::get_Epetra_Vector(*map, diag);
        RCP<Epetra_Vector> nceDiag                     = rcp_const_cast<Epetra_Vector>(eDiag);
        RCP<Xpetra::EpetraVectorT<int, Node>> xpEpDiag = rcp(new Xpetra::EpetraVectorT<int, Node>(nceDiag));
        xpDiag                                         = rcp_dynamic_cast<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(xpEpDiag, true);
      }
    }
#endif
    TEUCHOS_ASSERT(!xpDiag.is_null());
    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> M = Xpetra::MatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(xpDiag);
    return M;
  }

  static Teuchos::RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  toXpetra(const Teuchos::RCP<const Thyra::DiagonalLinearOpBase<Scalar>>& op) {
    return toXpetra(Teuchos::rcp_const_cast<Thyra::DiagonalLinearOpBase<Scalar>>(op));
  }

  static Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
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
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))
      Teuchos::RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>> tpetraMap = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
      if (tpetraMap == Teuchos::null)
        throw Exceptions::BadCast("Xpetra::ThyraUtils::toThyra: Cast from Xpetra::Map to Xpetra::TpetraMap failed");
      RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> tpMap                         = tpetraMap->getTpetra_Map();
      RCP<Thyra::TpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>> thyraTpetraMap = Thyra::tpetraVectorSpace<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpMap);
      thyraMap                                                                                = thyraTpetraMap;
#else
      throw Xpetra::Exceptions::RuntimeError("Problem DDD. Add TPETRA_INST_INT_LONG_LONG:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    if (map->lib() == Xpetra::UseEpetra) {
      Teuchos::RCP<const Xpetra::EpetraMapT<GlobalOrdinal, Node>> epetraMap = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GlobalOrdinal, Node>>(map);
      if (epetraMap == Teuchos::null)
        throw Exceptions::BadCast("Xpetra::ThyraUtils::toThyra: Cast from Xpetra::Map to Xpetra::EpetraMap failed");
      RCP<const Thyra::VectorSpaceBase<Scalar>> thyraEpetraMap = Thyra::create_VectorSpace(epetraMap->getEpetra_MapRCP());
      thyraMap                                                 = thyraEpetraMap;
    }
#endif

    return thyraMap;
  }

  static Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
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
      auto thyEpMap = Thyra::create_VectorSpace(Teuchos::rcp_dynamic_cast<const EpetraMapT<GlobalOrdinal, Node>>(vec->getMap())->getEpetra_MapRCP());
      auto epMV     = Teuchos::rcp_dynamic_cast<const EpetraMultiVectorT<GlobalOrdinal, Node>>(vec)->getEpetra_MultiVector();
      auto thyMV    = Thyra::create_MultiVector(epMV, thyEpMap);
      return thyMV;
    }
#endif

    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "MultiVector cannot be converted to Thyra.");
  }

  static Teuchos::RCP<const Thyra::VectorBase<Scalar>>
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
      auto thyEpMap = Thyra::create_VectorSpace(Teuchos::rcp_dynamic_cast<const EpetraMapT<GlobalOrdinal, Node>>(vec->getMap())->getEpetra_MapRCP());
      auto epVec    = rcp(Teuchos::rcp_dynamic_cast<const EpetraVectorT<GlobalOrdinal, Node>>(vec)->getEpetra_Vector(), false);
      auto thyVec   = Thyra::create_Vector(epVec, thyEpMap);
      return thyVec;
    }
#endif

    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError, "Vector cannot be converted to Thyra.");
  }

  static void updateThyra(Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> source, Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mapExtractor, const Teuchos::RCP<Thyra::MultiVectorBase<Scalar>>& target) {
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

  static Teuchos::RCP<const Thyra::LinearOpBase<Scalar>>
  toThyra(const Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat) {
    // create a Thyra operator from Xpetra::CrsMatrix
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> thyraOp = Teuchos::null;

#ifdef HAVE_XPETRA_TPETRA
    Teuchos::RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetraMat = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat);
    if (tpetraMat != Teuchos::null) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))

      Teuchos::RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpCrsMat = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat, true);
      Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpCrsMat        = xTpCrsMat->getTpetra_CrsMatrix();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));

      Teuchos::RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpRowMat  = Teuchos::rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpCrsMat, true);
      Teuchos::RCP<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpOperator = Teuchos::rcp_dynamic_cast<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpRowMat, true);

      thyraOp = Thyra::createConstLinearOp(tpOperator);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraOp));
#else
      throw Xpetra::Exceptions::RuntimeError("Add TPETRA_INST_INT_LONG_LONG:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> epetraMat = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(mat);
    if (epetraMat != Teuchos::null) {
      Teuchos::RCP<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> xEpCrsMat = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(mat, true);
      Teuchos::RCP<const Epetra_CrsMatrix> epCrsMat                               = xEpCrsMat->getEpetra_CrsMatrix();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epCrsMat));

      Teuchos::RCP<const Thyra::EpetraLinearOp> thyraEpOp = Thyra::epetraLinearOp(epCrsMat, "op");
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraEpOp));
      thyraOp = thyraEpOp;
    }
#endif
    return thyraOp;
  }

  static Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
  toThyra(const Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& mat) {
    // create a Thyra operator from Xpetra::CrsMatrix
    Teuchos::RCP<Thyra::LinearOpBase<Scalar>> thyraOp = Teuchos::null;

#ifdef HAVE_XPETRA_TPETRA
    Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat);
    if (tpetraMat != Teuchos::null) {
#if ((defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))

      Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> xTpCrsMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(mat, true);
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpCrsMat        = xTpCrsMat->getTpetra_CrsMatrixNonConst();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));

      Teuchos::RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpRowMat  = Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpCrsMat, true);
      Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpOperator = Teuchos::rcp_dynamic_cast<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tpRowMat, true);

      thyraOp = Thyra::createLinearOp(tpOperator);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraOp));
#else
      throw Xpetra::Exceptions::RuntimeError("Add TPETRA_INST_INT_LONG_LONG:BOOL=ON in your configuration.");
#endif
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> epetraMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(mat);
    if (epetraMat != Teuchos::null) {
      Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>> xEpCrsMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal, Node>>(mat, true);
      Teuchos::RCP<Epetra_CrsMatrix> epCrsMat                               = xEpCrsMat->getEpetra_CrsMatrixNonConst();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epCrsMat));

      Teuchos::RCP<Thyra::EpetraLinearOp> thyraEpOp = Thyra::nonconstEpetraLinearOp(epCrsMat, "op");
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraEpOp));
      thyraOp = thyraEpOp;
    }
#endif
    return thyraOp;
  }

  static Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
  toThyra(const Teuchos::RCP<Xpetra::BlockedCrsMatrix<double, int, long long, EpetraNode>>& mat);

};      // specialization on SC=double, LO=GO=int and NO=EpetraNode
#endif  // XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES

#endif  // HAVE_XPETRA_EPETRA

}  // end namespace Xpetra

#define XPETRA_THYRAUTILS_SHORT
#endif  // HAVE_XPETRA_THYRA

#endif  // XPETRA_THYRAUTILS_HPP
