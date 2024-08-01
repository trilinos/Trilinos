// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Xpetra_UnitTestHelpers.hpp>
#include <Teuchos_Comm.hpp>

#include "RTOpPack_ROpNorm1.hpp"

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_Vector.hpp"

#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"

#include "Xpetra_ThyraUtils.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_TpetraVector.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_EpetraVector.hpp"
#endif

#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"

namespace {

bool testMpi = true;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
}

Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm() {
  Teuchos::RCP<const Teuchos::Comm<int> > ret;
  if (testMpi) {
    ret = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
  } else {
    ret = rcp(new Teuchos::SerialComm<int>());
  }
  return ret;
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(Map, Create, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Scalar scalar_type;
  typedef Xpetra::Map<LO, GO, Node> map_type;
  typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
  typedef Xpetra::ThyraUtils<Scalar, LO, GO, Node> th_utils_type;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib

  Teuchos::Array<GlobalOrdinal> gids;
  gids.push_back(Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30);
  gids.push_back(Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30 + 10);
  gids.push_back(Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30 + 20);

  // create an Xpetra map
  Teuchos::RCP<const map_type> map = map_factory_type::Build(mylib,
                                                             3 * Teuchos::as<GlobalOrdinal>(comm->getSize()),
                                                             gids.view(0, 3),
                                                             0,
                                                             comm);

  // create Thyra vector space out of Xpetra Map
  Teuchos::RCP<const Thyra::VectorSpaceBase<scalar_type> > thMap = th_utils_type::toThyra(map);
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<Teuchos::Ordinal>(map->getGlobalNumElements()) != thMap->dim(), std::logic_error, "Global dimension of Xpetra map and Thyra VectorSpaceBase are different.");
  Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<scalar_type> > thSpmdMap = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<scalar_type> >(thMap);
  TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMap == Teuchos::null, std::logic_error, "Cannot cast VectorSpaceBase to SpmdVectorSpaceBase.");
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<Teuchos::Ordinal>(map->getLocalNumElements()) != thSpmdMap->localSubDim(), std::logic_error, "Local dimension of Xpetra map and Thyra VectorSpaceBase on one (or more) processor(s) are different.");

  Teuchos::RCP<const map_type> map2 = th_utils_type::toXpetra(thMap, comm);
  TEST_EQUALITY(map2->getGlobalNumElements(), 3 * Teuchos::as<Xpetra::global_size_t>(comm->getSize()));
  TEST_EQUALITY(map2->getLocalNumElements(), 3);
  TEST_EQUALITY(map2->getMinGlobalIndex(), Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30);
  TEST_EQUALITY(map2->getMaxGlobalIndex(), Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30 + 20);
  TEST_EQUALITY(map2->getGlobalElement(0), Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30);
  TEST_EQUALITY(map2->getGlobalElement(1), Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30 + 10);
  TEST_EQUALITY(map2->getGlobalElement(2), Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30 + 20);
  TEST_EQUALITY(map2->getLocalElement(Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30), 0);
  TEST_EQUALITY(map2->getLocalElement(Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30 + 10), 1);
  TEST_EQUALITY(map2->getLocalElement(Teuchos::as<GlobalOrdinal>(comm->getRank()) * 30 + 20), 2);
  TEST_EQUALITY(map2->isNodeGlobalElement(Teuchos::as<GlobalOrdinal>(1)), false);
  TEST_EQUALITY(map2->isNodeLocalElement(Teuchos::as<LocalOrdinal>(1)), true);
  TEST_EQUALITY(map2->isContiguous(), false);
  TEST_EQUALITY(map2->getIndexBase(), 0);
  TEST_EQUALITY(map2->isSameAs(*map), true);
  TEST_EQUALITY(map->isSameAs(*map2), true);
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, Create, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Scalar scalar_type;
  typedef Xpetra::Map<LO, GO, Node> map_type;
  typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
  typedef Xpetra::ThyraUtils<Scalar, LO, GO, Node> th_utils_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib

  // create an Xpetra map
  const LO numInd                  = 63;
  Teuchos::RCP<const map_type> map = map_factory_type::Build(mylib, numInd, 0, comm);

  // create Thyra vector space out of Xpetra Map
  Teuchos::RCP<const Thyra::VectorSpaceBase<scalar_type> > thMap = th_utils_type::toThyra(map);
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<Teuchos::Ordinal>(map->getGlobalNumElements()) != thMap->dim(), std::logic_error, "Global dimension of Xpetra map and Thyra VectorSpaceBase are different.");
  Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<scalar_type> > thSpmdMap = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<scalar_type> >(thMap);
  TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMap == Teuchos::null, std::logic_error, "Cannot cast VectorSpaceBase to SpmdVectorSpaceBase.");
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<Teuchos::Ordinal>(map->getLocalNumElements()) != thSpmdMap->localSubDim(), std::logic_error, "Local dimension of Xpetra map and Thyra VectorSpaceBase on one (or more) processor(s) are different.");

  // create Thyra MultiVector
  Teuchos::RCP<Thyra::MultiVectorBase<scalar_type> > thMVec         = Thyra::createMembers(thMap, 2);
  Teuchos::RCP<Thyra::SpmdMultiVectorBase<scalar_type> > thSpmdMVec = Teuchos::rcp_dynamic_cast<Thyra::SpmdMultiVectorBase<scalar_type> >(thMVec);
  TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMVec == Teuchos::null, std::logic_error, "Cannot cast MultiVectorBase to SpmdMultiVectorBase.");

  // fill multivector with some data
  const LocalOrdinal localOffset = (thSpmdMap != Teuchos::null ? thSpmdMap->localOffset() : 0);
  const LocalOrdinal localSubDim = (thSpmdMap != Teuchos::null ? thSpmdMap->localSubDim() : thMap->dim());
  Teuchos::RCP<Thyra::DetachedMultiVectorView<scalar_type> > thyData =
      Teuchos::rcp(new Thyra::DetachedMultiVectorView<scalar_type>(*thSpmdMVec, Teuchos::Range1D(localOffset, localOffset + localSubDim - 1)));

  // loop over all vectors in multivector
  for (LocalOrdinal j = 0; j < thSpmdMVec->domain()->dim(); ++j) {
    // loop over all local rows
    for (LocalOrdinal i = 0; i < localSubDim; ++i) {
      (*thyData)(i, 0) = 1;
      (*thyData)(i, 1) = 2;
    }
  }

  // calculate and check 1-norm of Thyra MultiVector
  RTOpPack::ROpNorm1<scalar_type> op;
  const LocalOrdinal numVec = thSpmdMVec->domain()->dim();
  Teuchos::Array<Teuchos::RCP<RTOpPack::ReductTarget> > rcp_op_targs(numVec);
  Teuchos::Array<Teuchos::Ptr<RTOpPack::ReductTarget> > op_targs(numVec);
  for (LocalOrdinal kc = 0; kc < numVec; ++kc) {
    rcp_op_targs[kc] = op.reduct_obj_create();
    op_targs[kc]     = rcp_op_targs[kc].ptr();
  }
  ::Thyra::applyOp<scalar_type>(op, Teuchos::tuple(Teuchos::ptrInArg(*(thMVec.ptr()))),
                                Teuchos::ArrayView<Teuchos::Ptr<Thyra::MultiVectorBase<scalar_type> > >(Teuchos::null),
                                op_targs);
  TEST_EQUALITY(op(*op_targs[0]), numInd);
  TEST_EQUALITY(op(*op_targs[1]), 2 * numInd);

  // create Xpetra multivector from Thyra multi vector
  Teuchos::RCP<mv_type> xpMVec = th_utils_type::toXpetra(thMVec, comm);
  TEUCHOS_TEST_FOR_EXCEPTION(xpMVec == Teuchos::null, std::logic_error, "Failed to convert Thyra::MultiVector to Xpetra::MultiVector.");
  TEST_EQUALITY(Teuchos::as<Teuchos::Ordinal>(xpMVec->getNumVectors()), numVec);
  TEST_EQUALITY(Teuchos::as<Teuchos::Ordinal>(xpMVec->getLocalLength()), localSubDim);
  TEST_EQUALITY(Teuchos::as<LO>(xpMVec->getGlobalLength()), numInd);

  std::vector<typename Teuchos::ScalarTraits<scalar_type>::magnitudeType> norms(numVec, STS::magnitude(STS::zero()));
  Teuchos::ArrayView<typename Teuchos::ScalarTraits<scalar_type>::magnitudeType> normsView(norms);
  xpMVec->norm1(normsView);
  TEST_EQUALITY(normsView[0], numInd);
  TEST_EQUALITY(normsView[1], 2 * numInd);
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, CreateProductMV, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Scalar scalar_type;
  typedef Xpetra::Map<LO, GO, Node> map_type;
  typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
  typedef Xpetra::ThyraUtils<Scalar, LO, GO, Node> th_utils_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib

  // create an Xpetra map
  const LO numA                     = 63;
  const LO numB                     = 24;
  Teuchos::RCP<const map_type> mapA = map_factory_type::Build(mylib, numA, 0, comm);
  Teuchos::RCP<const map_type> mapB = map_factory_type::Build(mylib, numB, 0, comm);

  // create Thyra vector space out of Xpetra Map
  Teuchos::RCP<const Thyra::VectorSpaceBase<scalar_type> > thMapA = th_utils_type::toThyra(mapA);
  Teuchos::RCP<const Thyra::VectorSpaceBase<scalar_type> > thMapB = th_utils_type::toThyra(mapB);

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<Teuchos::Ordinal>(mapA->getGlobalNumElements()) != thMapA->dim(), std::logic_error, "Global dimension of Xpetra map and Thyra VectorSpaceBase are different.");
  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<Teuchos::Ordinal>(mapB->getGlobalNumElements()) != thMapB->dim(), std::logic_error, "Global dimension of Xpetra map and Thyra VectorSpaceBase are different.");

  Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<scalar_type> > thSpmdMapA = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<scalar_type> >(thMapA);
  Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<scalar_type> > thSpmdMapB = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<scalar_type> >(thMapB);
  TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMapA == Teuchos::null, std::logic_error, "Cannot cast VectorSpaceBase to SpmdVectorSpaceBase.");
  TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMapB == Teuchos::null, std::logic_error, "Cannot cast VectorSpaceBase to SpmdVectorSpaceBase.");

  // create Thyra MultiVector
  Teuchos::RCP<Thyra::MultiVectorBase<scalar_type> > thMVecA         = Thyra::createMembers(thMapA, 1);
  Teuchos::RCP<Thyra::MultiVectorBase<scalar_type> > thMVecB         = Thyra::createMembers(thMapB, 1);
  Teuchos::RCP<Thyra::SpmdMultiVectorBase<scalar_type> > thSpmdMVecA = Teuchos::rcp_dynamic_cast<Thyra::SpmdMultiVectorBase<scalar_type> >(thMVecA);
  Teuchos::RCP<Thyra::SpmdMultiVectorBase<scalar_type> > thSpmdMVecB = Teuchos::rcp_dynamic_cast<Thyra::SpmdMultiVectorBase<scalar_type> >(thMVecB);
  TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMVecA == Teuchos::null, std::logic_error, "Cannot cast MultiVectorBase to SpmdMultiVectorBase.");
  TEUCHOS_TEST_FOR_EXCEPTION(thSpmdMVecB == Teuchos::null, std::logic_error, "Cannot cast MultiVectorBase to SpmdMultiVectorBase.");

  // fill multivector A with some data
  const LocalOrdinal localOffsetA = (thSpmdMapA != Teuchos::null ? thSpmdMapA->localOffset() : 0);
  const LocalOrdinal localSubDimA = (thSpmdMapA != Teuchos::null ? thSpmdMapA->localSubDim() : thMapA->dim());
  Teuchos::RCP<Thyra::DetachedMultiVectorView<scalar_type> > thyDataA =
      Teuchos::rcp(new Thyra::DetachedMultiVectorView<scalar_type>(*thSpmdMVecA, Teuchos::Range1D(localOffsetA, localOffsetA + localSubDimA - 1)));

  // loop over all vectors in multivector
  for (LocalOrdinal j = 0; j < thSpmdMVecA->domain()->dim(); ++j) {
    // loop over all local rows
    for (LocalOrdinal i = 0; i < localSubDimA; ++i) {
      (*thyDataA)(i, j) = 1;
    }
  }

  // fill multivector B with some data
  const LocalOrdinal localOffsetB = (thSpmdMapB != Teuchos::null ? thSpmdMapB->localOffset() : 0);
  const LocalOrdinal localSubDimB = (thSpmdMapB != Teuchos::null ? thSpmdMapB->localSubDim() : thMapB->dim());
  Teuchos::RCP<Thyra::DetachedMultiVectorView<scalar_type> > thyDataB =
      Teuchos::rcp(new Thyra::DetachedMultiVectorView<scalar_type>(*thSpmdMVecB, Teuchos::Range1D(localOffsetB, localOffsetB + localSubDimB - 1)));

  // loop over all vectors in multivector
  for (LocalOrdinal j = 0; j < thSpmdMVecB->domain()->dim(); ++j) {
    // loop over all local rows
    for (LocalOrdinal i = 0; i < localSubDimB; ++i) {
      (*thyDataB)(i, j) = 2;
    }
  }

  Teuchos::RCP<Thyra::DefaultProductVectorSpace<scalar_type> > thyProdVecSpace = Thyra::productVectorSpace(Teuchos::tuple(thMapA, thMapB));
  TEUCHOS_TEST_FOR_EXCEPTION(thyProdVecSpace == Teuchos::null, std::logic_error, "Failed to create product vector space.");

  // create product multi vector from multivectors A and B
  Teuchos::RCP<Thyra::DefaultProductMultiVector<scalar_type> > thyProdAB = Thyra::defaultProductMultiVector<scalar_type>(thyProdVecSpace, Teuchos::tuple(thMVecA, thMVecB));
  TEUCHOS_TEST_FOR_EXCEPTION(thyProdAB == Teuchos::null, std::logic_error, "Failed to create product multivector.");
  Teuchos::RCP<Thyra::MultiVectorBase<scalar_type> > thyProdMVec = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<scalar_type> >(thyProdAB);
  TEUCHOS_TEST_FOR_EXCEPTION(thyProdMVec == Teuchos::null, std::logic_error, "Downcast of product multivector to multivector failed.");

  // calculate and check 1-norm of Thyra MultiVector
  RTOpPack::ROpNorm1<scalar_type> op;
  const LocalOrdinal numVec = thSpmdMVecA->domain()->dim();
  Teuchos::Array<Teuchos::RCP<RTOpPack::ReductTarget> > rcp_op_targs(numVec);
  Teuchos::Array<Teuchos::Ptr<RTOpPack::ReductTarget> > op_targs(numVec);
  for (LocalOrdinal kc = 0; kc < numVec; ++kc) {
    rcp_op_targs[kc] = op.reduct_obj_create();
    op_targs[kc]     = rcp_op_targs[kc].ptr();
  }
  ::Thyra::applyOp<scalar_type>(op, Teuchos::tuple(Teuchos::ptrInArg(*(thyProdMVec.ptr()))),
                                Teuchos::ArrayView<Teuchos::Ptr<Thyra::MultiVectorBase<scalar_type> > >(Teuchos::null),
                                op_targs);
  TEST_EQUALITY(op(*op_targs[0]), numA + 2 * numB);

  // create Xpetra multivector from Thyra multi vector
  Teuchos::RCP<mv_type> xpMVec = th_utils_type::toXpetra(thyProdMVec, comm);
  TEUCHOS_TEST_FOR_EXCEPTION(xpMVec == Teuchos::null, std::logic_error, "Downcast of product multivector to multivector failed.");

  std::vector<typename Teuchos::ScalarTraits<scalar_type>::magnitudeType> norms(1, STS::magnitude(STS::zero()));
  Teuchos::ArrayView<typename Teuchos::ScalarTraits<scalar_type>::magnitudeType> normsView(norms);
  xpMVec->norm1(normsView);

  TEST_EQUALITY(normsView[0], numA + 2 * numB);

  // extract sub-blocks from Thyra multivector
  Teuchos::RCP<const Thyra::MultiVectorBase<scalar_type> > thybA = thyProdAB->getMultiVectorBlock(0);
  Teuchos::RCP<const Thyra::MultiVectorBase<scalar_type> > thybB = thyProdAB->getMultiVectorBlock(1);
  Teuchos::RCP<const mv_type> xpbA                               = th_utils_type::toXpetra(thybA, comm);
  Teuchos::RCP<const mv_type> xpbB                               = th_utils_type::toXpetra(thybB, comm);
  TEUCHOS_TEST_FOR_EXCEPTION(xpbA == Teuchos::null, std::logic_error, "Transformation from Thyra::MultiVector to Xpetra::MultiVector failed.");
  TEUCHOS_TEST_FOR_EXCEPTION(xpbB == Teuchos::null, std::logic_error, "Transformation from Thyra::MultiVector to Xpetra::MultiVector failed.");
  TEUCHOS_TEST_FOR_EXCEPTION(xpbA->getMap()->isSameAs(*mapA) == false, std::logic_error, "Map mismatch.");
  TEUCHOS_TEST_FOR_EXCEPTION(xpbB->getMap()->isSameAs(*mapB) == false, std::logic_error, "Map mismatch.");
  TEUCHOS_TEST_FOR_EXCEPTION(xpbA->getMap()->getMinAllGlobalIndex() != 0, std::logic_error, "Map inconsistency.");
  TEUCHOS_TEST_FOR_EXCEPTION(xpbB->getMap()->getMinAllGlobalIndex() != 0, std::logic_error, "Map inconsistency.");
  TEUCHOS_TEST_FOR_EXCEPTION(xpbA->getVector(0)->norm1() != numA, std::logic_error, "SubVector contains wrong entries.");
  TEUCHOS_TEST_FOR_EXCEPTION(xpbB->getVector(0)->norm1() != 2 * numB, std::logic_error, "SubVector contains wrong entries.");
}

//
// INSTANTIATIONS
//
#ifdef HAVE_XPETRA_TPETRA

#define XPETRA_TPETRA_TYPES(S, LO, GO, N)                                    \
  typedef typename Xpetra::TpetraMap<LO, GO, N> M##LO##GO##N;                \
  typedef typename Xpetra::TpetraMultiVector<S, LO, GO, N> MV##S##LO##GO##N; \
  typedef typename Xpetra::TpetraVector<S, LO, GO, N> V##S##LO##GO##N;

#endif

#ifdef HAVE_XPETRA_EPETRA

#define XPETRA_EPETRA_TYPES(S, LO, GO, N)                              \
  typedef typename Xpetra::EpetraMapT<GO, N> M##LO##GO##N;             \
  typedef typename Xpetra::EpetraMultiVectorT<GO, N> MV##S##LO##GO##N; \
  typedef typename Xpetra::EpetraVectorT<GO, N> V##S##LO##GO##N;

#endif

// list of all tests which run both with Epetra and Tpetra
#define XP_THYRAMULTIVECTOR_INSTANT(S, LO, GO, N)                                                                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(Map, Create, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, Create, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, CreateProductMV, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)

#if defined(HAVE_XPETRA_TPETRA)
#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>
TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_THYRAMULTIVECTOR_INSTANT)
#endif

#ifdef HAVE_XPETRA_EPETRA
typedef Xpetra::EpetraNode EpetraNode;

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_TYPES(double, int, int, EpetraNode)
XP_THYRAMULTIVECTOR_INSTANT(double, int, int, EpetraNode)
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
// Thyra has support for Epetra only but not for Epetra64
// EEE(double,int,LongLong,Xpetra::EpetraNode)
// XP_THYRAMULTIVECTOR_INSTANT(double,int,LongLong,Xpetra::EpetraNode)
#endif
#endif  // HAVE_TPETRA_SERIAL

}  // namespace
