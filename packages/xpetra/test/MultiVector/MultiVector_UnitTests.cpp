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
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Range1D.hpp>

#ifdef HAVE_XPETRA_TPETRA
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_TpetraVector.hpp"
#include "Tpetra_Details_Behavior.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_EpetraVector.hpp"
#endif  // HAVE_XPETRA_EPETRA

#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVectorFactory.hpp"  // taw: include MultiVectorFactory before VectorFactory (for BlockedMultiVector definition)
#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_MapExtractor.hpp"

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

// FINISH: add test for MultiVector with a node containing zero local entries
// FINISH: add tests for local MultiVectors

namespace Teuchos {
template <>
ScalarTraits<int>::magnitudeType
relErr(const int &s1, const int &s2) {
  typedef ScalarTraits<int> ST;
  return ST::magnitude(s1 - s2);
}

template <>
ScalarTraits<char>::magnitudeType
relErr(const char &s1, const char &s2) {
  typedef ScalarTraits<char> ST;
  return ST::magnitude(s1 - s2);
}
}  // namespace Teuchos

namespace {
using std::copy;
using std::endl;
using std::ostream_iterator;
using std::string;

using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::arrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::CONJ_TRANS;
using Teuchos::NO_TRANS;
using Teuchos::null;
using Teuchos::OrdinalTraits;
using Teuchos::Range1D;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ScalarTraits;
using Teuchos::SerialDenseMatrix;
using Teuchos::TRANS;
using Teuchos::Tuple;
using Teuchos::tuple;
using Teuchos::VERB_DEFAULT;
using Teuchos::VERB_EXTREME;
using Teuchos::VERB_HIGH;
using Teuchos::VERB_LOW;
using Teuchos::VERB_MEDIUM;
using Teuchos::VERB_NONE;

//   using Tpetra::Map;
//   using Tpetra::MultiVector;
using Xpetra::DefaultPlatform;
using Xpetra::global_size_t;
using Xpetra::GloballyDistributed;

//  using Tpetra::createContigMapWithNode;
//  using Tpetra::createLocalMapWithNode;
#ifdef HAVE_XPETRA_TPETRA
using Xpetra::useTpetra::createContigMapWithNode;
using Xpetra::useTpetra::createLocalMapWithNode;
#endif

bool testMpi         = true;
double errorTolSlack = 1.0e+2;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
  clp.setOption(
      "error-tol-slack", &errorTolSlack,
      "Slack off of machine epsilon used to check test results");
}

RCP<const Comm<int> > getDefaultComm() {
  RCP<const Comm<int> > ret;
  if (testMpi) {
    ret = DefaultPlatform::getDefaultPlatform().getComm();
  } else {
    ret = rcp(new Teuchos::SerialComm<int>());
  }
  return ret;
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, GetVector, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Scalar scalar_type;
  typedef Xpetra::Map<LO, GO, Node> map_type;
  typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> vec_type;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> mv_factory_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  RCP<const Comm<int> > comm = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib

  const size_t numVecs = 11;

  // Create a Map for the MultiVector X, and create X.
  const LO numInd         = 63;
  RCP<const map_type> map = map_factory_type::Build(mylib, numInd, 0, comm);
  RCP<mv_type> X          = mv_factory_type::Build(map, numVecs);
  X->putScalar(STS::zero());

  // Make sure that X has the correct number of columns.
  TEST_EQUALITY(X->getNumVectors(), numVecs);

  // Fill all entries of the j-th column of X with the value j.
  // Use X->getVectorNonConst(j) to get the j-th column of X.
  for (size_t j = 0; j < numVecs; ++j) {
    RCP<vec_type> X_j = X->getVectorNonConst(j);
    X_j->putScalar(as<scalar_type>(j));
  }

  const size_t numRows         = map->getGlobalNumElements();
  const magnitude_type normTol = as<magnitude_type>(numRows) * STM::eps();

  Teuchos::Array<magnitude_type> correctNorms(numVecs);
  Teuchos::Array<magnitude_type> allAtOnceNorms(numVecs);
  Teuchos::Array<magnitude_type> oneAtATimeNorms(numVecs);

  // Calculate what the 2-norm of each column of X should be.
  for (size_t j = 0; j < numVecs; ++j) {
    correctNorms[j] = as<magnitude_type>(j) *
                      STM::squareroot(as<magnitude_type>(numRows));
  }

  // Fill the Array with zeros, and use X->norm2(ArrayView) to
  // compute the 2-norm of each column of X.  Test the results.
  std::fill(allAtOnceNorms.begin(), allAtOnceNorms.end(), STM::zero());
  X->norm2(allAtOnceNorms);
  TEST_COMPARE_FLOATING_ARRAYS(correctNorms(), allAtOnceNorms(), normTol);

  // Use X_j->norm2() to compute the 2-norm of each column of X, and
  // use X->getVector(j) to get X_j.  Test the resulting norms.
  for (size_t j = 0; j < numVecs; ++j) {
    RCP<const vec_type> X_j = X->getVector(j);
    oneAtATimeNorms[j]      = X_j->norm2();
  }
  TEST_COMPARE_FLOATING_ARRAYS(correctNorms(), oneAtATimeNorms(), normTol);

  // Calculate what the 1-norm of each column of X should be.
  for (size_t j = 0; j < numVecs; ++j) {
    correctNorms[j] = as<magnitude_type>(j) * as<magnitude_type>(numRows);
  }

  // Fill the Array with zeros, and use X->norm1(ArrayView) to
  // compute the 1-norm of each column of X.  Test the results.
  std::fill(allAtOnceNorms.begin(), allAtOnceNorms.end(), STM::zero());
  X->norm1(allAtOnceNorms);
  TEST_COMPARE_FLOATING_ARRAYS(correctNorms(), allAtOnceNorms(), normTol);

  // Use X_j->norm1() to compute the 1-norm of each column of X, and
  // use X->getVector(j) to get X_j.  Test the resulting norms.
  for (size_t j = 0; j < numVecs; ++j) {
    RCP<const vec_type> X_j = X->getVector(j);
    oneAtATimeNorms[j]      = X_j->norm1();
  }
  TEST_COMPARE_FLOATING_ARRAYS(correctNorms(), oneAtATimeNorms(), normTol);

  // Calculate what the inf-norm of each column of X should be.
  for (size_t j = 0; j < numVecs; ++j) {
    correctNorms[j] = as<magnitude_type>(j);
  }

  // Fill the Array with zeros, and use X->normInf(ArrayView) to
  // compute the inf-norm of each column of X.  Test the results.
  std::fill(allAtOnceNorms.begin(), allAtOnceNorms.end(), STM::zero());
  X->normInf(allAtOnceNorms);
  TEST_COMPARE_FLOATING_ARRAYS(correctNorms(), allAtOnceNorms(), normTol);

  // Use X_j->normInf() to compute the inf-norm of each column of X,
  // and use X->getVector(j) to get X_j.  Test the resulting norms.
  for (size_t j = 0; j < numVecs; ++j) {
    RCP<const vec_type> X_j = X->getVector(j);
    oneAtATimeNorms[j]      = X_j->normInf();
  }
  TEST_COMPARE_FLOATING_ARRAYS(correctNorms(), oneAtATimeNorms(), normTol);
}

//
// Bug 6115 test: Ensure that Xpetra::Vector::operator= does a deep copy.
//
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(Vector, AssignmentDeepCopies, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Scalar scalar_type;
  typedef Xpetra::Map<LO, GO, Node> map_type;
  typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
  typedef Xpetra::Vector<Scalar, LO, GO, Node> vec_type;
  typedef Xpetra::VectorFactory<Scalar, LO, GO, Node> vec_factory_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  RCP<const Comm<int> > comm = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib

  // Create a Map, which will be the row, domain, and range Map of the matrix A.
  const LO numInd         = 63;
  RCP<const map_type> map = map_factory_type::Build(mylib, numInd, 0, comm);

  RCP<vec_type> v = vec_factory_type::Build(map);
  v->putScalar(STS::one());

  // Remember the norm of v, to make sure that neither apply() call
  // changes it.  Remember both the 2-norm and the 1-norm.
  const magnitude_type v_norm    = v->norm2();
  const magnitude_type v_oneNorm = v->norm1();

  // Test the computed 1-norm of v against what we know it should
  // be.  This test ensures that we aren't just subtracting zero
  // from zero in all the tests below.
  const magnitude_type v_expectedOneNorm =
      as<magnitude_type>(map->getGlobalNumElements());

  TEST_FLOATING_EQUALITY(
      v_expectedOneNorm,
      v_oneNorm,
      STM::squareroot(v_expectedOneNorm) * STM::eps());

  // Keep a copy of v, to test that neither apply() call changes it.
  RCP<vec_type> vcopy = vec_factory_type::Build(map);

  // Xpetra's operator= does a deep copy, like Epetra, but unlike
  // Tpetra (as of early 2014).
  *vcopy = *v;

  // Make sure that vcopy and v have the same norm.  It's OK for the
  // norms to be slightly different, due to nondeterminism in
  // parallel collectives.
  const magnitude_type vcopy_norm    = vcopy->norm2();
  const magnitude_type vcopy_oneNorm = vcopy->norm1();

  const magnitude_type norm_tol =
      static_cast<magnitude_type>(map->getGlobalNumElements()) * STM::eps();

  TEST_FLOATING_EQUALITY(v_norm, vcopy_norm, norm_tol);
  TEST_FLOATING_EQUALITY(v_oneNorm, vcopy_oneNorm, norm_tol);

  // Make sure that if you change vcopy, v doesn't change.  That is,
  // vcopy must be a true deep copy of v.  Do this by setting all
  // entries of the Vector to 2, using putScalar(2).
  {
    vcopy->putScalar(as<scalar_type>(2));
    // Changing all the entries from 1 to 2 should double the
    // 1-norm, but give a little wiggle room for rounding error.
    const magnitude_type two = STM::one() + STM::one();
    // First make sure that the 1-norm of vcopy is as expected.
    TEST_FLOATING_EQUALITY(two * v_oneNorm, vcopy->norm1(), norm_tol);
    // Now make sure that the 1-norm of v has not changed.
    TEST_FLOATING_EQUALITY(v_oneNorm, v->norm1(), norm_tol);

    // Restore vcopy, using v.
    *vcopy = *v;
  }

  // Make sure that if you change vcopy, v doesn't change.  That is,
  // vcopy must be a true deep copy of v.  Do this by setting the
  // first local entry of the Vector to 10000, using
  // getDataNonConst(0).
  {
    Teuchos::ArrayRCP<Scalar> vcopy_data = vcopy->getDataNonConst(0);
    if (map->getLocalNumElements() != 0) {
      vcopy_data[0] += static_cast<magnitude_type>(10000.0);
    }
    // Destroy the view, so that the changes get written back to the Vector.
    vcopy_data = Teuchos::null;

    // Adding 10000 to an entry had better change the 2-norm by at least sqrt(10000) = 100.
    const magnitude_type minChange = static_cast<magnitude_type>(100.0);
    TEUCHOS_TEST_COMPARE(
        STM::magnitude(vcopy_norm - vcopy->norm2()), >, minChange,
        out, success);

    // Restore the original vcopy, by doing a deep copy again.
    // Xpetra's operator= does a deep copy, like Epetra, but unlike
    // Tpetra (as of early 2014).
    *vcopy = *v;

    // Make sure the original copy got restored.
    TEST_FLOATING_EQUALITY(vcopy_norm, vcopy->norm2(), norm_tol);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, AssignmentDeepCopies, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Scalar scalar_type;
  typedef Xpetra::Map<LO, GO, Node> map_type;
  typedef Xpetra::MapFactory<LO, GO, Node> map_factory_type;
  typedef Xpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
  typedef Xpetra::MultiVectorFactory<Scalar, LO, GO, Node> mv_factory_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  RCP<const Comm<int> > comm = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib

  // Create a Map, which will be the row, domain, and range Map of the matrix A.
  const LO numInd         = 63;
  RCP<const map_type> map = map_factory_type::Build(mylib, numInd, 0, comm);
  const size_t numVecs    = 11;

  RCP<mv_type> X = mv_factory_type::Build(map, numVecs);
  Teuchos::Array<magnitude_type> X_correctOneNorms(numVecs);
  const GO numRows = map->getGlobalNumElements();
  for (size_t j = 0; j < numVecs; ++j) {
    X->getVectorNonConst(j)->putScalar(as<scalar_type>(j));
    X_correctOneNorms[j] = as<magnitude_type>(numRows) * as<magnitude_type>(j);
  }

  // Remember the norms of the columns of X, to make sure that
  // MultiVector::operator= really does a deep copy.  Remember both
  // the two-norms and the one-norms.
  Teuchos::Array<magnitude_type> X_twoNorms(numVecs);
  X->norm2(X_twoNorms());
  Teuchos::Array<magnitude_type> X_oneNorms(numVecs);
  X->norm1(X_oneNorms());

  // Test the computed 1-norms of the columns of X against what we
  // know they should be.  This test ensures that we aren't just
  // subtracting zero from zero in all the tests below.
  const magnitude_type normTol =
      STM::squareroot(as<magnitude_type>(numRows)) * STM::eps();
  TEST_COMPARE_FLOATING_ARRAYS(X_correctOneNorms(), X_oneNorms(), normTol);

  // Keep a copy of X, to test that operator= really does a deep copy.
  RCP<mv_type> X_copy = mv_factory_type::Build(map, numVecs);
  *X_copy             = *X;

  // Make sure that the columns of X and X_copy have the same norms.
  Teuchos::Array<magnitude_type> X_copy_twoNorms(numVecs);
  X_copy->norm2(X_copy_twoNorms());
  Teuchos::Array<magnitude_type> X_copy_oneNorms(numVecs);
  X_copy->norm1(X_copy_oneNorms());

  TEST_COMPARE_FLOATING_ARRAYS(X_oneNorms(), X_copy_oneNorms(), normTol);
  TEST_COMPARE_FLOATING_ARRAYS(X_twoNorms(), X_copy_twoNorms(), normTol);

  // Make sure that if you change X_copy, X doesn't change.  That
  // is, X_copy must be a true deep copy of X.  Do this by setting
  // all entries of X_copy to 2, using putScalar(2).
  {
    X_copy->putScalar(as<scalar_type>(2));
    Teuchos::Array<magnitude_type> newOneNorms(numVecs);
    X_copy->norm1(newOneNorms());

    Teuchos::Array<magnitude_type> expectedOneNorms(numVecs);
    for (size_t j = 0; j < numVecs; ++j) {
      expectedOneNorms[j] =
          as<magnitude_type>(2) * as<magnitude_type>(numRows);
    }

    // First make sure that the 1-norms of the columns X_copy are as expected.
    TEST_COMPARE_FLOATING_ARRAYS(newOneNorms(), expectedOneNorms(), normTol);
    // Now make sure that the 1-norms of the columns of X have not changed.
    X->norm1(newOneNorms());
    TEST_COMPARE_FLOATING_ARRAYS(newOneNorms(), X_oneNorms(), normTol);

    // Restore X_copy to be a deep copy of X.
    *X_copy = *X;
  }
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, NonMemberConstructorsEpetra, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_EPETRA
  // typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib

  // create a Map
  const size_t numLocal = 13;
  const size_t numVecs  = 7;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, numLocal, comm);
  // Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(Xpetra::UseEpetra, INVALID, numLocal, 0, comm);
  if (mylib == Xpetra::UseEpetra) {
    RCP<const Xpetra::EpetraMapT<GlobalOrdinal, Node> > emap = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraMapT<GlobalOrdinal, Node> >(map);
    RCP<Epetra_MultiVector> mvec                             = Teuchos::rcp(new Epetra_MultiVector(emap->getEpetra_Map(), numVecs));
    RCP<Epetra_Vector> vec                                   = Teuchos::rcp(new Epetra_Vector(emap->getEpetra_Map()));
    RCP<MV> xmv                                              = Teuchos::rcp_dynamic_cast<MV>(Xpetra::toXpetra<GlobalOrdinal, Node>(mvec));
    // RCP<V>  xv  = Teuchos::rcp_dynamic_cast<V >(Xpetra::toXpetra<GlobalOrdinal,Node>(vec)); // there is no toXpetra for Vectors!
    TEST_EQUALITY(xmv->getNumVectors(), numVecs);
    // TEST_EQUALITY_CONST(xv->getNumVectors(), 1);
  }
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, NonMemberConstructorsTpetra, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib

  // create a Map
  const size_t numLocal = 13;
  const size_t numVecs  = 7;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, numLocal, comm);

#ifdef HAVE_XPETRA_TPETRA
  if (mylib == Xpetra::UseTpetra) {
    RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> > tmap     = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
    RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > mvec = Tpetra::createMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmap->getTpetra_Map(), numVecs);
    RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > vec       = Tpetra::createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmap->getTpetra_Map());
    TEST_EQUALITY(mvec->getNumVectors(), numVecs);
    TEST_EQUALITY_CONST(vec->getNumVectors(), 1);
    RCP<const MV> xmv = Teuchos::rcp_dynamic_cast<const MV>(Xpetra::toXpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mvec));
    RCP<const V> xv   = Teuchos::rcp_dynamic_cast<const V>(Xpetra::toXpetra<Scalar, LocalOrdinal, GlobalOrdinal, Node>(vec));
    TEST_EQUALITY(xmv->getNumVectors(), numVecs);
    TEST_EQUALITY_CONST(xv->getNumVectors(), 1);
  }
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, basic, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  const int numRanks          = comm->getSize();
  const int myRank            = comm->getRank();
  EXTRACT_LIB(comm, M)  // returns mylib
  // create a Map
  const size_t numLocal = 13;
  const size_t numVecs  = 7;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, numLocal, comm);
  MV mvec(map, numVecs, true);
  TEST_EQUALITY(mvec.getNumVectors(), numVecs);
  TEST_EQUALITY(mvec.getLocalLength(), numLocal);
  TEST_EQUALITY(mvec.getGlobalLength(), numRanks * numLocal);

// Norms are not computed by Epetra_IntMultiVector so far
#ifdef HAVE_XPETRA_EPETRA
  if (!std::is_same_v<typename MV::node_type, Xpetra::EpetraNode>)
#endif
  {
    if (!(std::is_same_v<typename MV::scalar_type, int> || std::is_same_v<typename MV::scalar_type, long long int>)) {
      out << "Running the norm tests!" << std::endl;
      // we zeroed it out in the constructor; all norms should be zero
      Array<Magnitude> norms(numVecs), zeros(numVecs);
      std::fill(zeros.begin(), zeros.end(), ScalarTraits<Magnitude>::zero());
      mvec.norm2(norms);
      TEST_COMPARE_ARRAYS(norms, zeros);
      mvec.norm1(norms);
      TEST_COMPARE_ARRAYS(norms, zeros);
      mvec.normInf(norms);
      TEST_COMPARE_ARRAYS(norms, zeros);
    }
  }

  Scalar testValue = 2, sumValue = 3;
  LocalOrdinal testLID  = 7;
  GlobalOrdinal testGID = myRank * numLocal + testLID;
  out << "myRank: " << myRank << ", testGID=" << testGID << std::endl;
  mvec.replaceLocalValue(testLID, 3, testValue);
  mvec.replaceLocalValue(testLID, 4, testValue);
  mvec.sumIntoLocalValue(testLID, 4, sumValue);
  mvec.replaceGlobalValue(testGID, 5, testValue);
  mvec.replaceGlobalValue(testGID, 6, testValue);
  mvec.sumIntoGlobalValue(testGID, 6, sumValue);
  ArrayRCP<const Scalar> replaceLocalData  = mvec.getData(3);
  ArrayRCP<const Scalar> sumIntoLocalData  = mvec.getData(4);
  ArrayRCP<const Scalar> replaceGlobalData = mvec.getData(5);
  ArrayRCP<const Scalar> sumIntoGlobalData = mvec.getData(6);

  if (std::is_same_v<typename MV::scalar_type, int> || std::is_same_v<typename MV::scalar_type, long long int>) {
    TEST_EQUALITY(replaceLocalData[testLID], testValue);
    TEST_EQUALITY(sumIntoLocalData[testLID], testValue + sumValue);
    TEST_EQUALITY(replaceGlobalData[testLID], testValue);
    TEST_EQUALITY(sumIntoGlobalData[testLID], testValue + sumValue);
  } else {
    TEST_FLOATING_EQUALITY(replaceLocalData[testLID], testValue, 1.0e-10);
    TEST_FLOATING_EQUALITY(sumIntoLocalData[testLID], testValue + sumValue, 1.0e-10);
    TEST_FLOATING_EQUALITY(replaceGlobalData[testLID], testValue, 1.0e-10);
    TEST_FLOATING_EQUALITY(sumIntoGlobalData[testLID], testValue + sumValue, 1.0e-10);
  }
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, BadConstNumVecs, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib
  // create a Map
  const size_t numLocal = 13;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, numLocal, comm);
  TEST_THROW(MV mvec(map, 0), std::invalid_argument);
  if (std::numeric_limits<size_t>::is_signed) {
    TEST_THROW(MV mvec(map, INVALID), std::invalid_argument);
  }
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, BadConstLDA, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  // numlocal > LDA
  // ergo, the arrayview doesn't contain enough data to specify the entries
  // also, if bounds checking is enabled, check that bad bounds are caught

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib
  const size_t numLocal = 2;
  const size_t numVecs  = 2;
  // multivector has two vectors, each proc having two values per vector
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, numLocal, comm);

  // we need 4 scalars to specify values on each proc
  Array<Scalar> values(4);
#ifdef HAVE_TPETRA_DEBUG
  // too small an ArrayView (less than 4 values) is met with an exception, if debugging is on
  TEST_THROW(MV mvec(map, values(0, 3), 2, numVecs), std::runtime_error);
  // it could also be too small for the given LDA:
  TEST_THROW(MV mvec(map, values(), 2 + 1, numVecs), std::runtime_error);
  // too small for number of entries in a Vector
  TEST_THROW(V vec(map, values(0, 1)), std::runtime_error);
#endif
  // LDA < numLocal throws an exception anytime
  TEST_THROW(MV mvec(map, values(0, 4), 1, numVecs), std::runtime_error);
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, NonContigView, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  if (ScalarTraits<Scalar>::isOrdinal) return;

  typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
  const Mag tol               = errorTolSlack * errorTolSlack * ScalarTraits<Mag>::eps();  // extra slack on this test; dots() seem to be a little sensitive for single precision types
  const Mag M0                = ScalarTraits<Mag>::zero();
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib
  // create a Map
  const size_t numLocal = 53;  // making this larger reduces the change that A below will have no non-zero entries, i.e., that C = abs(A) is still equal to A (we assume it is not)
  const size_t numVecs  = 7;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, numLocal, comm);

  //
  // we will create a non-contig subview of the vector; un-viewed vectors should not be changed
  Tuple<size_t, 4> inView1 = tuple<size_t>(1, 4, 3, 2);
  Tuple<size_t, 3> exView1 = tuple<size_t>(0, 5, 6);
  Tuple<size_t, 4> inView2 = tuple<size_t>(6, 0, 4, 3);
  Tuple<size_t, 4> exView2 = tuple<size_t>(1, 2, 5, 7);
  const size_t numView     = 4;
  TEUCHOS_TEST_FOR_EXCEPTION(numView != as<size_t>(inView1.size()), std::logic_error, "Someone ruined a test invariant.");
  TEUCHOS_TEST_FOR_EXCEPTION(numView != as<size_t>(inView1.size()), std::logic_error, "Someone ruined a test invariant.");  // TODO: is it a duplication ? an error ?
  TEUCHOS_TEST_FOR_EXCEPTION(numView != as<size_t>(inView2.size()), std::logic_error, "Someone ruined a test invariant.");
  {
    // test dot, all norms, randomize
    MV mvOrig1(map, numVecs), mvOrig2(map, numVecs + 1), mvWeights(map, numVecs);
    mvWeights.randomize();
    RCP<const MV> mvW1         = mvWeights.subView(tuple<size_t>(0));
    RCP<const MV> mvSubWeights = mvWeights.subView(inView1);
    mvOrig1.randomize();
    mvOrig2.randomize();
    //
    Array<Mag> nOrig2(numVecs), nOrig1(numVecs), nOrigI(numVecs), nOrigW(numVecs), nOrigW1(numVecs);
    Array<Scalar> meansOrig(numVecs), dotsOrig(numView);
    mvOrig1.norm1(nOrig1());
    mvOrig1.norm2(nOrig2());
    mvOrig1.normInf(nOrigI());
    mvOrig1.meanValue(meansOrig());
    for (size_t j = 0; j < numView; ++j) {
      RCP<const V> v1 = mvOrig1.getVector(inView1[j]),
                   v2 = mvOrig2.getVector(inView2[j]);
      dotsOrig[j]     = v1->dot(*v2);
    }
    // create the views, compute and test
    RCP<MV> mvView1       = mvOrig1.subViewNonConst(inView1);
    RCP<const MV> mvView2 = mvOrig2.subView(inView2);
    Array<Mag> nView2(numView), nView1(numView), nViewI(numView), nViewW(numView), nViewW1(numView);
    Array<Scalar> meansView(numView), dotsView(numView);
    mvView1->norm1(nView1());
    mvView1->norm2(nView2());
    mvView1->normInf(nViewI());
    mvView1->meanValue(meansView());
    mvView1->dot(*mvView2, dotsView());
    for (size_t j = 0; j < numView; ++j) {
      TEST_FLOATING_EQUALITY(nOrig1[inView1[j]], nView1[j], tol);
      TEST_FLOATING_EQUALITY(nOrig2[inView1[j]], nView2[j], tol);
      TEST_FLOATING_EQUALITY(nOrigI[inView1[j]], nViewI[j], tol);
      TEST_FLOATING_EQUALITY(meansOrig[inView1[j]], meansView[j], tol);
      TEST_FLOATING_EQUALITY(dotsOrig[j], dotsView[j], tol);
    }
    // randomize the view, compute view one-norms, test difference
    mvView2 = Teuchos::null;
    mvView1->randomize();
    Array<Mag> nView1_aft(numView);
    mvView1->norm1(nView1_aft());
    for (size_t j = 0; j < numView; ++j) {
      TEST_INEQUALITY(nView1[j], nView1_aft[j]);
    }
    // release the view, test that viewed columns changed, others didn't
    mvView1 = Teuchos::null;
    Array<Mag> nOrig1_aft(numVecs);
    mvOrig1.norm1(nOrig1_aft());
    for (size_t j = 0; j < as<size_t>(inView1.size()); ++j) {
      TEST_INEQUALITY(nOrig1[inView1[j]], nOrig1_aft[inView1[j]]);
    }
    for (size_t j = 0; j < as<size_t>(exView1.size()); ++j) {
      TEST_FLOATING_EQUALITY(nOrig1[exView1[j]], nOrig1_aft[exView1[j]], tol);
    }
  }
  {
    MV mvOrigA(map, numVecs), mvOrigB(map, numVecs), mvOrigC(map, numVecs + 1);
    mvOrigA.randomize();
    mvOrigB.randomize();
    mvOrigC.randomize();
    Array<Mag> nrmOrigA(numVecs), nrmOrigB(numVecs), nrmOrigC(numVecs + 1);
    mvOrigA.norm2(nrmOrigA());
    mvOrigB.norm2(nrmOrigB());
    mvOrigC.norm2(nrmOrigC());
    RCP<MV> mvViewA = mvOrigA.subViewNonConst(inView1);
    RCP<MV> mvViewB = mvOrigB.subViewNonConst(inView1);
    RCP<MV> mvViewC = mvOrigC.subViewNonConst(inView2);
    // set C = abs(A)
    {
      Array<Scalar> mnA_bef(inView1.size()), mnC_bef(inView1.size()),
          mnA_aft(inView1.size()), mnC_aft(inView1.size());
      mvViewA->meanValue(mnA_bef());
      mvViewC->meanValue(mnC_bef());
      mvViewC->abs(*mvViewA);
      mvViewA->meanValue(mnA_aft());
      mvViewC->meanValue(mnC_aft());
      for (size_t j = 0; j < as<size_t>(inView1.size()); ++j) {
        TEST_FLOATING_EQUALITY(mnA_bef[j], mnA_aft[j], tol);
        TEST_INEQUALITY(mnC_bef[j], mnC_aft[j]);
      }
    }
    // then set A = B = C
    // good excuse for some double views
    // use full views of C and B for this, check means before and after
    // to make sure that only A and B change.
    {
      Array<Scalar> A_bef(inView1.size()), B_bef(inView1.size()), C_bef(inView2.size());
      mvViewA->meanValue(A_bef());
      mvViewB->meanValue(B_bef());
      mvViewC->meanValue(C_bef());
      RCP<MV> doubleViewA       = mvViewA->subViewNonConst(Range1D(0, inView1.size() - 1));
      RCP<MV> doubleViewB       = mvViewB->subViewNonConst(Range1D(0, inView1.size() - 1));
      RCP<const MV> doubleViewC = mvViewC->subView(Range1D(0, inView1.size() - 1));
      (*doubleViewA) = (*doubleViewB) = (*doubleViewC);
      doubleViewA                     = Teuchos::null;
      doubleViewB                     = Teuchos::null;
      doubleViewC                     = Teuchos::null;
      Array<Scalar> A_aft(inView1.size()), B_aft(inView1.size()), C_aft(inView2.size());
      mvViewA->meanValue(A_aft());
      mvViewB->meanValue(B_aft());
      mvViewC->meanValue(C_aft());
      for (size_t j = 0; j < as<size_t>(inView1.size()); ++j) {
        TEST_FLOATING_EQUALITY(C_bef[j], C_aft[j], tol);
        TEST_FLOATING_EQUALITY(C_bef[j], B_aft[j], tol);
        TEST_FLOATING_EQUALITY(C_bef[j], A_aft[j], tol);
        TEST_INEQUALITY(A_bef[j], A_aft[j]);
        TEST_INEQUALITY(B_bef[j], B_aft[j]);
      }
    }
    {
      TEUCHOS_TEST_FOR_EXCEPTION(inView1.size() != 4, std::logic_error, "Someone ruined a test invariant.");
      Tuple<size_t, 4> reorder = tuple<size_t>(3, 1, 0, 2);
      RCP<MV> dvA              = mvViewA->subViewNonConst(reorder);
      RCP<MV> dvB              = mvViewB->subViewNonConst(reorder);
      RCP<MV> dvC              = mvViewC->subViewNonConst(reorder);
      // C == B == A
      //   C *= 2                ->  C == 2*A == 2*B            scale(alpha)
      dvC->scale(as<Scalar>(2));
      //   A = -C + 2*A          ->  C == 2*B, A == 0           update(alpha,mv,beta)
      dvA->update(as<Scalar>(-1), *dvC, as<Scalar>(2));
      //   C = 2*A + 2*B - .5*C ->   C == B, A == 0,            update(alpha,mv,beta,mv,gamma)
      dvC->update(as<Scalar>(2), *dvA, as<Scalar>(2), *dvB, as<Scalar>(-.5));
      //   B = 0.5              ->   B = 0.5, A == 0,           putScalar(alpha)
      dvB->putScalar(as<Scalar>(0.5));
      //   C.recip(B)           ->   C = 2, B == 0.5, A == 0,   reciprocal(mv)
      dvC->reciprocal(*dvB);
      //   B = C/2              ->   A == 0, B == 1, C == 2
      dvB->scale(as<Mag>(0.5), *dvC);
      dvA = Teuchos::null;
      dvB = Teuchos::null;
      dvC = Teuchos::null;
      Array<Mag> nrmA(4), nrmB(4), nrmC(4);
      mvViewA->norm1(nrmA());  // norm1(0)   = 0
      mvViewB->norm1(nrmB());  // norm1(1.0) = N
      mvViewC->norm1(nrmC());  // norm1(2.0) = 2 * N
      const Mag OneN = as<Mag>(mvViewA->getGlobalLength());
      const Mag TwoN = OneN + OneN;
      for (size_t j = 0; j < 4; ++j) {
        TEST_FLOATING_EQUALITY(nrmA[j], M0, tol);
        TEST_FLOATING_EQUALITY(nrmB[j], OneN, tol);
        TEST_FLOATING_EQUALITY(nrmC[j], TwoN, tol);
      }
    }
    // done with these views; clear them, ensure that only the viewed
    // vectors changed in the original multivectors
    mvViewA = Teuchos::null;
    mvViewB = Teuchos::null;
    mvViewC = Teuchos::null;
    Array<Mag> nrmOrigA_aft(numVecs), nrmOrigB_aft(numVecs), nrmOrigC_aft(numVecs + 1);
    mvOrigA.norm2(nrmOrigA_aft());
    mvOrigB.norm2(nrmOrigB_aft());
    mvOrigC.norm2(nrmOrigC_aft());
    for (size_t j = 0; j < as<size_t>(inView1.size()); ++j) {
      TEST_INEQUALITY(nrmOrigA[inView1[j]], nrmOrigA_aft[inView1[j]]);
      TEST_INEQUALITY(nrmOrigB[inView1[j]], nrmOrigB_aft[inView1[j]]);
      TEST_INEQUALITY(nrmOrigC[inView2[j]], nrmOrigC_aft[inView2[j]]);
    }
    for (size_t j = 0; j < as<size_t>(exView1.size()); ++j) {
      TEST_FLOATING_EQUALITY(nrmOrigA[exView1[j]], nrmOrigA_aft[exView1[j]], tol);
      TEST_FLOATING_EQUALITY(nrmOrigB[exView1[j]], nrmOrigB_aft[exView1[j]], tol);
    }
    for (size_t j = 0; j < as<size_t>(exView1.size()); ++j) {
      TEST_FLOATING_EQUALITY(nrmOrigC[exView2[j]], nrmOrigC_aft[exView2[j]], tol);
    }
  }
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, Describable, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  const LocalOrdinal INVALID = OrdinalTraits<LocalOrdinal>::invalid();
  RCP<const Comm<int> > comm = getDefaultComm();
  const int myImageID        = comm->getRank();
  EXTRACT_LIB(comm, M)  // returns mylib
  // create Map
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, 3, comm);

  // test labeling
  const string lbl("mvecA");
  MV mvecA(map, 2);
  string desc1 = mvecA.description();
  if (myImageID == 0) out << desc1 << endl;
#ifdef XPETRA_NOT_IMPLEMENTED
  mvecA.setObjectLabel(lbl);
#endif
  string desc2 = mvecA.description();
  if (myImageID == 0) out << desc2 << endl;
  if (myImageID == 0) {
#ifdef XPETRA_NOT_IMPLEMENTED
    TEST_EQUALITY(mvecA.getObjectLabel(), lbl);
#endif
  }
  // test describing at different verbosity levels
  if (myImageID == 0) out << "Describing with verbosity VERB_DEFAULT..." << endl;
  mvecA.describe(out);
  comm->barrier();
  comm->barrier();
  if (myImageID == 0) out << "Describing with verbosity VERB_NONE..." << endl;
  mvecA.describe(out, VERB_NONE);
  comm->barrier();
  comm->barrier();
  if (myImageID == 0) out << "Describing with verbosity VERB_LOW..." << endl;
  mvecA.describe(out, VERB_LOW);
  comm->barrier();
  comm->barrier();
  if (myImageID == 0) out << "Describing with verbosity VERB_MEDIUM..." << endl;
  mvecA.describe(out, VERB_MEDIUM);
  comm->barrier();
  comm->barrier();
  if (myImageID == 0) out << "Describing with verbosity VERB_HIGH..." << endl;
  mvecA.describe(out, VERB_HIGH);
  comm->barrier();
  comm->barrier();
  if (myImageID == 0) out << "Describing with verbosity VERB_EXTREME..." << endl;
  mvecA.describe(out, VERB_EXTREME);
  comm->barrier();
  comm->barrier();
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, BadMultiply, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
#ifdef XPETRA_NOT_IMPLEMENTED
  const Scalar S1 = ScalarTraits<Scalar>::one(),
               S0 = ScalarTraits<Scalar>::zero();
#endif
  // case 1: C(local) = A^X(local) * B^X(local)  : four of these
  {
    // create local Maps
    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map3l = createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(3, comm),
                                                               map2l = createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(2, comm);
    MV mvecA(map3l, 2),
        mvecB(map2l, 3),
        mvecD(map2l, 2);
    // failures, 8 combinations:
    // [NTC],[NTC]: A,B don't match
    // [NTC],[NTC]: C doesn't match A,B
#ifdef XPETRA_NOT_IMPLEMENTED
    TEST_THROW(mvecD.multiply(NO_TRANS, NO_TRANS, S1, mvecA, mvecA, S0), std::runtime_error);      // 2x2: 3x2 x 3x2
    TEST_THROW(mvecD.multiply(NO_TRANS, CONJ_TRANS, S1, mvecA, mvecB, S0), std::runtime_error);    // 2x2: 3x2 x 3x2
    TEST_THROW(mvecD.multiply(CONJ_TRANS, NO_TRANS, S1, mvecB, mvecA, S0), std::runtime_error);    // 2x2: 3x2 x 3x2
    TEST_THROW(mvecD.multiply(CONJ_TRANS, CONJ_TRANS, S1, mvecB, mvecB, S0), std::runtime_error);  // 2x2: 3x2 x 3x2
    TEST_THROW(mvecD.multiply(NO_TRANS, NO_TRANS, S1, mvecA, mvecB, S0), std::runtime_error);      // 2x2: 3x2 x 2x3
    TEST_THROW(mvecD.multiply(NO_TRANS, CONJ_TRANS, S1, mvecA, mvecA, S0), std::runtime_error);    // 2x2: 3x2 x 2x3
    TEST_THROW(mvecD.multiply(CONJ_TRANS, NO_TRANS, S1, mvecB, mvecB, S0), std::runtime_error);    // 2x2: 3x2 x 2x3
    TEST_THROW(mvecD.multiply(CONJ_TRANS, CONJ_TRANS, S1, mvecB, mvecA, S0), std::runtime_error);  // 2x2: 3x2 x 2x3
#endif
  }
  // case 2: C(local) = A^T(distr) * B  (distr)  : one of these
  {
    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map3n = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 3, comm),
                                                               map2n = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 2, comm);
    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map2l = createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(2, comm),
                                                               map3l = createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(3, comm);
    MV mv3nx2(map3n, 2),
        mv2nx2(map2n, 2),
        mv2lx2(map2l, 2),
        mv2lx3(map2l, 3),
        mv3lx2(map3l, 2),
        mv3lx3(map3l, 3);
#ifdef XPETRA_NOT_IMPLEMENTED
    // non-matching input lengths
    TEST_THROW(mv2lx2.multiply(CONJ_TRANS, NO_TRANS, S1, mv3nx2, mv2nx2, S0), std::runtime_error);  // (2 x 3n) x (2n x 2) not compat
    TEST_THROW(mv2lx2.multiply(CONJ_TRANS, NO_TRANS, S1, mv2nx2, mv3nx2, S0), std::runtime_error);  // (2 x 2n) x (3n x 2) not compat
    // non-matching output size
    TEST_THROW(mv3lx3.multiply(CONJ_TRANS, NO_TRANS, S1, mv3nx2, mv3nx2, S0), std::runtime_error);  // (2 x 3n) x (3n x 2) doesn't fit 3x3
    TEST_THROW(mv3lx2.multiply(CONJ_TRANS, NO_TRANS, S1, mv3nx2, mv3nx2, S0), std::runtime_error);  // (2 x 3n) x (3n x 2) doesn't fit 3x2
    TEST_THROW(mv2lx3.multiply(CONJ_TRANS, NO_TRANS, S1, mv3nx2, mv3nx2, S0), std::runtime_error);  // (2 x 3n) x (3n x 2) doesn't fit 2x3
#endif
  }
  // case 3: C(distr) = A  (distr) * B^X(local)  : two of these
  {
    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map3n = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 3, comm),
                                                               map2n = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 2, comm);
    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map2l = createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(2, comm),
                                                               map3l = createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(3, comm);
    MV mv3nx2(map3n, 2),
        mv2nx2(map2n, 2),
        mv2x3(map2l, 3),
        mv3x2(map3l, 2);
#ifdef XPETRA_NOT_IMPLEMENTED
    // non-matching input lengths
    TEST_THROW(mv3nx2.multiply(NO_TRANS, CONJ_TRANS, S1, mv3nx2, mv2x3, S0), std::runtime_error);  // (3n x 2) x (3 x 2) (trans) not compat
    TEST_THROW(mv3nx2.multiply(NO_TRANS, NO_TRANS, S1, mv3nx2, mv3x2, S0), std::runtime_error);    // (3n x 2) x (3 x 2) (nontrans) not compat
    // non-matching output sizes
    TEST_THROW(mv3nx2.multiply(NO_TRANS, CONJ_TRANS, S1, mv3nx2, mv3x2, S0), std::runtime_error);  // (3n x 2) x (2 x 3) doesn't fit 3nx2
    TEST_THROW(mv3nx2.multiply(NO_TRANS, NO_TRANS, S1, mv3nx2, mv2x3, S0), std::runtime_error);    // (3n x 2) x (2 x 3) doesn't fit 3nx2
#endif
  }
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, Multiply, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  using Teuchos::View;
  // typedef typename ScalarTraits<Scalar>::magnitudeType Mag;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
#ifdef XPETRA_NOT_IMPLEMENTED
  const int numImages = comm->getSize();
#endif
  // create a Map
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map3n = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 3, comm),
                                                             map2n = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 2, comm);
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > lmap3 = createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(3, comm),
                                                             lmap2 = createLocalMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(2, comm);
#ifdef XPETRA_NOT_IMPLEMENTED
  const Scalar S1 = ScalarTraits<Scalar>::one(),
               S0 = ScalarTraits<Scalar>::zero();
  const Mag M0    = ScalarTraits<Mag>::zero();
#endif
  // case 1: C(local) = A^X(local) * B^X(local)  : four of these
  // deterministic input/output
  {
    MV mv3x2l(lmap3, 2),
        mv2x3l(lmap2, 3),
        mv2x2l(lmap2, 2),
        mv3x3l(lmap3, 3);
    // fill multivectors with ones
    mv3x2l.putScalar(ScalarTraits<Scalar>::one());
    mv2x3l.putScalar(ScalarTraits<Scalar>::one());
    // fill expected answers Array
    Teuchos::Array<Scalar> check2(4, 3);  // each entry (of four) is the product [1 1 1]*[1 1 1]' = 3
    Teuchos::Array<Scalar> check3(9, 2);  // each entry (of nine) is the product [1 1]*[1 1]' = 2
    // test
    ArrayRCP<const Scalar> tmpView;
#ifdef XPETRA_NOT_IMPLEMENTED
    mv3x3l.multiply(NO_TRANS, NO_TRANS, S1, mv3x2l, mv2x3l, S0);
    tmpView = mv3x3l.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView(0, 9), check3, M0);
    mv2x2l.multiply(NO_TRANS, CONJ_TRANS, S1, mv2x3l, mv2x3l, S0);
    tmpView = mv2x2l.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView(0, 4), check2, M0);
    mv2x2l.multiply(CONJ_TRANS, NO_TRANS, S1, mv3x2l, mv3x2l, S0);
    tmpView = mv2x2l.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView(0, 4), check2, M0);
    mv3x3l.multiply(CONJ_TRANS, CONJ_TRANS, S1, mv2x3l, mv3x2l, S0);
    tmpView = mv3x3l.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView(0, 9), check3, M0);
#endif
  }
  // case 1: C(local) = A^X(local) * B^X(local)  : four of these
  // random input/output
  {
#ifdef XPETRA_NOT_IMPLEMENTED
    Array<Scalar> tmvCopy1(6), tmvCopy2(6);
    ArrayView<Scalar> sdmView(Teuchos::null);
    MV tmv3x2(lmap3, 2),
        tmv2x3(lmap2, 3),
        tmv2x2(lmap2, 2),
        tmv3x3(lmap3, 3);
    // fill multivectors with random, get copy of contents
    tmv3x2.randomize();
    tmv3x2.get1dCopy(tmvCopy1(), 3);
    tmv2x3.randomize();
    tmv2x3.get1dCopy(tmvCopy2(), 2);
    // point SerialDenseMatrices at copies
    SerialDenseMatrix<int, Scalar> sdm3x2(View, tmvCopy1.getRawPtr(), 3, 3, 2);
    SerialDenseMatrix<int, Scalar> sdm2x3(View, tmvCopy2.getRawPtr(), 2, 2, 3);
    // space for answers
    SerialDenseMatrix<int, Scalar> sdm2x2(2, 2), sdm3x3(3, 3);
    // test: perform local Tpetra::MultiVector multiply and Teuchos::SerialDenseMatrix multiply, then check that answers are equivalent
    ArrayRCP<const Scalar> tmpView;

    {
      tmv3x3.multiply(NO_TRANS, NO_TRANS, S1, tmv3x2, tmv2x3, S0);
      sdm3x3.multiply(NO_TRANS, NO_TRANS, S1, sdm3x2, sdm2x3, S0);
      tmpView = tmv3x3.get1dView();
      sdmView = arrayView(sdm3x3.values(), sdm3x3.numRows() * sdm3x3.numCols());
      TEST_COMPARE_FLOATING_ARRAYS(tmpView, sdmView, ScalarTraits<Mag>::eps() * 10.);
    }
    {
      tmv2x2.multiply(NO_TRANS, CONJ_TRANS, S1, tmv2x3, tmv2x3, S0);
      sdm2x2.multiply(NO_TRANS, CONJ_TRANS, S1, sdm2x3, sdm2x3, S0);
      tmpView = tmv2x2.get1dView();
      sdmView = arrayView(sdm2x2.values(), sdm2x2.numRows() * sdm2x2.numCols());
      TEST_COMPARE_FLOATING_ARRAYS(tmpView, sdmView, ScalarTraits<Mag>::eps() * 10.);
    }
    {
      tmv2x2.multiply(CONJ_TRANS, NO_TRANS, S1, tmv3x2, tmv3x2, S0);
      sdm2x2.multiply(CONJ_TRANS, NO_TRANS, S1, sdm3x2, sdm3x2, S0);
      tmpView = tmv2x2.get1dView();
      sdmView = arrayView(sdm2x2.values(), sdm2x2.numRows() * sdm2x2.numCols());
      TEST_COMPARE_FLOATING_ARRAYS(tmpView, sdmView, ScalarTraits<Mag>::eps() * 10.);
    }
    {
      tmv3x3.multiply(CONJ_TRANS, CONJ_TRANS, S1, tmv2x3, tmv3x2, S0);
      sdm3x3.multiply(CONJ_TRANS, CONJ_TRANS, S1, sdm2x3, sdm3x2, S0);
      tmpView = tmv3x3.get1dView();
      sdmView = arrayView(sdm3x3.values(), sdm3x3.numRows() * sdm3x3.numCols());
      TEST_COMPARE_FLOATING_ARRAYS(tmpView, sdmView, ScalarTraits<Mag>::eps() * 10.);
    }
#endif
  }
#ifdef XPETRA_NOT_IMPLEMENTED
  // case 2: C(local) = A^T(distr) * B  (distr)  : one of these
  {
    MV mv3nx2(map3n, 2),
        mv3nx3(map3n, 3),
        // locals
        mv2x2(lmap2, 2),
        mv2x3(lmap2, 3),
        mv3x2(lmap3, 2),
        mv3x3(lmap3, 3);
    // fill multivectors with ones
    mv3nx3.putScalar(ScalarTraits<Scalar>::one());
    mv3nx2.putScalar(ScalarTraits<Scalar>::one());
    // fill expected answers Array
    ArrayRCP<const Scalar> tmpView;
    Teuchos::Array<Scalar> check(9, 3 * numImages);
    // test
    mv2x2.multiply(CONJ_TRANS, NO_TRANS, S1, mv3nx2, mv3nx2, S0);
    tmpView = mv2x2.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView, check(0, tmpView.size()), M0);
    mv2x3.multiply(CONJ_TRANS, NO_TRANS, S1, mv3nx2, mv3nx3, S0);
    tmpView = mv2x3.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView, check(0, tmpView.size()), M0);
    mv3x2.multiply(CONJ_TRANS, NO_TRANS, S1, mv3nx3, mv3nx2, S0);
    tmpView = mv3x2.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView, check(0, tmpView.size()), M0);
    mv3x3.multiply(CONJ_TRANS, NO_TRANS, S1, mv3nx3, mv3nx3, S0);
    tmpView = mv3x3.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView, check(0, tmpView.size()), M0);
  }
  // case 3: C(distr) = A  (distr) * B^X(local)  : two of these
  {
    MV mv3nx2(map3n, 2),
        mv3nx3(map3n, 3),
        // locals
        mv2x3(lmap2, 3);
    // fill multivectors with ones
    mv2x3.putScalar(S1);
    // fill expected answers Array
    ArrayRCP<const Scalar> tmpView;
    Teuchos::Array<Scalar> check2(9, 2), check3(6, 3);
    // test
    mv3nx3.putScalar(S1);
    mv3nx2.putScalar(S1);
    mv3nx3.multiply(NO_TRANS, NO_TRANS, S1, mv3nx2, mv2x3, S0);
    tmpView = mv3nx3.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView, check2, M0);
    mv3nx3.putScalar(S1);
    mv3nx2.putScalar(S1);
    mv3nx2.multiply(NO_TRANS, CONJ_TRANS, S1, mv3nx3, mv2x3, S0);
    tmpView = mv3nx2.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView, check3, M0);
  }
#endif
#endif  // HAVE_XPETRA_TPETRA
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, ElementWiseMultiply, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  using Teuchos::View;
  typedef typename ScalarTraits<Scalar>::magnitudeType Mag;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  EXTRACT_LIB(comm, M)  // returns mylib
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map3n =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, 3, comm);

  const Mag M0    = ScalarTraits<Mag>::zero();
  const Scalar S1 = ScalarTraits<Scalar>::one();
  const Scalar S0 = ScalarTraits<Scalar>::zero();
  {
    // case 1: C = S1*A@B ('@' denotes element-wise multiplication)
    // C has 2 vectors, A has 1 vector, B has 2 vectors.
    // A and B will be filled with 1s, so C should get filled with 1s.
    V A(map3n, 1);
    MV B(map3n, 2),
        C(map3n, 2);
    // fill multivectors with ones
    A.putScalar(ScalarTraits<Scalar>::one());
    B.putScalar(ScalarTraits<Scalar>::one());
    // fill expected answers Array
    Teuchos::Array<Scalar> check2(6, 1);  // each entry (of six) is 1
    // test
    ArrayRCP<const Scalar> tmpView;
    C.elementWiseMultiply(S1, A, B, S0);
    tmpView = C.get1dView();
    TEST_COMPARE_FLOATING_ARRAYS(tmpView(0, 6), check2, M0);
  }
#endif  // HAVE_XPETRA_TPETRA
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, BadConstAA, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // constructor takes ArrayView<ArrayView<Scalar> A, NumVectors
  // A.size() == NumVectors
  // A[i].size() >= MyLength

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();

  EXTRACT_LIB(comm, M)  // returns mylib

  // create a Map
  // multivector has two vectors, each proc having two values per vector
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map2 =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, 2, comm);
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map3 =
      Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(mylib, INVALID, 3, comm);

  // we need 4 scalars to specify values on each proc
  Array<Scalar> values(4);
  Array<ArrayView<const Scalar> > arrOfarr(2, ArrayView<const Scalar>(Teuchos::null));
  Array<ArrayView<const Scalar> > emptyArr;
  arrOfarr[0] = values(0, 2);
  arrOfarr[1] = values(2, 2);
  // arrOfarr.size() == 0
  TEST_THROW(MV mvec(map2, emptyArr(), 0), std::runtime_error);
#ifdef HAVE_TPETRA_DEBUG
  // individual ArrayViews could be too small
  TEST_THROW(MV mvec(map3, arrOfarr(), 2), std::runtime_error);
#endif
}

#if 0  // not compiling in Epetra only mode. not sure why
  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL( MultiVector, BadDot, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal , Node )
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    EXTRACT_LIB(comm,M); // returns mylib
    // create a Map
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map1 =
        Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(mylib, INVALID,1,comm);
    RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map2 =
        Xpetra::UnitTestHelpers::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(mylib, INVALID,2,comm);
    {
      MV mv12(map1,1),
         mv21(map2,1),
         mv22(map2,2);
      Array<Scalar> dots(2);
#ifdef HAVE_XPETRA_TPETRA
      if(mylib==Xpetra::UseTpetra) {
        // incompatible maps
        TEST_THROW(mv12.dot(mv21,dots()),std::runtime_error);
        // incompatible numvecs
        TEST_THROW(mv22.dot(mv21,dots()),std::runtime_error);
        // too small output array
#ifdef TEUCHOS_DEBUG
        TEST_THROW(mv22.dot(mv22,dots(0,1)),std::runtime_error);
#endif
      }
#endif
    }
    {
      V v1(map1),
        v2(map2);
#ifdef HAVE_XPETRA_TPETRA
      if (mylib == Xpetra::UseTpetra) {
        // incompatible maps
        TEST_THROW(v1.dot(v2),std::runtime_error);
        TEST_THROW(v2.dot(v1),std::runtime_error);
        // wrong size output array through MultiVector interface
        Array<Scalar> dots(2);
#ifdef TEUCHOS_DEBUG
        TEST_THROW(v1.dot(v2,dots()),std::runtime_error);
        TEST_THROW(v2.dot(v1,dots()),std::runtime_error);
#endif
#endif
      }
    }
  }
#endif  // if 0 // not compiling in Epetra only mode

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, OrthoDot, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef typename ScalarTraits<Scalar>::magnitudeType Mag;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  // TODO FAILED: const Scalar S0 = ScalarTraits<Scalar>::zero();
  const Mag M0 = ScalarTraits<Mag>::zero();

  RCP<const Comm<int> > comm = getDefaultComm();
  const int numImages        = comm->getSize();

  const size_t numLocal                                          = 2;
  const size_t numVectors                                        = 3;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  const bool zeroOut                                             = true;
  MV mvec1(map, numVectors, zeroOut),
      mvec2(map, numVectors, zeroOut);
  Array<Scalar> dots1(numVectors), dots2(numVectors), zeros(numVectors);
  Array<Mag> norms1(numVectors), norms2(numVectors), ans(numVectors);
  std::fill(zeros.begin(), zeros.end(), ScalarTraits<Scalar>::zero());
  // these should be numerically orthogonal even in finite arithmetic, because both are zero. 1-norms are zero.
  mvec1.dot(mvec2, dots1());
  mvec2.dot(mvec1, dots2());
  TEST_COMPARE_FLOATING_ARRAYS(dots2, zeros, M0);
  TEST_COMPARE_FLOATING_ARRAYS(dots1, zeros, M0);
  // TODO FAILED: TEST_EQUALITY_CONST( mvec1.getVector(0)->dot(*mvec2.getVector(0)), S0);
  mvec1.norm1(norms1());
  mvec2.norm1(norms2());
  std::fill(ans.begin(), ans.end(), M0);
  TEST_COMPARE_FLOATING_ARRAYS(norms1, ans, M0);
  TEST_COMPARE_FLOATING_ARRAYS(norms1, ans, M0);
  // replace local entries s.t.
  // mvec1 = [1 1]  and  mvec2 = [0 0]
  //         [0 0]               [1 1]
  // still numerically orthogonal even in finite arithmetic. norms are numImages.
  for (size_t j = 0; j < numVectors; ++j) {
    mvec1.replaceLocalValue(0, j, ScalarTraits<Scalar>::one());
    mvec2.replaceGlobalValue(map->getGlobalElement(1), j, ScalarTraits<Scalar>::one());
  }
  mvec1.dot(mvec2, dots1());
  mvec2.dot(mvec1, dots2());
  TEST_COMPARE_FLOATING_ARRAYS(dots2, zeros, M0);
  TEST_COMPARE_FLOATING_ARRAYS(dots1, zeros, M0);
  // TODO FAILED: TEST_EQUALITY_CONST( mvec1.getVector(0)->dot(*mvec2.getVector(0)), S0);
  mvec1.norm1(norms1());
  mvec2.norm1(norms2());
  std::fill(ans.begin(), ans.end(), as<Mag>(numImages));
  TEST_COMPARE_FLOATING_ARRAYS(norms1, ans, M0);
  TEST_COMPARE_FLOATING_ARRAYS(norms2, ans, M0);
  // sum into local entries s.t.
  // mvec1 = [1 1]  and  mvec2 = [-1 -1]
  //         [1 1]               [ 1  1]
  // still numerically orthogonal even in finite arithmetic. norms are 2*numImages.
  for (size_t j = 0; j < numVectors; ++j) {
    mvec1.sumIntoLocalValue(1, j, ScalarTraits<Scalar>::one());
    mvec2.sumIntoGlobalValue(map->getGlobalElement(0), j, -ScalarTraits<Scalar>::one());
  }
  mvec1.dot(mvec2, dots1());
  mvec2.dot(mvec1, dots2());
  TEST_COMPARE_FLOATING_ARRAYS(dots2, zeros, M0);
  TEST_COMPARE_FLOATING_ARRAYS(dots1, zeros, M0);
  // TODO FAILED: TEST_EQUALITY_CONST( mvec1.getVector(0)->dot(*mvec2.getVector(0)), S0);
  mvec1.norm1(norms1());
  mvec2.norm1(norms2());
  std::fill(ans.begin(), ans.end(), as<Mag>(2 * numImages));
  TEST_COMPARE_FLOATING_ARRAYS(norms1, ans, M0);
  TEST_COMPARE_FLOATING_ARRAYS(norms2, ans, M0);
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, CopyView, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef typename ScalarTraits<Scalar>::magnitudeType Mag;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const Scalar S0             = ScalarTraits<Scalar>::zero();
  const Mag M0                = ScalarTraits<Mag>::zero();
  const Mag tol               = errorTolSlack * ScalarTraits<Mag>::eps();
  RCP<const Comm<int> > comm  = getDefaultComm();
  // create a Map
  const size_t numLocal                                          = 7;
  const size_t numVectors                                        = 13;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  MV A(map, numVectors, false);
  {
    A.randomize();
    TEUCHOS_TEST_FOR_EXCEPT(numVectors != 13);
    Range1D inds1(8, 12);
    // get a subview and a subcopy of certain vectors of A
    // check that the norms are the same
    // change the view, delete it, verify that the copy doesn't change but that A does
    A.randomize();
    Array<Mag> A_bef(numVectors),
        A_aft(numVectors),
        Av_bef(inds1.size()),
        Av_aft(inds1.size()),
        Ac_bef(inds1.size()),
        Ac_aft(inds1.size());
    A.norm2(A_bef());
    // get view and its norms
#ifdef XPETRA_NOT_IMPLEMENTED
    RCP<MV> Av = A.subViewNonConst(inds1);
    Av->norm2(Av_bef());
    // get copy and its norms
    RCP<MV> Ac = A.subCopy(inds1);
    Ac->norm2(Ac_bef());
    // set view to zero
    Av->putScalar(ScalarTraits<Scalar>::zero());
    // get norms of view
    Av->norm2(Av_aft());
    // free the view, copying data back to A
    Av = Teuchos::null;
    // get norms of A and copy
    Ac->norm2(Ac_aft());
    A.norm2(A_aft());
    // norms of copy and view before should match norms of A
    for (size_t i = 0; i < as<size_t>(inds1.size()); ++i) {
      TEST_FLOATING_EQUALITY(A_bef[inds1.lbound() + i], Ac_bef[i], tol);
    }
    TEST_COMPARE_FLOATING_ARRAYS(Ac_bef, Av_bef, tol);
    // norms of copy (before and after) should match
    TEST_COMPARE_FLOATING_ARRAYS(Ac_bef, Ac_aft, tol);
    // norms of view after should be zero, as should corresponding A norms
    for (size_t i = 0; i < as<size_t>(inds1.size()); ++i) {
      TEST_EQUALITY_CONST(Av_aft[i], M0);
      TEST_EQUALITY_CONST(A_aft[inds1.lbound() + i], M0);
    }
#endif
  }
  {
    A.randomize();
    TEUCHOS_TEST_FOR_EXCEPT(numVectors != 13);
    Tuple<size_t, 5> inds = tuple<size_t>(0, 5, 6, 7, 12);
    // get a subview and a subcopy of certain vectors of A
    // check that the norms are the same
    // change the view, delete it, verify that the copy doesn't change but that A does
    Array<Mag> A_bef(numVectors),
        A_aft(numVectors),
        Av_bef(inds.size()),
        Av_aft(inds.size()),
        Ac_bef(inds.size()),
        Ac_aft(inds.size());
    A.norm2(A_bef());
#ifdef XPETRA_NOT_IMPLEMENTED
    // get view and its norms
    RCP<MV> Av = A.subViewNonConst(inds);
    Av->norm2(Av_bef());
    // get copy and its norms
    RCP<MV> Ac = A.subCopy(inds);
    Ac->norm2(Ac_bef());
    // set view to zero
    Av->putScalar(ScalarTraits<Scalar>::zero());
    // get norms of view
    Av->norm2(Av_aft());
    // free the view, copying data back to A
    Av = Teuchos::null;
    // get norms of A and copy
    Ac->norm2(Ac_aft());
#endif
    A.norm2(A_aft());
    // norms of copy and view before should match norms of A
    for (size_t i = 0; i < as<size_t>(inds.size()); ++i) {
      // TODO FAILED: TEST_FLOATING_EQUALITY( A_bef[inds[i]], Ac_bef[i], tol );
    }
    TEST_COMPARE_FLOATING_ARRAYS(Ac_bef, Av_bef, tol);
    // norms of copy (before and after) should match
    TEST_COMPARE_FLOATING_ARRAYS(Ac_bef, Ac_aft, tol);
    // norms of view after should be zero, as should corresponding A norms
    for (size_t i = 0; i < as<size_t>(inds.size()); ++i) {
      // TODO FAILED: TEST_EQUALITY_CONST( Av_aft[i], M0 );
      // TODO FAILED: TEST_EQUALITY_CONST( A_aft[inds[i]], M0 );
    }
  }
  {
    A.randomize();
    Array<Mag> Anorms(numVectors);
    A.norm2(Anorms());
    TEUCHOS_TEST_FOR_EXCEPT(numVectors != 13);
    for (size_t vc = 0; vc < 2; ++vc) {
      // vc == 0 -> view
      // vc == 1 -> copy
#ifdef XPETRA_NOT_IMPLEMENTED
      for (size_t t = 0; t < 4; ++t) {
        //  t |   outer   |   inner
        // ---|-----------|-----------
        //  0 | ArrayView | ArrayView
        //  1 |  Range1D  | ArrayView
        //  2 | ArrayView |  Range1D
        //  3 |  Range1D  |  Range1D
        //
        // outer grabs 5-9
        // inner grabs 1-3 of those, corresponding to 6-8
        RCP<const MV> sub1, sub2;
        if ((t & 1) == 0) {
          Tuple<size_t, 5> inds = tuple<size_t>(5, 6, 7, 8, 9);
          if (vc == 0)
            sub1 = A.subView(inds);
          else
            sub1 = A.subCopy(inds);
        } else {
          Range1D inds(5, 9);
          if (vc == 0)
            sub1 = A.subView(inds);
          else
            sub1 = A.subCopy(inds);
        }
        TEST_EQUALITY_CONST(sub1->getNumVectors(), 5);
        if ((t & 2) == 0) {
          Tuple<size_t, 3> inds = tuple<size_t>(1, 2, 3);
          if (vc == 0)
            sub2 = sub1->subView(inds);
          else
            sub2 = sub1->subCopy(inds);
        } else {
          Range1D inds(1, 3);
          if (vc == 0)
            sub2 = sub1->subView(inds);
          else
            sub2 = sub1->subCopy(inds);
        }
        TEST_EQUALITY_CONST(sub2->getNumVectors(), 3);
        Array<Mag> subnorms(3);
        sub2->norm2(subnorms());
        TEST_COMPARE_FLOATING_ARRAYS(Anorms(6, 3), subnorms(), tol);
      }
#endif
    }
  }
  {
    A.randomize();
    {
      // check that 1dView and 1dCopy have the same values
      ArrayRCP<const Scalar> view;
      Array<Scalar> copy(numLocal * numVectors);
      view = A.get1dView();
      A.get1dCopy(copy(), numLocal);
      TEST_COMPARE_FLOATING_ARRAYS(view, copy, M0);
    }
    {
      // check that 1dView and 1dCopy have the same values
      ArrayRCP<Scalar> view;
      Array<Scalar> copy(numLocal * numVectors);
      view = A.get1dViewNonConst();
      A.get1dCopy(copy(), numLocal);
      TEST_COMPARE_FLOATING_ARRAYS(view, copy, M0);
      // clear view, ensure that A is zero
      std::fill(view.begin(), view.end(), S0);
      view = Teuchos::null;
      Array<Mag> norms(numVectors), zeros(numVectors, M0);
      A.norm2(norms());
      TEST_COMPARE_FLOATING_ARRAYS(norms, zeros, M0);
    }
    A.randomize();
    {
      // check that 1dView and 1dCopy have the same values
      ArrayRCP<ArrayRCP<const Scalar> > views;
      Array<Scalar> copyspace(numLocal * numVectors);
      Array<ArrayView<Scalar> > copies(numVectors);
      for (size_t j = 0; j < numVectors; ++j) {
        copies[j] = copyspace(numLocal * j, numLocal);
      }
      views = A.get2dView();
      A.get2dCopy(copies());
      for (size_t j = 0; j < numVectors; ++j) {
        TEST_COMPARE_FLOATING_ARRAYS(views[j], copies[j], M0);
      }
    }
    {
      // check that 1dView and 1dCopy have the same values
      Array<Scalar> copyspace(numLocal * numVectors);
      Array<ArrayView<Scalar> > copies(numVectors);
      for (size_t j = 0; j < numVectors; ++j) {
        copies[j] = copyspace(numLocal * j, numLocal);
      }
      ArrayRCP<ArrayRCP<Scalar> > views = A.get2dViewNonConst();
      A.get2dCopy(copies());
      for (size_t j = 0; j < numVectors; ++j) {
        TEST_COMPARE_FLOATING_ARRAYS(views[j], copies[j], M0);
      }
      // clear view, ensure that A is zero
      for (size_t j = 0; j < numVectors; ++j) {
        std::fill(views[j].begin(), views[j].end(), S0);
      }
      views = Teuchos::null;
      Array<Mag> norms(numVectors), zeros(numVectors, M0);
      A.norm2(norms());
      TEST_COMPARE_FLOATING_ARRAYS(norms, zeros, M0);
    }
  }
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, OffsetView, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  // typedef typename ScalarTraits<Scalar>::magnitudeType Mag;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  // TODO const Scalar S0 = ScalarTraits<Scalar>::zero();
  // TODO const Mag M0 = ScalarTraits<Mag>::zero();
  // TODO const Mag tol = errorTolSlack * ScalarTraits<Mag>::eps();
  RCP<const Comm<int> > comm = getDefaultComm();
  // create a Map
  const size_t numLocal1  = 3;
  const size_t numLocal2  = 4;
  const size_t numLocal   = numLocal1 + numLocal2;
  const size_t numVectors = 6;
  Array<size_t> even(tuple<size_t>(1, 3, 5));
  Array<size_t> odd(tuple<size_t>(0, 2, 4));
  TEUCHOS_TEST_FOR_EXCEPTION(even.size() != odd.size(), std::logic_error, "Test setup assumption violated.");
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > fullMap = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map1    = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal1, comm);
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map2    = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal2, comm);
  RCP<MV> A                                                          = rcp(new MV(fullMap, numVectors, false));
#ifdef XPETRA_NOT_IMPLEMENTED
  {
    // contig source multivector
    RCP<MV> A1 = A->offsetViewNonConst(map1, 0);
    RCP<MV> A2 = A->offsetViewNonConst(map2, numLocal1);
    TEST_EQUALITY(A1->getLocalLength(), numLocal1);
    TEST_EQUALITY(A2->getLocalLength(), numLocal2);
    TEST_EQUALITY(A1->getNumVectors(), numVectors);
    TEST_EQUALITY(A2->getNumVectors(), numVectors);
    Array<Mag> A_befr(numVectors),
        A1_befr(numVectors),
        A2_befr(numVectors),
        A_aft1(numVectors),
        A1_aft1(numVectors),
        A2_aft1(numVectors),
        A_aft2(numVectors),
        A1_aft2(numVectors),
        A2_aft2(numVectors);
    // compute norms of A, A1 and A2
    A->randomize();
    A->norm2(A_befr());
    A1->norm2(A1_befr());
    A2->norm2(A2_befr());
    // set A1 = zeros, compute norms of A, A1 and A2
    A1->putScalar(S0);
    A->norm2(A_aft1());
    A1->norm2(A1_aft1());
    A2->norm2(A2_aft1());
    // set A2 = zeros, compute norms of A, A1 and A2
    A2->putScalar(S0);
    A->norm2(A_aft2());
    A1->norm2(A1_aft2());
    A2->norm2(A2_aft2());
    // change to A1 should not affect A2
    // change to A2 should not affect A1
    // change to A1 or A2 should change A
    // A should be zero after setting A1 to zero and A2 to zero
    for (size_t i = 0; i < numVectors; ++i) {
      TEST_EQUALITY_CONST(A_aft1[i] < A_befr[i] + tol, true);  // shrunk as A1 = 0
      TEST_EQUALITY_CONST(A_aft2[i] < A_aft1[i] + tol, true);  // shurnk as A2 = 0
      TEST_EQUALITY_CONST(A_aft2[i], M0);                      // ... to zero
      TEST_EQUALITY_CONST(A1_aft1[i], M0);                     // was set to zero
      TEST_EQUALITY_CONST(A1_aft2[i], M0);                     // should not have been changed
      TEST_FLOATING_EQUALITY(A2_befr[i], A2_aft1[i], tol);     // should not have been changed
      TEST_EQUALITY_CONST(A2_aft2[i], M0);                     // was set to zero
    }
  }
#endif
#ifdef XPETRA_NOT_IMPLEMENTED
  {
    // non-contig source multivector
    RCP<MV> A1e = A->subViewNonConst(even)->offsetViewNonConst(map1, 0);
    RCP<MV> A2e = A->subViewNonConst(even)->offsetViewNonConst(map2, numLocal1);
    RCP<MV> A1o = A->subViewNonConst(odd)->offsetViewNonConst(map1, 0);
    RCP<MV> A2o = A->subViewNonConst(odd)->offsetViewNonConst(map2, numLocal1);
    TEST_EQUALITY(A1e->getLocalLength(), numLocal1);
    TEST_EQUALITY(A1o->getLocalLength(), numLocal1);
    TEST_EQUALITY(A2e->getLocalLength(), numLocal2);
    TEST_EQUALITY(A2o->getLocalLength(), numLocal2);
    const size_t numSubVecs = (size_t)even.size();
    TEST_EQUALITY(A1e->getNumVectors(), numSubVecs);
    TEST_EQUALITY(A2e->getNumVectors(), numSubVecs);
    TEST_EQUALITY(A1o->getNumVectors(), numSubVecs);
    TEST_EQUALITY(A2o->getNumVectors(), numSubVecs);
    A->randomize();
    Array<Mag> b1(numSubVecs), b2(numSubVecs), b3(numSubVecs), bw(numVectors);  // before putScalar(): unchanged 1, 2, 3; whole
    Array<Mag> a1(numSubVecs), a2(numSubVecs), a3(numSubVecs), aw(numVectors);  // after putScalar(): ...
    Array<Mag> changed(numSubVecs), zeros(numSubVecs, M0);
    for (int i = 0; i < 4; ++i) {
      ArrayView<RCP<MV> > allMVs;  // (changed,three unchanged)
      switch (i) {
        case 0:
          allMVs = tuple<RCP<MV> >(A1e, A2e, A1o, A2o);
          break;
        case 1:
          allMVs = tuple<RCP<MV> >(A2e, A1o, A2o, A1e);
          break;
        case 2:
          allMVs = tuple<RCP<MV> >(A1o, A2o, A1e, A2e);
          break;
        case 3:
          allMVs = tuple<RCP<MV> >(A2o, A1e, A2e, A1o);
          break;
      }
      allMVs[1]->norm2(b1());
      allMVs[2]->norm2(b2());
      allMVs[3]->norm2(b3());
      A->norm2(bw());
      allMVs[0]->putScalar(S0);
      allMVs[0]->norm2(changed());
      allMVs[1]->norm2(a1());
      allMVs[2]->norm2(a2());
      allMVs[3]->norm2(a3());
      A->norm2(aw());
      TEST_COMPARE_FLOATING_ARRAYS(b1, a1, tol);
      TEST_COMPARE_FLOATING_ARRAYS(b2, a2, tol);
      TEST_COMPARE_FLOATING_ARRAYS(b3, a3, tol);
      TEST_COMPARE_ARRAYS(changed(), zeros());
      for (size_t i = 0; i < numVectors; ++i) {
        TEST_EQUALITY_CONST(aw[i] < bw[i] + tol, true);  // shrunk
      }
    }
  }
#endif
#ifdef XPETRA_NOT_IMPLEMENTED
  {
    RCP<const MV> A1 = A->offsetView(map1, 0);
    RCP<const MV> A2 = A->offsetView(map2, numLocal1);
    TEST_EQUALITY(A1->getLocalLength(), numLocal1);
    TEST_EQUALITY(A2->getLocalLength(), numLocal2);
    TEST_EQUALITY(A1->getNumVectors(), numVectors);
    TEST_EQUALITY(A2->getNumVectors(), numVectors);
    Array<Mag> A_bef(numVectors),
        A1_bef(numVectors),
        A2_bef(numVectors),
        A_aft(numVectors),
        A1_aft(numVectors),
        A2_aft(numVectors);
    // compute norms of A, A1 and A2
    A->randomize();
    A->norm2(A_bef());
    A1->norm2(A1_bef());
    A2->norm2(A2_bef());
    A->putScalar(S0);
    A->norm2(A_aft());
    A1->norm2(A1_aft());
    A2->norm2(A2_aft());
    for (size_t i = 0; i < numVectors; ++i) {
      TEST_EQUALITY_CONST(A_bef[i] < A1_bef[i] + A2_bef[i] + tol, true);
      TEST_EQUALITY_CONST(A_aft[i], S0);
      TEST_EQUALITY_CONST(A1_aft[i], S0);
      TEST_EQUALITY_CONST(A2_aft[i], S0);
    }
  }
#endif
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, ZeroScaleUpdate, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef typename ScalarTraits<Scalar>::magnitudeType Mag;

  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const Mag M0                = ScalarTraits<Mag>::zero();
  RCP<const Comm<int> > comm  = getDefaultComm();
  // create a Map
  const size_t numLocal                                          = 2;
  const size_t numVectors                                        = 2;
  const size_t LDA                                               = 2;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  Array<Scalar> values(6);
  // values = {1, 1, 2, 2, 4, 4}
  // values(0,4) = {1, 1, 2, 2} = [1 2]
  //                            = [1 2]
  // values(2,6) = {2, 2, 4, 4} = [2 4]
  //                            = [2 4]
  // a multivector A constructed from the first
  // has values .5 of a multivector B constructed from the second
  // then 2*A - B = 0
  // we test both scale(), both update(), and norm()
  values[0] = as<Scalar>(1);
  values[1] = as<Scalar>(1);
  values[2] = as<Scalar>(2);
  values[3] = as<Scalar>(2);
  values[4] = as<Scalar>(4);
  values[5] = as<Scalar>(4);
  MV A(map, values(0, 4), LDA, numVectors),
      B(map, values(2, 4), LDA, numVectors);
  Array<Mag> norms(numVectors), zeros(numVectors);
  std::fill(zeros.begin(), zeros.end(), M0);
  //
  //      [.... ....]
  // A == [ones ones]
  //      [.... ....]
  //
  //      [.... ....]
  // B == [twos twos]
  //      [.... ....]
  //
  //   set A2 = A
  //   scale it by 2 in situ
  //   check that it equals B: subtraction in situ
  {
    MV A2(A);
    A2.scale(as<Scalar>(2));
    A2.update(as<Scalar>(-1), B, as<Scalar>(1));
    A2.norm2(norms);
    TEST_COMPARE_FLOATING_ARRAYS(norms, zeros, M0);
  }
  //   set A2 = A
  //   check that it equals B: scale,subtraction in situ
  {
    MV A2(A);
    A2.update(as<Scalar>(-1), B, as<Scalar>(2));
    A2.norm2(norms);
    // TODO:FAILED TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
  }
  //   set C random
  //   set it to zero by combination with A,B
  {
    MV C(map, numVectors);
    C.randomize();
    C.update(as<Scalar>(-1), B, as<Scalar>(2), A, as<Scalar>(0));
    C.norm2(norms);
    // TODO:FAILED      TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
  }
  //   set C random
  //   scale it ex-situ
  //   check that it equals B: subtraction in situ
  // TODO this is only available with Tpetra???
  /*{
    MV C(map,numVectors);
    C.scale(as<Scalar>(2),A);
    C.update(as<Scalar>(1),B,as<Scalar>(-1));
    C.norm2(norms);
    //TODO:FAILED  TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,M0);
  }*/
#endif  // HAVE_XPETRA_TPETRA
}

////
#if 0  // TAW fix me
  TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL( MultiVector, ScaleAndAssign, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node )
  {
#ifdef HAVE_XPETRA_TPETRA
    using std::endl;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef typename STS::magnitudeType Mag;
    typedef Teuchos::ScalarTraits<Mag> STM;
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
    typedef Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vec_type;

    if (STS::isOrdinal) {
      return;
    }

    STS::seedrandom (0); // consistent seed
    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    const Mag tol = errorTolSlack * STM::eps ();
    const Mag M0 = STM::zero ();

    RCP<const Comm<int> > comm = getDefaultComm ();

    // create a Map
    const size_t numLocal = 23;
    const size_t numVectors = 11;
    RCP<const map_type> map =
      createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node> (INVALID, numLocal, comm);

    // Use random multivector A
    // Set B = A * 2 manually.
    // Therefore, if C = 2*A, then C == B
    // If C = A and C *= 2, then C == B
    // This test operator= and all of our scale ops
    // We'll do Vector and MultiVector variations
    // Also, ensure that other vectors aren't changed
    MV A (map, numVectors, false);
    MV B (map, numVectors, false);
    A.randomize();
    Array<Mag> Anrms(numVectors);
    A.norm2(Anrms());
    // Set B = A * 2, using different techniques, depending on the value of j % 4.
    //
    // 0: getVector, B_j = A_j (Vector::operator=), B_j(i) *= 2
    // 1: getVector, B_j(i) = 2 * A_j(i)
    // 2: getVector, B_j->update (2, *A_j, 0)
    // 3. A_j = A.getData(j), B_j = B.getDataNonConst(j), B_j(i) = A_j(i) * 2
    TEUCHOS_TEST_FOR_EXCEPT(numVectors < 4);

    for (size_t j = 0; j < numVectors; ++j) {
      // assign j-th vector of B to 2 * j-th vector of A
      switch (j % 4) {
      case 0:
        {
          // B(:,j) := A(:,j), and B(i,j) *= 2 for all i.
          RCP<vec_type> bj = B.getVectorNonConst (j);
          RCP<const vec_type> aj = A.getVector (j);
          (*bj) = (*aj);
          ArrayRCP<Scalar> bjview = bj->getDataNonConst (0); // zero-th column of bj
          for (size_t i = 0; i < numLocal; ++i) {
            bjview[i] *= as<Scalar> (2);
          }
        }
        break;
      case 1:
        {
          // B(i,j) := 2 * A(i,j) for all i.
          RCP<vec_type> B_j = B.getVectorNonConst (j);
          RCP<const vec_type> A_j = A.getVector (j);
          // View of the zero-th column of A_j is a view of the j-th column of A.
          ArrayRCP<const Scalar> A_j_view = A_j->getData (0);
          // View of the zero-th column of B_j is a view of the j-th column of B.
          ArrayRCP<Scalar> B_j_view = B_j->getDataNonConst (0);
          for (size_t i = 0; i < numLocal; ++i) {
            B_j_view[i] = as<Scalar> (2) * A_j_view[i];
          }
        }
        break;
      case 2:
        {
          RCP<vec_type> B_j = B.getVectorNonConst (j);
          RCP<const vec_type> A_j = A.getVector (j);
          B_j->update (as<Scalar> (2), *A_j, STS::zero ());
        }
        break;
      case 3:
      default:
        {
          ArrayRCP<Scalar>       bjview = B.getDataNonConst(j);
          ArrayRCP<const Scalar> ajview = A.getData(j);
          for (size_t i=0; i < numLocal; ++i) {
            bjview[i] = as<Scalar>(2) * ajview[i];
          }
        }
        break;
      }
    }

    // Check that A wasn't modified
    out << "Check that A wasn't modified" << endl;
    {
      Teuchos::OSTab tab1 (out);
      Array<Mag> Anrms_aft(numVectors);
      A.norm2(Anrms_aft());
      TEST_COMPARE_FLOATING_ARRAYS(Anrms(),Anrms_aft(),tol);
    }
    // Check that C.Scale(2, A) results in C == B
    out << "Check that C.scale(2, A) results in C == B" << endl;
    {
      Teuchos::OSTab tab1 (out);
      MV C (map, numVectors, false);
      C.scale (as<Scalar> (2), A);
      C.update (-1.0, B, 1.0); // C := C - B
      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm2(Cnorms());
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms(),zeros,tol);
    }
    out << "Check that A wasn't modified" << endl;
    {
      Teuchos::OSTab tab1 (out);
      Array<Mag> Anrms_aft (numVectors);
      A.norm2 (Anrms_aft ());
      TEST_COMPARE_FLOATING_ARRAYS( Anrms (), Anrms_aft (), tol);
    }
    // Check that C = A, C.scale(2) results in C == B
    out << "Check that C = A, C.scale(2) results in C == B" << endl;
    {
      Teuchos::OSTab tab1 (out);
      MV C (map, numVectors, false);
      C = A;
      C.scale (as<Scalar> (2));
      C.update (-1.0, B, 1.0);
      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm2(Cnorms());
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms(),zeros,tol);
    }
    out << "Check that A wasn't modified" << endl;
    {
      Teuchos::OSTab tab1 (out);
      Array<Mag> Anrms_aft (numVectors);
      A.norm2 (Anrms_aft ());
      TEST_COMPARE_FLOATING_ARRAYS( Anrms (), Anrms_aft (), tol);
    }
    // Check that C = A, C_j.scale(2) for all j, results in C == B
    out << "Check that C = A, C_j.scale(2) for all j, results in C == B" << endl;
    {
      Teuchos::OSTab tab1 (out);
      MV C(map,numVectors,false);
      C = A;
      for (size_t j = 0; j < numVectors; ++j) {
        RCP<vec_type> C_j = C.getVectorNonConst (j);
        C_j->scale (as<Scalar> (2));
      }
      C.update (-1.0, B, 1.0);
      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm2(Cnorms());
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms(),zeros,tol);
    }
    out << "Check that A wasn't modified" << endl;
    {
      Teuchos::OSTab tab1 (out);
      Array<Mag> Anrms_aft (numVectors);
      A.norm2 (Anrms_aft ());
      TEST_COMPARE_FLOATING_ARRAYS( Anrms (), Anrms_aft (), tol);
    }
    // Check that C = A, C.scale([2, 2, ..., 2]) results in C == B
    out << "Check that C = A, C.scale([2, 2, ..., 2]) results in C == B" << endl;
    {
      Teuchos::OSTab tab1 (out);
      MV C(map,numVectors,false);
      C = A;
      Array<Scalar> twos(numVectors,as<Scalar>(2));
      C.scale(twos());
      C.update(-1.0,B,1.0);
      Array<Mag> Cnorms(numVectors), zeros(numVectors,M0);
      C.norm2(Cnorms());
      TEST_COMPARE_FLOATING_ARRAYS(Cnorms(),zeros,tol);
    }
    out << "Check that A wasn't modified" << endl;
    {
      Teuchos::OSTab tab1 (out);
      Array<Mag> Anrms_aft (numVectors);
      A.norm2 (Anrms_aft ());
      TEST_COMPARE_FLOATING_ARRAYS( Anrms (), Anrms_aft (), tol);
    }
#endif  // HAVE_XPETRA_TPETRA
  }
#endif  // fix me!

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(Vector, ZeroScaleUpdate, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const Mag M0                = ScalarTraits<Mag>::zero();
  RCP<const Comm<int> > comm  = getDefaultComm();
  // create a Map
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 2, comm);
  Array<Scalar> values(6);
  // values = {1, 1, 2, 2}
  // values(0,2) = {1, 1} = [1]
  //                      = [1]
  // values(2,2) = {2, 2} = [2]
  //                      = [2]
  // a vector A constructed from the first
  // has values .5 of a vector B constructed from the second
  // thus 2*A - B = 0
  // we test both scale(), both update(), and norm()
  values[0] = as<Scalar>(1);
  values[1] = as<Scalar>(1);
  values[2] = as<Scalar>(2);
  values[3] = as<Scalar>(2);
  V A(map, values(0, 2)),
      B(map, values(2, 2));
  Mag norm;
  Array<Mag> norms(1);
  //
  //      [....]
  // A == [ones]
  //      [....]
  //
  //      [....]
  // B == [twos]
  //      [....]
  //
  //   set A2 = A
  //   scale it by 2 in situ
  //   check that it equals B: subtraction in situ
  {
    V A2(A);
    A2.scale(as<Scalar>(2));
    A2.update(as<Scalar>(-1), B, as<Scalar>(1));
    norm = A2.norm2();
    A2.norm2(norms());
    TEST_EQUALITY(norm, M0);
    TEST_EQUALITY(norm, norms[0]);
  }
  //   set A2 = A
  //   check that it equals B: scale,subtraction in situ
  {
    V A2(A);
    A2.update(as<Scalar>(-1), B, as<Scalar>(2));
    norm = A2.norm2();
    A2.norm2(norms());
    // TODO FAILED: TEST_EQUALITY(norm,M0);
    // TODO FAILED: TEST_EQUALITY(norm,norms[0]);
  }
  //   set C random
  //   set it to zero by combination with A,B
  {
    V C(map);
    C.randomize();
    C.update(as<Scalar>(-1), B, as<Scalar>(2), A, as<Scalar>(0));
    norm = C.norm2();
    C.norm2(norms());
    // TODO FAILED: TEST_EQUALITY(norm,M0);
    // TODO FAILED: TEST_EQUALITY(norm,norms[0]);
  }
  //   set C random
  //   scale it ex-situ
  //   check that it equals B: subtraction in situ
  {  // TODO only available with Tpetra??
    V C(map);
    C.randomize();
    C.scale(as<Scalar>(2), A);
    C.update(as<Scalar>(1), B, as<Scalar>(-1));
    norm = C.norm2();
    C.norm2(norms());
    // TODO FAILED: TEST_EQUALITY(norm,M0);
    // TODO FAILED: TEST_EQUALITY(norm,norms[0]);
  }
#endif  // HAVE_XPETRA_TPETRA
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, CopyConst, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA

  typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const Mag M0                = ScalarTraits<Mag>::zero();
  RCP<const Comm<int> > comm  = getDefaultComm();
  // create a Map
  const size_t numLocal                                          = 13;
  const size_t numVectors                                        = 7;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  {
    // create random MV
    MV mvorig(map, numVectors);
    mvorig.randomize();
    // create non-const subview, test copy constructor
    TEUCHOS_TEST_FOR_EXCEPT(numVectors != 7);
    Tuple<size_t, 3> inds = tuple<size_t>(1, 3, 5);
#ifdef XPETRA_NOT_IMPLEMENTED
    RCP<MV> mvview = mvorig.subViewNonConst(inds);
#endif
    Array<Mag> norig(numVectors), nsub(inds.size()), ncopy(inds.size());
    mvorig.normInf(norig());
    for (size_t j = 0; j < as<size_t>(inds.size()); ++j) {
      nsub[j] = norig[inds[j]];
    }
#ifdef XPETRA_NOT_IMPLEMENTED
    MV mvcopy(*mvview);
    mvcopy.normInf(ncopy());
    TEST_COMPARE_FLOATING_ARRAYS(ncopy, nsub, M0);
    // reset both the view and the copy of the view, ensure that they are independent
    Teuchos::Array<Mag> nsub_aft(inds.size()), ones(inds.size(), as<Mag>(1));
    Teuchos::Array<Mag> ncopy_aft(inds.size()), twos(inds.size(), as<Mag>(2));
    mvview->putScalar(as<Scalar>(1));
    mvcopy.putScalar(as<Scalar>(2));
    mvview->normInf(nsub_aft());
    mvcopy.normInf(ncopy_aft());
    TEST_COMPARE_FLOATING_ARRAYS(nsub_aft, ones, M0);
    TEST_COMPARE_FLOATING_ARRAYS(ncopy_aft, twos, M0);
#endif
  }
  {
    // create random MV
    MV morig(map, numVectors);
    morig.randomize();
    // test copy constructor with
    // copy it
    MV mcopy1(morig), mcopy2(morig);
    // verify that all three have identical values
    Array<Mag> norig(numVectors), ncopy1(numVectors), ncopy2(numVectors);
    morig.normInf(norig);
    mcopy1.normInf(ncopy1);
    mcopy2.normInf(ncopy2);
    TEST_COMPARE_FLOATING_ARRAYS(norig, ncopy1, M0);
    TEST_COMPARE_FLOATING_ARRAYS(norig, ncopy2, M0);
    // modify all three
    morig.putScalar(as<Scalar>(0));
    mcopy1.putScalar(as<Scalar>(1));
    mcopy2.putScalar(as<Scalar>(2));
    // compute norms, check
    Array<Mag> zeros(numVectors, as<Mag>(0)), ones(numVectors, as<Mag>(1)), twos(numVectors, as<Mag>(2));
    morig.normInf(norig);
    mcopy1.normInf(ncopy1);
    mcopy2.normInf(ncopy2);
    // TODO:FAILED TEST_COMPARE_FLOATING_ARRAYS(norig,zeros,M0);
    //             TEST_COMPARE_FLOATING_ARRAYS(ncopy1,ones,M0);
    //             TEST_COMPARE_FLOATING_ARRAYS(ncopy2,twos,M0);
  }
#endif  // HAVE_XPETRA_TPETRA
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(Vector, CopyConst, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  // create a Map
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 2, comm);
  // create random MV
  V morig(map);
  morig.randomize();
  // copy it
  V mcopy1(morig), mcopy2(morig);
  // verify that all three have identical values
  Magnitude norig, ncopy1, ncopy2;
  norig  = morig.normInf();
  ncopy1 = mcopy1.normInf();
  ncopy2 = mcopy2.normInf();
  TEST_EQUALITY(norig, ncopy1);
  TEST_EQUALITY(norig, ncopy2);
  TEST_EQUALITY(ncopy1, ncopy2);
  // modify all three
  morig.putScalar(as<Scalar>(0));
  mcopy1.putScalar(as<Scalar>(1));
  mcopy2.putScalar(as<Scalar>(2));
  // compute norms
  norig  = morig.normInf();
  ncopy1 = mcopy1.normInf();
  ncopy2 = mcopy2.normInf();
  // check them
  // TODO FAILED: TEST_EQUALITY(norig, as<Scalar>(0));
  // TODO FAILED: TEST_EQUALITY(ncopy1,as<Scalar>(1));
  // TODO FAILED: TEST_EQUALITY(ncopy2,as<Scalar>(2));
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(Vector, Indexing, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef ScalarTraits<Scalar> SCT;
  typedef typename SCT::magnitudeType Magnitude;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  // create a Map
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 100, comm);
  // create two random Vector objects
  V v1(map), v2(map);
  v1.randomize();
  v2.randomize();
  // set values in both to 1.0
  // for the first, do via putScalar()
  // the the second, do manually, looping over all elements
  // verify that both have identical values
  v1.putScalar(SCT::one());
  {
    ArrayRCP<Scalar> view = v2.get1dViewNonConst();
    for (typename ArrayRCP<Scalar>::iterator v = view.begin(); v != view.end(); ++v) {
      *v = SCT::one();
    }
    view = Teuchos::null;
  }
  Magnitude err;
  // subtract v2 from v1; this should result in v1 = zeros
  v1.update(-1.0, v2, 1.0);
  err = v1.norm2();
  TEST_EQUALITY_CONST(err, SCT::zero());
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, SingleVecNormalize, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  // this documents a usage case in Anasazi::SVQBOrthoManager, which was failing
  // error turned out to be a neglected return in both implementations of update(),
  // after passing the buck to scale() in the case of alpha==0 or beta==0 or gamma=0
  if (ScalarTraits<Scalar>::isOrdinal) return;

  typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const Magnitude M1          = ScalarTraits<Magnitude>::one();
  // TODO unused: const Magnitude M0 = ScalarTraits<Magnitude>::zero();
  RCP<const Comm<int> > comm = getDefaultComm();
  // create a Map
  const size_t numLocal                                          = 10;
  const size_t numVectors                                        = 6;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  // create random MV
  MV mv(map, numVectors);
  mv.randomize();
  // compute the norms
  Array<Magnitude> norms(numVectors);
  mv.norm2(norms());
#ifdef XPETRA_NOT_IMPLEMENTED
  for (size_t j = 0; j < numVectors; ++j) {
    // get a view of column j, normalize it using update()
    RCP<MV> mvj = mv.subViewNonConst(tuple<size_t>(j));
    switch (j) {
      case 0:
        mvj->scale(M1 / norms[j]);
        break;
      case 1:
        mvj->update(M1 / norms[j], *mvj, M0);
        break;
      case 2:
        mvj->update(M0, *mvj, M1 / norms[j]);
        break;
      case 3:
        mvj->update(M0, *mvj, M1 / norms[j], *mvj, M0);
        break;
      case 4:
        mvj->update(M1 / norms[j], *mvj, M0, *mvj, M0);
        break;
      case 5:
        mvj->update(M0, *mvj, M0, *mvj, M1 / norms[j]);
        break;
    }
  }
#endif
  mv.norm2(norms());  // should be all one now
  Array<Magnitude> ones(numVectors, M1);
  // TODO FAILED: TEST_COMPARE_FLOATING_ARRAYS(norms,ones,ScalarTraits<Magnitude>::eps()*as<Magnitude>(10.));
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, CountDot, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA

  typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const Magnitude M0          = ScalarTraits<Magnitude>::zero();
  RCP<const Comm<int> > comm  = getDefaultComm();
  const int numImages         = comm->getSize();
  // create a Map
  const size_t numLocal                                          = 2;
  const size_t numVectors                                        = 3;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  Array<Scalar> values(6);
  // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
  //                               [0 1 2]
  // dot(values,values) = [0*0+0*0 1+1*1 + 2*2+2*2] = [0 2 8]
  // summed over all procs, this is [0 2*nprocs 8*nprocs]
  values[0] = as<Scalar>(0);
  values[1] = as<Scalar>(0);
  values[2] = as<Scalar>(1);
  values[3] = as<Scalar>(1);
  values[4] = as<Scalar>(2);
  values[5] = as<Scalar>(2);
  MV mvec1(map, values(), 2, numVectors),
      mvec2(map, values(), 2, numVectors);
  Array<Scalar> dots1(numVectors), dots2(numVectors), answer(numVectors);
  answer[0] = as<Scalar>(0);
  answer[1] = as<Scalar>(2 * numImages);
  answer[2] = as<Scalar>(8 * numImages);
  // do the dots
  mvec1.dot(mvec2, dots1());
  mvec2.dot(mvec1, dots2());
  // check the answers
  TEST_COMPARE_FLOATING_ARRAYS(dots1, dots2, M0);
  TEST_COMPARE_FLOATING_ARRAYS(dots1, answer, M0);
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, CountDotNonTrivLDA, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  // same as CountDot, but the A,LDA has a non-trivial LDA (i.e., LDA != myLen)

  typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const Magnitude M0          = ScalarTraits<Magnitude>::zero();
  RCP<const Comm<int> > comm  = getDefaultComm();
  const int numImages         = comm->getSize();
  // create a Map
  const size_t numLocal                                          = 2;
  const size_t numVectors                                        = 3;
  const size_t LDA                                               = 3;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  Array<Scalar> values(9);
  // A = {0, 0, -1, 1, 1, -1, 2, 2, -1} = [0   1  2]
  //                                      [0   1  2]
  //                                      [-1 -1 -1]
  // processed as a 2 x 3 with LDA==3, the result it
  //            values =       [0 1 2]
  //                           [0 1 2]
  // dot(values,values) = [0*0+0*0 1+1*1 + 2*2+2*2] = [0 2 8]
  // summed over all procs, this is [0 2*nprocs 8*nprocs]
  values[0] = as<Scalar>(0);
  values[1] = as<Scalar>(0);
  values[2] = as<Scalar>(-1);
  values[3] = as<Scalar>(1);
  values[4] = as<Scalar>(1);
  values[5] = as<Scalar>(-1);
  values[6] = as<Scalar>(2);
  values[7] = as<Scalar>(2);
  values[8] = as<Scalar>(-1);
  MV mvec1(map, values(), LDA, numVectors),
      mvec2(map, values(), LDA, numVectors);
  Array<Scalar> dots1(numVectors), dots2(numVectors), answer(numVectors);
  answer[0] = as<Scalar>(0);
  answer[1] = as<Scalar>(2 * numImages);
  answer[2] = as<Scalar>(8 * numImages);
  // do the dots
  mvec1.dot(mvec2, dots1());
  mvec2.dot(mvec1, dots2());
  // check the answers
  TEST_COMPARE_FLOATING_ARRAYS(dots1, dots2, M0);
  TEST_COMPARE_FLOATING_ARRAYS(dots1, answer, M0);
#endif  // HAVE_XPETRA_TPETRA
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, CountNorm1, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA

  typedef typename ScalarTraits<Scalar>::magnitudeType MT;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const MT M0                 = ScalarTraits<MT>::zero();
  RCP<const Comm<int> > comm  = getDefaultComm();
  const int numImages         = comm->getSize();
  // create a Map
  const size_t numLocal                                          = 2;
  const size_t numVectors                                        = 3;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  Array<Scalar> values(6);
  // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
  //                               [0 1 2]
  // norm1(values) = [0 2 4]
  // over all procs, this is [0 2*nprocs 4*nprocs]
  // mean is [0 1 2]
  values[0] = as<Scalar>(0);
  values[1] = as<Scalar>(0);
  values[2] = as<Scalar>(1);
  values[3] = as<Scalar>(1);
  values[4] = as<Scalar>(2);
  values[5] = as<Scalar>(2);
  MV mvec(map, values(), 2, numVectors);
  // compute, check norms
  {
    Array<MT> norms(numVectors), answer(numVectors);
    answer[0] = as<MT>(0);
    answer[1] = as<MT>(2 * numImages);
    answer[2] = as<MT>(4 * numImages);
    mvec.norm1(norms());
    TEST_COMPARE_FLOATING_ARRAYS(norms, answer, M0);
  }
  {
    // compute, check means
    Array<Scalar> means(numVectors), answer(numVectors);
    mvec.meanValue(means());
    answer[0] = as<Scalar>(0);
    answer[1] = as<Scalar>(1);
    answer[2] = as<Scalar>(2);
    TEST_COMPARE_FLOATING_ARRAYS(means, answer, M0);
    for (size_t j = 0; j < numVectors; ++j) {
      // TODO FAILED: TEST_EQUALITY( mvec.getVector(j)->meanValue(), answer[j] );
    }
  }
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, CountNormInf, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA

  typedef typename ScalarTraits<Scalar>::magnitudeType MT;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const MT M0                 = ScalarTraits<MT>::zero();
  RCP<const Comm<int> > comm  = getDefaultComm();
  // create a Map
  const size_t numLocal                                          = 2;
  const size_t numVectors                                        = 3;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  Array<Scalar> values(6);
  // values = {0, 0, 1, 1, 2, 2} = [0 1 2]
  //                               [0 1 2]
  // normInf(values) = [0 1 2]
  // over all procs, this is [0 1 2]
  values[0] = as<Scalar>(0);
  values[1] = as<Scalar>(0);
  values[2] = as<Scalar>(1);
  values[3] = as<Scalar>(1);
  values[4] = as<Scalar>(2);
  values[5] = as<Scalar>(2);
  MV mvec(map, values(), 2, numVectors);
  Array<MT> norms(numVectors), answer(numVectors);
  answer[0] = as<MT>(0);
  answer[1] = as<MT>(1);
  answer[2] = as<MT>(2);
  // do the dots
  mvec.normInf(norms());
  // check the answers
  TEST_COMPARE_FLOATING_ARRAYS(norms, answer, M0);
#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, Norm2, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef typename ScalarTraits<Scalar>::magnitudeType MT;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  const MT M0                 = ScalarTraits<MT>::zero();
  RCP<const Comm<int> > comm  = getDefaultComm();
  // create a Map
  const size_t numLocal                                          = 13;
  const size_t numVectors                                        = 7;
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, numLocal, comm);
  MV mvec(map, numVectors);
  // randomize the multivector
  mvec.randomize();
  // take norms; they should not be zero
  Array<MT> normsRand(numVectors), normsZero(numVectors);
  mvec.norm2(normsRand());
  // zero the vector
  mvec.putScalar(ScalarTraits<Scalar>::zero());
  // take norms; they should be zero
  mvec.norm2(normsZero());
  // check the answers
  bool local_success = true;
  for (size_t i = 0; i < numVectors; ++i) {
    TEST_ARRAY_ELE_INEQUALITY(normsRand, i, M0);
    TEST_ARRAY_ELE_EQUALITY(normsZero, i, M0);
  }
  success &= local_success;
#endif
}

//// TODO this code should be generalized
TEUCHOS_UNIT_TEST_TEMPLATE_6_DECL(MultiVector, BadCombinations, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
  const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
  RCP<const Comm<int> > comm  = getDefaultComm();
  const int myImageID         = comm->getRank();
  // create a Map
  const Scalar rnd = ScalarTraits<Scalar>::random();
  // two maps: one has two entries per process, the other disagrees on process 0
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map1 = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, 2, comm),
                                                             map2 = createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(INVALID, myImageID == 0 ? 1 : 2, comm);
  // multivectors from different maps are incompatible for all ops
  // multivectors from the same map are compatible only if they have the same number of
  //    columns
  MV m1n1(map1, 1), m1n2(map1, 2), m2n2(map2, 2), m1n2_2(map1, 2);
  Array<Scalar> dots(1);
  Array<Mag> norms(1);
  // TODO: what exceptions are on the other processors thrown?
  if (myImageID == 0) {
    // FINISH: test multiply (both), reciprocalMultiply
    TEST_THROW(m1n2.dot(m1n1, dots()), std::runtime_error);  // dot
    TEST_THROW(m1n2.dot(m2n2, dots()), std::runtime_error);
    TEST_THROW(m1n2.abs(m1n1), std::runtime_error);  // abs
    TEST_THROW(m1n2.abs(m2n2), std::runtime_error);
    TEST_THROW(m1n2.abs(m1n1), std::runtime_error);  // abs
    TEST_THROW(m1n2.abs(m2n2), std::runtime_error);
    // TEST_THROW(m1n2.scale(rnd,m1n1), std::runtime_error); // abs  // TODO only available with Tpetra??
    // TEST_THROW(m1n2.scale(rnd,m2n2), std::runtime_error);
    TEST_THROW(m1n2.update(rnd, m1n1, rnd), std::runtime_error);  // update(alpha,A,beta)
    if (::Tpetra::Details::Behavior::debug()) {
      TEST_THROW(m1n2.update(rnd, m2n2, rnd), std::runtime_error);
      TEST_THROW(m1n2.update(rnd, m2n2, rnd, m1n2_2, rnd), std::runtime_error);  // update(alpha,A,beta,B,gamma) // A incompat
      TEST_THROW(m1n2.update(rnd, m2n2, rnd, m1n2_2, rnd), std::runtime_error);  // incompt is length            // A incompat
      TEST_THROW(m1n2.update(rnd, m1n2_2, rnd, m2n2, rnd), std::runtime_error);  // B incompat
      TEST_THROW(m1n2.update(rnd, m1n2_2, rnd, m2n2, rnd), std::runtime_error);  // B incompat
      TEST_THROW(m1n2.update(rnd, m2n2, rnd, m2n2, rnd), std::runtime_error);    // A,B incompat
      TEST_THROW(m1n2.update(rnd, m2n2, rnd, m2n2, rnd), std::runtime_error);    // A,B incompat
      TEST_THROW(m1n2.update(rnd, m1n1, rnd, m1n2_2, rnd), std::runtime_error);  // incompt is numVecs           // A incompat
      TEST_THROW(m1n2.update(rnd, m1n1, rnd, m1n2_2, rnd), std::runtime_error);  // A incompat
    }
    TEST_THROW(m1n2.update(rnd, m1n2_2, rnd, m1n1, rnd), std::runtime_error);  // B incompat
    TEST_THROW(m1n2.update(rnd, m1n2_2, rnd, m1n1, rnd), std::runtime_error);  // B incompat
    TEST_THROW(m1n2.update(rnd, m1n1, rnd, m1n1, rnd), std::runtime_error);    // A,B incompat
    TEST_THROW(m1n2.update(rnd, m1n1, rnd, m1n1, rnd), std::runtime_error);    // A,B incompat
    TEST_THROW(m1n2.reciprocal(m1n1), std::runtime_error);                     // reciprocal
    TEST_THROW(m1n2.reciprocal(m2n2), std::runtime_error);
  }
#endif
}

TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, Constructor_Epetra, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_EPETRA
  RCP<const Teuchos::Comm<int> > comm = getDefaultComm();

  {
    TEST_NOTHROW(M(10, 0, comm));
    TEST_NOTHROW(MV(Teuchos::rcp(new M(10, 0, comm)), 3));
  }

#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_PTHREAD)
  {
    typedef Xpetra::EpetraMapT<GlobalOrdinal, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> mm;
    TEST_THROW(mm(10, 0, comm), Xpetra::Exceptions::RuntimeError);
    typedef Xpetra::EpetraMultiVectorT<GlobalOrdinal, Tpetra::KokkosCompat::KokkosThreadsWrapperNode> mx;
    TEST_THROW(mx(Teuchos::null, 3), Xpetra::Exceptions::RuntimeError);
  }
#endif
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_CUDA)
  {
    typedef Xpetra::EpetraMapT<GlobalOrdinal, Tpetra::KokkosCompat::KokkosCudaWrapperNode> mm;
    TEST_THROW(mm(10, 0, comm), Xpetra::Exceptions::RuntimeError);
    typedef Xpetra::EpetraMultiVectorT<GlobalOrdinal, Tpetra::KokkosCompat::KokkosCudaWrapperNode> mx;
    TEST_THROW(mx(Teuchos::null, 3), Xpetra::Exceptions::RuntimeError);
  }
#endif
#if defined(HAVE_XPETRA_TPETRA) && defined(HAVE_TPETRA_INST_HIP)
  {
    typedef Xpetra::EpetraMapT<GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> mm;
    TEST_THROW(mm(10, 0, comm), Xpetra::Exceptions::RuntimeError);
    typedef Xpetra::EpetraMultiVectorT<GlobalOrdinal, Tpetra::KokkosCompat::KokkosHIPWrapperNode> mx;
    TEST_THROW(mx(Teuchos::null, 3), Xpetra::Exceptions::RuntimeError);
  }
#endif

#endif
}

////
TEUCHOS_UNIT_TEST_TEMPLATE_7_DECL(MultiVector, Typedefs, M, MV, V, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#ifdef HAVE_XPETRA_TPETRA
  typedef typename MV::scalar_type scalar_type;
  typedef typename MV::local_ordinal_type local_ordinal_type;
  typedef typename MV::global_ordinal_type global_ordinal_type;
  typedef typename MV::node_type node_type;
  TEST_EQUALITY_CONST((std::is_same_v<scalar_type, Scalar>) == true, true);
  TEST_EQUALITY_CONST((std::is_same_v<local_ordinal_type, LocalOrdinal>) == true, true);
  TEST_EQUALITY_CONST((std::is_same_v<global_ordinal_type, GlobalOrdinal>) == true, true);
  TEST_EQUALITY_CONST((std::is_same_v<node_type, Node>) == true, true);
#endif
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

#define XPETRA_EPETRA_NO_ORDINAL_SCALAR_TYPES(S, LO, GO, N)            \
  typedef typename Xpetra::EpetraMapT<GO, N> M##LO##GO##N;             \
  typedef typename Xpetra::EpetraMultiVectorT<GO, N> MV##S##LO##GO##N; \
  typedef typename Xpetra::EpetraVectorT<GO, N> V##S##LO##GO##N;

#define XPETRA_EPETRA_ORDINAL_SCALAR_TYPES(S, LO, GO, N)                  \
  typedef typename Xpetra::EpetraIntMultiVectorT<GO, N> MV##S##LO##GO##N; \
  typedef typename Xpetra::EpetraIntVectorT<GO, N> V##S##LO##GO##N;

#endif

// List of tests which run only with Tpetra
#define XP_TPETRA_MULTIVECTOR_INSTANT(S, LO, GO, N)                                                                                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, BadConstLDA, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, Describable, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, ElementWiseMultiply, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, CopyConst, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                 \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(Vector, CopyConst, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(Vector, Indexing, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                       \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, OrthoDot, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, CountDot, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, CountDotNonTrivLDA, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, CountNorm1, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, CountNormInf, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, Norm2, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, CopyView, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, OffsetView, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, ZeroScaleUpdate, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(Vector, ZeroScaleUpdate, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, BadMultiply, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                               \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, SingleVecNormalize, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, Multiply, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_6_INSTANT(MultiVector, BadCombinations, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)                           \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, NonMemberConstructorsTpetra, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(Vector, AssignmentDeepCopies, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)

// List of tests which run only with Epetra
#define XP_EPETRA_MULTIVECTOR_INSTANT(S, LO, GO, N)                                                                                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, NonMemberConstructorsEpetra, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, Constructor_Epetra, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)

// list of all tests which run both with Epetra and Tpetra
// TODO: move more lists from the upper list to this list
#define XP_MULTIVECTOR_NO_ORDINAL_INSTANT(S, LO, GO, N)                                                                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, AssignmentDeepCopies, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, GetVector, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)            \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, BadConstNumVecs, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, BadConstAA, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)

#define XP_MULTIVECTOR_INSTANT(S, LO, GO, N)                                                                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, basic, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N) \
  TEUCHOS_UNIT_TEST_TEMPLATE_7_INSTANT(MultiVector, Typedefs, M##LO##GO##N, MV##S##LO##GO##N, V##S##LO##GO##N, S, LO, GO, N)

// can we relax the INT INT?
#if defined(HAVE_XPETRA_TPETRA)

#include <TpetraCore_config.h>
#include <TpetraCore_ETIHelperMacros.h>

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLGN(XPETRA_TPETRA_TYPES)
TPETRA_INSTANTIATE_SLGN(XP_MULTIVECTOR_INSTANT)
// no ordinal types as scalar for testing as some tests use ScalarTraits::eps...
// TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR ( XPETRA_TPETRA_NO_ORDINAL_TYPES )
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_MULTIVECTOR_NO_ORDINAL_INSTANT)
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(XP_TPETRA_MULTIVECTOR_INSTANT)

#endif

#if defined(HAVE_XPETRA_EPETRA)

#include "Xpetra_Map.hpp"  // defines EpetraNode
typedef Xpetra::EpetraNode EpetraNode;
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
XPETRA_EPETRA_NO_ORDINAL_SCALAR_TYPES(double, int, int, EpetraNode)
XPETRA_EPETRA_ORDINAL_SCALAR_TYPES(int, int, int, EpetraNode)
XP_MULTIVECTOR_NO_ORDINAL_INSTANT(double, int, int, EpetraNode)
XP_MULTIVECTOR_INSTANT(int, int, int, EpetraNode)
XP_EPETRA_MULTIVECTOR_INSTANT(double, int, int, EpetraNode)
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
XPETRA_EPETRA_NO_ORDINAL_SCALAR_TYPES(double, int, LongLong, EpetraNode)
XPETRA_EPETRA_ORDINAL_SCALAR_TYPES(int, int, LongLong, EpetraNode)
XP_MULTIVECTOR_INSTANT(int, int, LongLong, EpetraNode)
XP_MULTIVECTOR_NO_ORDINAL_INSTANT(double, int, LongLong, EpetraNode)
XP_EPETRA_MULTIVECTOR_INSTANT(double, int, LongLong, EpetraNode)
#endif
#endif

}  // namespace
