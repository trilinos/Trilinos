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
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "Xpetra_ConfigDefs.hpp" //TODO
#ifdef HAVE_XPETRA_TPETRA
#  include <Tpetra_ConfigDefs.hpp>
#endif // HAVE_XPETRA_TPETRA
#include "Xpetra_DefaultPlatform.hpp" //TODO
#include "Teuchos_as.hpp"

//#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_Cloner.hpp"

#include <type_traits> // NOTE (mfh 20 Jul 2016) part of fix for #508

namespace {
  using std::sort;
  using std::find;

  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  using Xpetra::DefaultPlatform;

  using Xpetra::Map;
  using Xpetra::MapFactory;
  using Xpetra::Matrix;
  using Xpetra::CrsMatrixWrap;

  using Xpetra::viewLabel_t;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  RCP<const Comm<int> > getDefaultComm() {
    if (testMpi)
      return DefaultPlatform::getDefaultPlatform().getComm();
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  //
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Cloner, MapCloneTpetra, LO, GO, N2 )
  {
#ifdef HAVE_XPETRA_TPETRA
    typedef typename KokkosClassic::DefaultNode::DefaultNodeType N1;

    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();

    const size_t numLocal = 10;
    const size_t INVALID = OrdinalTraits<size_t>::invalid();

    RCP<N1> n1(new N1());
    RCP<N2> n2(new N2());

    typedef Map<LO,GO,N1> Map1;
    typedef Map<LO,GO,N2> Map2;

    // create a contiguous uniform distributed map with numLocal entries per node
    RCP<const Map1> map1 = MapFactory<LO,GO,N1>::createContigMap(Xpetra::UseTpetra, INVALID, numLocal, comm);
    RCP<const Map2> map2  = Xpetra::clone(*map1, n2);
    RCP<const Map1> map1b = Xpetra::clone(*map2, n1);
    TEST_EQUALITY_CONST(map1->isCompatible(*map1b), true);
    TEST_EQUALITY_CONST(map1->isSameAs(*map1b),     true);
#endif
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Cloner, MapCloneEpetra, LO, GO, N2 )
  {
#ifdef HAVE_XPETRA_EPETRA
    typedef typename KokkosClassic::DefaultNode::DefaultNodeType N1;

    // NOTE (mfh 20 Jul 2016) This fixes #508.  Epetra only works for
    // Node = Xpetra::EpetraNode.  This typedef is currently defined
    // in Xpetra_Map.hpp.
    if (! std::is_same<N1, Xpetra::EpetraNode>::value ||
        ! std::is_same<N2, Xpetra::EpetraNode>::value) {
      out << "MapCloneEpetra test only makes sense if N1 and N2 are both "
        "Xpetra::EpetraNode.  This is not the case, so I'm leaving the test "
        "early." << std::endl;
      return;
    }
    else {
      // create a comm
      RCP<const Comm<int> > comm = getDefaultComm();

      const size_t numLocal = 10;
      const size_t INVALID = OrdinalTraits<size_t>::invalid();

      RCP<N1> n1(new N1());
      RCP<N2> n2(new N2());

      typedef Map<LO,GO,N1> Map1;
      typedef Map<LO,GO,N2> Map2;
      RCP<const Map1> map1 = MapFactory<LO,GO,N1>::createContigMap(Xpetra::UseEpetra, INVALID, numLocal, comm);

      bool thrown = false;
      try {
        RCP<const Map2> map2 = clone(*map1, n2);
      } catch (...) {
        thrown = true;
      }
      // Check that Epetra throws
      TEST_EQUALITY(thrown, true);
    }
#endif
  }

  // FIXME (mfh 28 Sep 2013) This test does not compile when I run the
  // check-in test script, disable forward packages, and explicitly
  // enable Ifpack2, Xpetra, and Zoltan2.  See the FIXME in the body
  // of the test.  I didn't change anything in Epetra, Tpetra, Xpetra,
  // or Zoltan2, and I only changed Ifpack2::AdditiveSchwarz, so I
  // have no idea why the test is failing.  I'm disabling the test for
  // now so that I can get my work done.
#if 0
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Cloner, MatrixCloneTpetra, Scalar, LO, GO, N2 )
  {
#ifdef HAVE_XPETRA_TPETRA
    typedef typename KokkosClassic::DefaultNode::DefaultNodeType N1;

    RCP<N1> n1(new N1());
    RCP<N2> n2(new N2());

    typedef Matrix<Scalar, LO, GO, N1> Matrix1;
    typedef Matrix<Scalar, LO, GO, N1> Matrix2;

    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 10;
    const size_t INVALID = OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t

    RCP<Matrix1> matrix1(new CrsMatrixWrap<Scalar,LO,GO,N1>(MapFactory<LO,GO,N1>::createContigMap(Xpetra::UseTpetra, INVALID, numLocal, comm), 1));
    matrix1->fillComplete();

    // FIXME (mfh 28 Sep 2013) This line does not compile when I run
    // the check-in test script, disable forward packages, and
    // explicitly enable Ifpack2, Xpetra, and Zoltan2.  The compiler
    // says that the result of the right-hand side has a different
    // Node type than the left-hand side.
    RCP<Matrix2> matrix2  = clone(*matrix1, n2);
    RCP<Matrix1> matrix1b = clone(*matrix2, n1);
    TEST_EQUALITY_CONST(matrix1b->getRowMap()->isCompatible(*matrix1->getRowMap()), true);
    TEST_EQUALITY_CONST(matrix1b->getRowMap()->isSameAs(*matrix1->getRowMap()),     true);
#endif
  }
#endif // 0

  // FIXME (mfh 28 Sep 2013) This test does not compile when I run the
  // check-in test script, disable forward packages, and explicitly
  // enable Ifpack2, Xpetra, and Zoltan2.  See the FIXME in the body
  // of the test.  I didn't change anything in Epetra, Tpetra, Xpetra,
  // or Zoltan2, and I only changed Ifpack2::AdditiveSchwarz, so I
  // have no idea why the test is failing.  I'm disabling the test for
  // now so that I can get my work done.
#if 0
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( Cloner, MatrixCloneEpetra, Scalar, LO, GO, N2 )
  {
#ifdef HAVE_XPETRA_EPETRA
    typedef typename KokkosClassic::DefaultNode::DefaultNodeType N1;

    RCP<N1> n1(new N1());
    RCP<N2> n2(new N2());

    typedef Matrix<Scalar, LO, GO, N1> Matrix1;
    typedef Matrix<Scalar, LO, GO, N1> Matrix2;

    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numLocal = 10;
    const size_t INVALID = OrdinalTraits<size_t>::invalid(); // TODO: global_size_t instead of size_t

    RCP<Matrix1> matrix1(new CrsMatrixWrap<Scalar,LO,GO,N1>(MapFactory<LO,GO,N1>::createContigMap(Xpetra::UseEpetra, INVALID, numLocal, comm), 1));
    matrix1->fillComplete();
    bool thrown = false;
    try {
      // FIXME (mfh 28 Sep 2013) This line does not compile when I run
      // the check-in test script, disable forward packages, and
      // explicitly enable Ifpack2, Xpetra, and Zoltan2.  The compiler
      // says that the result of the right-hand side has a different
      // Node type than the left-hand side.
      RCP<Matrix2> matrix2 = clone(*matrix1, n2);
    } catch (...) {
      thrown = true;
    }
    TEST_EQUALITY(thrown, true);
#endif
  }
#endif // 0

  //
  // INSTANTIATIONS
  //

  // mfh 15 Oct 2014: In order to exercise clone(), we have to pick
  // NodeType different than the default Node type.  We can use the
  // macros defined in KokkosClassic to figure out some Node type that
  // differs from the default Node type.  Of course, if only one Node
  // type is defined, we have no choice but to use it.

#if defined(HAVE_TPETRA_DEFAULTNODE_OPENMPNODE)
#if defined(HAVE_TPETRA_INST_PTHREAD)
  typedef Kokkos::Compat::KokkosThreadsWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_OPENMP)
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_CUDA)
  typedef Kokkos::Compat::KokkosCudaWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_SERIAL)
  typedef Kokkos::Compat::KokkosSerialWrapperNode NodeType;
#  else
  // There's only one Node type defined, so we have no choice but to use it.
  typedef KokkosClassic::DefaultNode::DefaultNodeType NodeType;
#  endif
#elif defined(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE)
#  if defined(HAVE_TPETRA_INST_PTHREAD)
  typedef Kokkos::Compat::KokkosThreadsWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_OPENMP)
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode NodeType;
// #  elif defined(HAVE_TPETRA_INST_CUDA)
//   typedef Kokkos::Compat::KokkosCudaWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_SERIAL)
  typedef Kokkos::Compat::KokkosSerialWrapperNode NodeType;
#  else
  // There's only one Node type defined, so we have no choice but to use it.
  typedef KokkosClassic::DefaultNode::DefaultNodeType NodeType;
#  endif
#elif defined(HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE)
#  if defined(HAVE_TPETRA_INST_PTHREAD)
  typedef Kokkos::Compat::KokkosThreadsWrapperNode NodeType;
// #  elif defined(HAVE_TPETRA_INST_OPENMP)
//   typedef Kokkos::Compat::KokkosOpenMPWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_CUDA)
  typedef Kokkos::Compat::KokkosCudaWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_SERIAL)
  typedef Kokkos::Compat::KokkosSerialWrapperNode NodeType;
#  else
  // There's only one Node type defined, so we have no choice but to use it.
  typedef KokkosClassic::DefaultNode::DefaultNodeType NodeType;
#  endif
#elif defined(HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE)
#  if defined(HAVE_TPETRA_INST_OPENMP)
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_CUDA)
  typedef Kokkos::Compat::KokkosCudaWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_SERIAL)
  typedef Kokkos::Compat::KokkosSerialWrapperNode NodeType;
// #  elif defined(HAVE_TPETRA_INST_PTHREAD)
//   typedef Kokkos::Compat::KokkosThreadsWrapperNode NodeType;
#  else
  // There's only one Node type defined, so we have no choice but to use it.
  typedef KokkosClassic::DefaultNode::DefaultNodeType NodeType;
#  endif
#elif defined(HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE)
#  if defined(HAVE_TPETRA_INST_PTHREAD)
  typedef Kokkos::Compat::KokkosThreadsWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_OPENMP)
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_CUDA)
  typedef Kokkos::Compat::KokkosCudaWrapperNode NodeType;
// #  elif defined(HAVE_TPETRA_INST_SERIAL)
//   typedef Kokkos::Compat::KokkosSerialWrapperNode NodeType;
#  else
  // There's only one Node type defined, so we have no choice but to use it.
  typedef KokkosClassic::DefaultNode::DefaultNodeType NodeType;
#  endif
#elif defined(HAVE_TPETRA_DEFAULTNODE_SERIALNODE)
#  if defined(HAVE_TPETRA_INST_PTHREAD)
  typedef Kokkos::Compat::KokkosThreadsWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_OPENMP)
  typedef Kokkos::Compat::KokkosOpenMPWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_CUDA)
  typedef Kokkos::Compat::KokkosCudaWrapperNode NodeType;
#  elif defined(HAVE_TPETRA_INST_SERIAL)
  typedef Kokkos::Compat::KokkosSerialWrapperNode NodeType;
#  else
  // There's only one Node type defined, so we have no choice but to use it.
  typedef KokkosClassic::DefaultNode::DefaultNodeType NodeType;
#  endif
#else
  // Pick some reasonable default.
  typedef KokkosClassic::DefaultNode::DefaultNodeType NodeType;
#endif


#if defined(HAVE_XPETRA_TPETRA)
#ifdef HAVE_TPETRA_INST_INT_INT
        TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Cloner, MapCloneTpetra, int, int, NodeType )
#endif
#ifdef HAVE_TPETRA_INST_INT_LONG
        TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Cloner, MapCloneTpetra, int, long, NodeType )
#endif
#ifdef HAVE_TPETRA_INST_INT_LONG_LONG
        typedef long long LongLong;
        TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Cloner, MapCloneTpetra, int, LongLong, NodeType )
#endif
#endif

#if defined(HAVE_XPETRA_EPETRA)
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Cloner, MapCloneEpetra, int, int, NodeType )
#endif
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
typedef long long LongLong;
TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Cloner, MapCloneEpetra, int, LongLong, NodeType )
#endif
#endif

  // FIXME (mfh 28 Sep 2013) I disabled this test.  Please uncomment
  // the line below if you want to reenable the test.
  //TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Cloner, MatrixCloneTpetra, double, int, int, NodeType )

  // FIXME (mfh 28 Sep 2013) I disabled this test.  Please uncomment
  // the line below if you want to reenable the test.
  //TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Cloner, MatrixCloneEpetra, double, int, int, NodeType )
// #else
// #warning Skipping Cloner tests as KokkosClassic is not enabled
// #endif
}
