// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Details_KokkosCounter.hpp"

// TODO: add test where some nodes have zero rows
// TODO: add test where non-"zero" graph is used to build matrix; if no values are added to matrix, the operator effect should be zero. This tests that matrix values are initialized properly.
// TODO: add test where dynamic profile initially has no allocation, then entries are added. this will test new view functionality.

namespace Teuchos {
  template <>
    ScalarTraits<int>::magnitudeType
    relErr( const int &s1, const int &s2 )
    {
      typedef ScalarTraits<int> ST;
      return ST::magnitude(s1-s2);
    }

  template <>
    ScalarTraits<char>::magnitudeType
    relErr( const char &s1, const char &s2 )
    {
      typedef ScalarTraits<char> ST;
      return ST::magnitude(s1-s2);
    }
}

namespace {

  // no ScalarTraits<>::eps() for integer types

  template <class Scalar, bool hasMachineParameters> struct TestingTolGuts {};

  template <class Scalar>
  struct TestingTolGuts<Scalar, true> {
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType testingTol()
      { return Teuchos::ScalarTraits<Scalar>::eps(); }
  };

  template <class Scalar>
  struct TestingTolGuts<Scalar, false> {
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType testingTol()
      { return 0; }
  };

  template <class Scalar>
  static typename Teuchos::ScalarTraits<Scalar>::magnitudeType testingTol()
  {
    return TestingTolGuts<Scalar, Teuchos::ScalarTraits<Scalar>::hasMachineParameters>::
      testingTol();
  }

  using Tpetra::TestingUtilities::getDefaultComm;

  using std::endl;
  using std::swap;

  using std::string;

  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::arcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Teuchos::ETransp;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::EDiag;
  using Teuchos::UNIT_DIAG;
  using Teuchos::NON_UNIT_DIAG;
  using Teuchos::EUplo;
  using Teuchos::UPPER_TRI;
  using Teuchos::LOWER_TRI;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  using Tpetra::Map;
  using Tpetra::MultiVector;
  using Tpetra::Vector;
  using Tpetra::Operator;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::RowMatrix;
  using Tpetra::Import;
  using Tpetra::global_size_t;
  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;
  using Tpetra::createVector;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, NonSquare, LO, GO, Scalar, Node )
  {

    // skip test if Scalar is not (float or double)
    if constexpr (!(std::is_same_v<Scalar, float> || std::is_same_v<Scalar, double>)) {
      out << "SKIP: unsupported scalar type" << std::endl;
      TEST_EQUALITY_CONST(1,1); // SKIP
      return;
    }
    // skip test if LO != int
    if constexpr (!std::is_same_v<LO, int>) {
      out << "SKIP: unsupported local ordinal type" << std::endl;
      TEST_EQUALITY_CONST(1,1); // SKIP
      return;
    }
    // skip test if CUDA enables and not CUDA space
  #if defined(HAVE_TPETRA_CUDA) || defined(HAVE_TPETRACORE_CUDA)
    if constexpr (!std::is_same_v<typename Node::execution_space, Kokkos::Cuda>) {
      out << "SKIP: non-CUDA exec space" << std::endl;
      TEST_EQUALITY_CONST(1,1); // SKIP
      return;
    }
  #endif
  // skip test if HIP enabled and not HIP space
  #if defined(HAVE_TPETRA_HIP) || defined(HAVE_TPETRACORE_HIP)
    if constexpr (!std::is_same_v<typename Node::execution_space, Kokkos::HIP>) {
      out << "SKIP: non-HIP exec space" << std::endl;
      TEST_EQUALITY_CONST(1,1); // SKIP
      return;
    }
  #endif

    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef Map<LO,GO,Node> map_type;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int M = 3;
    const int P = 5;
    const int N = comm->getSize();
    const int myImageID = comm->getRank();
    // create Maps
    // matrix is M*N-by-P
    //                  col
    //            0        1                  P-1
    //    0  [0        MN              ... (P-1)MN     ]
    //    .  [...      ...                 ...         ]
    //    0  [M-1      MN+M-1              (P-1)MN+M-1 ]
    //p   1  [M        MN+M                            ]
    //r   .  [...      ...                             ] = [A_ij], where A_ij = i+jMN
    //o   1  [2M-1     MN+2M-1                         ]
    //c   .  [...                                      ]
    //   N-1 [(N-1)M   MN+(N-1)(M-1)                   ]
    //    .  [...      ...                             ]
    //   N-1 [MN-1     MN+MN-1                         ]
    //
    // row map, range map is [0,M-1] [M,2M-1] [2M,3M-1] ... [MN-M,MN-1]
    // domain map will be map for X (lclmap)
    //
    // input multivector X is not distributed:
    //
    //   X = [  0    P    ...  (numVecs-1)P ]
    //       [ ...  ....  ...       ...     ] = [X_ji], where X_ij = i+jP
    //       [ P-1  2P-1  ...   numVecs*P-1 ]
    //
    // the result of the non-transpose multiplication should be
    //                              P-1
    // (A*X)_ij = sum_k A_ik X_kj = sum (i+kMN)(k+jP) = jiP^2 + (i+jMNP)(P^2-P)/2 + MNP(P-1)(2P-1)/6
    //                              k=0
    //
    //
    //
    RCP<const map_type> rowmap (new map_type (INVALID, M, 0, comm));
    RCP<const map_type> lclmap = createLocalMapWithNode<LO,GO,Node> (P, comm);

    // create the matrix
    MAT A(rowmap,P);
    for (GO i=0; i<static_cast<GO>(M); ++i) {
      for (GO j=0; j<static_cast<GO>(P); ++j) {
        A.insertGlobalValues( M*myImageID+i, tuple<GO>(j), tuple<Scalar>(M*myImageID+i + j*M*N) );
      }
    }
    // call fillComplete()
    A.fillComplete(lclmap, rowmap);

    {
      // build the input multivector X with 1 vector
      const int numVecs  = 1;
      MV X(lclmap,numVecs);
      for (GO i=0; i<static_cast<GO>(P); ++i) {
        for (GO j=0; j<static_cast<GO>(numVecs); ++j) {
          X.replaceGlobalValue(i,j,static_cast<Scalar>(i+j*P));
        }
      }
      // allocate output multivec
      MV Bout(rowmap,numVecs);
      // test the action
      Bout.randomize();
      Tpetra::Details::KokkosRegionCounter::reset();
      Tpetra::Details::KokkosRegionCounter::start();
      A.apply(X,Bout);
      Tpetra::Details::KokkosRegionCounter::stop();

      TEST_COMPARE(Tpetra::Details::KokkosRegionCounter::get_count_region_contains("spmv[TPL_"), ==, 1);
    }

    {
      // build the input multivector X with 3 vectors
      const int numVecs  = 3;
      MV X(lclmap,numVecs);
      for (GO i=0; i<static_cast<GO>(P); ++i) {
        for (GO j=0; j<static_cast<GO>(numVecs); ++j) {
          X.replaceGlobalValue(i,j,static_cast<Scalar>(i+j*P));
        }
      }
      // allocate output multivec
      MV Bout(rowmap,numVecs);
      // test the action
      Bout.randomize();
      Tpetra::Details::KokkosRegionCounter::reset();
      Tpetra::Details::KokkosRegionCounter::start();
      A.apply(X,Bout);
      Tpetra::Details::KokkosRegionCounter::stop();

#if defined(HAVE_TPETRA_CUDA) || defined(HAVE_TPETRACORE_CUDA)

      // We don't use cuSPARSE SpMM for all versions.
      // Below is the exact condition used to enable cuSPARSE for spmv_mv in KokkosKernels.
      // (see KokkosSparse_spmv_mv_tpl_spec_avail.hpp for more details)
      //
      // - 10300 and below produced incorrect results and failed KK tests
      // - 11702 also produced incorrect results for very small inputs, causing a Tpetra test to fail
#if defined(CUSPARSE_VERSION) && (10301 <= CUSPARSE_VERSION) && (CUSPARSE_VERSION != 11702)
      if constexpr (std::is_same_v<typename Node::execution_space, Kokkos::Cuda>) {
        const size_t numTplCalls =
          Tpetra::Details::KokkosRegionCounter::get_count_region_contains("spmv[TPL_") +
          Tpetra::Details::KokkosRegionCounter::get_count_region_contains("spmv_mv[TPL_"); // added in Kernels 4.3.1, okay to look for before
        TEST_COMPARE(numTplCalls, ==, 1);
      }
#endif // Specific cuSparse versions work
#endif // CUDA

#if defined(HAVE_TPETRA_HIP) || defined(HAVE_TPETRACORE_HIP)
      // The multivector case is not yet hooked up in Kokkos Kernels.
      if constexpr (std::is_same_v<typename Node::execution_space, Kokkos::HIP>) {
        const size_t numTplCalls =
          Tpetra::Details::KokkosRegionCounter::get_count_region_contains("spmv[TPL_") +
          Tpetra::Details::KokkosRegionCounter::get_count_region_contains("spmv_mv[TPL_"); // added in Kernels 4.3.1, okay to look for before
#if (KOKKOSKERNELS_VERSION >= 40399)
        // rocSparse spmv_mv has been added to KokkosKernels develop
        // Delete the #else branch at the next release 
        TEST_COMPARE(numTplCalls, ==, 1);
#else
        TEST_COMPARE(numTplCalls, ==, 0);
#endif
      }
#endif // HIP
    }

    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0; // output argument
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "KokkosKernels TPL use was not detected where it was expected!" << endl;
    }
  }

//
// INSTANTIATIONS
//
#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, NonSquare, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )
}
