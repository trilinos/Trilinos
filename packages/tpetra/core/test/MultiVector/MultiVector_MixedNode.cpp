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
#include "Tpetra_Vector.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include <iterator>

// FINISH: add test for MultiVector with a node containing zero local entries
// FINISH: add tests for local MultiVectors


// Macro that marks a function as "possibly unused," in order to
// suppress build warnings.
#if ! defined(TRILINOS_UNUSED_FUNCTION)
#  if defined(__GNUC__) || (defined(__INTEL_COMPILER) && !defined(_MSC_VER))
#    define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#  elif defined(__clang__)
#    if __has_attribute(unused)
#      define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#    else
#      define TRILINOS_UNUSED_FUNCTION
#    endif // Clang has 'unused' attribute
#  elif defined(__IBMCPP__)
// IBM's C++ compiler for Blue Gene/Q (V12.1) implements 'used' but not 'unused'.
//
// http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp
#    define TRILINOS_UNUSED_FUNCTION
#  else // some other compiler
#    define TRILINOS_UNUSED_FUNCTION
#  endif
#endif // ! defined(TRILINOS_UNUSED_FUNCTION)


namespace Teuchos {
  template <>
  ScalarTraits<int>::magnitudeType
  relErr( const int &s1, const int &s2 ) {
    typedef ScalarTraits<int> ST;
    return ST::magnitude(s1-s2);
  }

  template <>
  ScalarTraits<char>::magnitudeType
  relErr( const char &s1, const char &s2 ) {
    typedef ScalarTraits<char> ST;
    return ST::magnitude(s1-s2);
  }
}

namespace {

  using Tpetra::TestingUtilities::getDefaultComm;

  using std::endl;
  using std::copy;
  using std::string;

  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::arrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::null;
  using Teuchos::Range1D;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::OrdinalTraits;
  using Teuchos::outArg;
  using Teuchos::ScalarTraits;
  using Teuchos::SerialDenseMatrix;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  using Tpetra::Map;
  using Tpetra::MultiVector;
  using Tpetra::global_size_t;
  using Tpetra::GloballyDistributed;

  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;

  double errorTolSlack = 1.0e+2;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
    clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
  }

  // no ScalarTraits<>::eps() for integer types
  template <class Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType testingTol() { return Teuchos::ScalarTraits<Scalar>::eps(); }
  template <>
  TRILINOS_UNUSED_FUNCTION int testingTol<int>() { return 0; }
  template <>
  TRILINOS_UNUSED_FUNCTION long testingTol<long>() { return 0; }

  //
  // UNIT TESTS
  //

  template<class Scalar,class MapClass>
  bool mv_test_basic(RCP<const MapClass> & map, size_t numVecs, Teuchos::FancyOStream & out) {
    bool success=true;
    using LO = typename MapClass::local_ordinal_type;
    using GO = typename MapClass::global_ordinal_type;
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    using MV = MultiVector<Scalar,LO,GO,typename MapClass::node_type>;

   Teuchos::RCP<MV> mvec;
    TEST_NOTHROW( mvec = rcp (new MV(map, numVecs, true)) );
    if (mvec.is_null ()) {
      out << "MV constructor threw an exception: returning" << endl;
      return false;
    }
    TEST_EQUALITY( mvec->getNumVectors(), numVecs );
    TEST_EQUALITY( mvec->getLocalLength(), map->getLocalNumElements() );
    TEST_EQUALITY( mvec->getGlobalLength(), map->getGlobalNumElements() );
    
    out << "Test that all norms are zero" << endl;
    Teuchos::Array<Magnitude> norms(numVecs), zeros(numVecs);
    std::fill(zeros.begin(),zeros.end(),ScalarTraits<Magnitude>::zero());
    TEST_NOTHROW( mvec->norm2(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    TEST_NOTHROW( mvec->norm1(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    TEST_NOTHROW( mvec->normInf(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,zeros,ScalarTraits<Magnitude>::zero());
    out << *mvec << endl;

    Magnitude onesNorm = map->getGlobalNumElements() * ScalarTraits<Magnitude>::one();
    Teuchos::Array<Magnitude> checks(numVecs);      
    std::fill(checks.begin(),checks.end(),onesNorm);
    out << "Fill with all ones and check norm" <<endl;
    mvec->putScalar(ScalarTraits<Scalar>::one());
    TEST_NOTHROW( mvec->norm1(norms) );
    TEST_COMPARE_FLOATING_ARRAYS(norms,checks,ScalarTraits<Magnitude>::zero());
    out << *mvec << endl;
    return success;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVector, mixed_node_basic, Scalar, LO, GO)
  {
    // Here we use CudaNode on Rank 0 and SerialNode on everything else
    typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
    constexpr bool debug = true;
    using SerialNode         = Tpetra::KokkosCompat::KokkosSerialWrapperNode;
    using SerialMap          = Map<LO,GO,SerialNode>;
    using CudaNode           = Tpetra::KokkosCompat::KokkosCudaWrapperNode;
    using CudaMap            = Map<LO,GO,CudaNode>;

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test: MultiVector, mixed_node_basic" << endl;
    Teuchos::OSTab tab0 (myOut);

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid ();
    RCP<const Comm<int> > comm = getDefaultComm ();
    const int numImages = comm->getSize ();
    const int myRank = comm->getRank();

    myOut << "Create Map" << endl;
    const size_t numLocal = 13;
    const size_t numVecs  = 7;
    const GO indexBase = 0;
    RCP<const CudaMap> map0;
    RCP<const SerialMap> map1;
    if(myRank==0) map0=rcp (new CudaMap (INVALID, numLocal, indexBase, comm));
    else map1=rcp (new SerialMap (INVALID, numLocal, indexBase, comm));

    myOut << "Test MultiVector's usual constructor" << endl;  
    if(myRank==0) success=mv_test_basic<Scalar,CudaMap>(map0,numVecs,myOut);
    else success=mv_test_basic<Scalar,SerialMap>(map1,numVecs,myOut);

    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }


  ///
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVector, mixed_node_basic, SCALAR, LO, GO ) 


  typedef Tpetra::Map<>::local_ordinal_type default_local_ordinal_type;
  typedef Tpetra::Map<>::global_ordinal_type default_global_ordinal_type;


  TPETRA_ETI_MANGLING_TYPEDEFS()


  TPETRA_INSTANTIATE_SLG_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

}
