// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_Core.hpp>

#include "Amesos2_MultiVecAdapter.hpp"
#include "Amesos2_Util.hpp"
#include "Amesos2_Meta.hpp"

namespace {

  using std::cout;
  using std::endl;
  using std::string;

  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::ptrInArg;
  using Teuchos::outArg;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::ScalarTraits;
  using Teuchos::OrdinalTraits;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;

  using Tpetra::global_size_t;
  using Tpetra::MultiVector;
  using Tpetra::Map;

  using Amesos2::MultiVecAdapter;
  using Amesos2::createMultiVecAdapter;

  using Amesos2::Util::get_1d_copy_helper;
  using Amesos2::Util::put_1d_data_helper;
  using Amesos2::ROOTED;
  using Amesos2::DISTRIBUTED;
  using Amesos2::GLOBALLY_REPLICATED;


  typedef Tpetra::Map<>::node_type Node;

  bool testMpi = false;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-mpi", "test-serial", &testMpi,
                  "Test Serial by default (for now) or force MPI test.  In a serial build,"
                  " this option is ignored and a serial comm is always used." );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    RCP<const Comm<int> > ret;
    if( testMpi ){
      ret = Tpetra::getDefaultComm();
    } else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  /*
  RCP<FancyOStream> getDefaultOStream()
  {
    return( VerboseObjectBase::getDefaultOStream() );
  }
  */

  /*
   * UNIT TESTS
   */

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Initialization, SCALAR, LO, GO )
  {
    /* Test correct initialization of the MultiVecAdapter
     *
     * - All Constructors
     * - Correct initialization of class members
     * - Correct typedefs
     */
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = getDefaultComm();
    // const size_t numprocs = comm->getSize();
    // const size_t rank     = comm->getRank();
    // create a Map
    const size_t numLocal = 10;

    out << "Creating Map" << std::endl;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

    out << "Creating MultiVector" << std::endl;
    RCP<MV> mv = rcp(new MV(map,11));
    mv->randomize();
    // RCP<FancyOStream> os = getDefaultOStream();
    // mv->describe(*os,Teuchos::VERB_EXTREME);

    out << "Creating adapter" << std::endl;
    RCP<ADAPT> adapter = createMultiVecAdapter(mv);

    // The following should all pass at compile time
    out << "Running compile-time tests" << std::endl;
    TEST_ASSERT( (std::is_same_v<SCALAR,        typename ADAPT::scalar_t>) );
    TEST_ASSERT( (std::is_same_v<LO,            typename ADAPT::local_ordinal_t>) );
    TEST_ASSERT( (std::is_same_v<GO,            typename ADAPT::global_ordinal_t>) );
    TEST_ASSERT( (std::is_same_v<Node,          typename ADAPT::node_t>) );
    TEST_ASSERT( (std::is_same_v<global_size_t, typename ADAPT::global_size_t>) );
    TEST_ASSERT( (std::is_same_v<MV,            typename ADAPT::multivec_t>) );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Dimensions, SCALAR, LO, GO )
  {
    // Test that the dimensions reported by the adapter match those as reported
    // by the Tpetra::MultiVector
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = getDefaultComm();
    // const size_t numprocs = comm->getSize();
    // const size_t rank     = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

    RCP<MV> mv = rcp(new MV(map,11));
    mv->randomize();
    // RCP<FancyOStream> os = getDefaultOStream();
    // mv->describe(*os,Teuchos::VERB_EXTREME);

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);

    TEST_EQUALITY( mv->getLocalLength(),  adapter->getLocalLength()      );
    TEST_EQUALITY( mv->getNumVectors(),   adapter->getLocalNumVectors()  );
    TEST_EQUALITY( mv->getNumVectors(),   adapter->getGlobalNumVectors() );
    TEST_EQUALITY( mv->getGlobalLength(), adapter->getGlobalLength()     );
    TEST_EQUALITY( mv->getStride(),       adapter->getStride()           );

  }

  // Test the get1dCopy() method of MultiVecAdapter.  Check against a
  // known multivector, and also check against what is returned by the
  // Tpetra::MultiVector.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Copy, SCALAR, LO, GO )
  {
    using std::endl;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

    out << "Getting communicator" << endl;

    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numprocs = comm->getSize();
    const size_t rank     = comm->getRank();

    // create a Map
    const size_t numVectors = 7;
    const size_t numLocal = 13;
    out << "Creating Map" << endl;
    RCP<Map<LO,GO,Node> > map = rcp (new Map<LO,GO,Node> (INVALID,numLocal,0,comm));

    out << "Creating MultiVector" << endl;
    RCP<MV> mv = rcp (new MV (map, numVectors));
    mv->randomize ();

    out << "Creating adapter" << endl;
    RCP<ADAPT> adapter = createMultiVecAdapter (mv);
    TEST_EQUALITY_CONST( adapter.is_null (), false );

    Array<SCALAR> original (numVectors * numLocal * numprocs);
    Array<SCALAR> copy (numVectors * numLocal * numprocs);

    ///////////////////////////////////
    // Check a global copy at rank=0 //
    ///////////////////////////////////

    out << "Checking global copy" << endl;
    get_1d_copy_helper<ADAPT,SCALAR>::do_get (ptrInArg (*adapter), copy (),
                                              numLocal*numprocs, ROOTED);

    // Only rank=0 process has global copy of the mv data, check against an import
    size_t my_num_elems = OrdinalTraits<size_t>::zero();
    if (rank == 0) {
      my_num_elems = numLocal*numprocs;
    }
    RCP<const Map<LO,GO,Node> > root_map =
      rcp (new Map<LO,GO,Node> (numLocal*numprocs, my_num_elems, 0, comm));
    MV root_mv (root_map, numVectors);
    Tpetra::Import<LO,GO,Node> importer (map, root_map);
    root_mv.doImport (*mv, importer, Tpetra::REPLACE);

    root_mv.get1dCopy (original (), numLocal * numprocs);

    TEST_COMPARE_ARRAYS( original, copy );

    /////////////////////////////////////
    // Check a local copy at all ranks //
    /////////////////////////////////////

    out << "Checking local copy" << endl;

    TEST_EQUALITY_CONST( mv.is_null (), false );
    TEST_EQUALITY_CONST( adapter.is_null (), false );

    original.clear();
    original.resize(numVectors*numLocal);
    copy.clear();
    copy.resize(numVectors*numLocal);
    mv->randomize();

    out << "About to call get1dCopy" << endl;
    // FIXME (mfh 22 Jan 2014) Should use local stride, not local length.
    mv->get1dCopy (original (), mv->getLocalLength ());
    out << "Successfully called get1dCopy" << endl;

    get_1d_copy_helper<ADAPT,SCALAR>::do_get(ptrInArg(*adapter), copy(), numLocal, DISTRIBUTED);

    out << "Successfully called get_1d_copy_helper" << endl;

    // Check that the values remain the same
    TEST_COMPARE_ARRAYS( original, copy );
    out << "Done!" << endl;
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Copy_Map, SCALAR, LO, GO )
  {
    /* Test the get1dCopy() method of MultiVecAdapter.  This test
       checks the map-based interface to the get1dCopy method.  We
       create an initial multivector which is distributed across all
       MPI processes and then get a copy of the multivector data on
       only the first 2 ranks.
     */
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numprocs = comm->getSize();
    const size_t rank     = comm->getRank();

    // create a Map
    const size_t numVectors = 7;
    const size_t numLocal = 13;
    const size_t total_rows = numLocal * numprocs;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

    RCP<MV> mv = rcp(new MV(map,numVectors));
    mv->randomize();

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);
    Array<SCALAR> original(numVectors*numLocal*numprocs);
    Array<SCALAR> global_copy(numVectors*numLocal*numprocs);

    /*
     * We divide the numLocal * numprocs rows of the multivector
     * amongst the first two rank (or just the first rank in case we
     * don't have more than one) and get the copies.  Each processor
     * checks against a globally-replicated copy of the multivector
     * data.
     */
    size_t my_num_rows = OrdinalTraits<size_t>::zero();
    if ( numprocs > 1 ){
      if ( rank < 2 ){
        my_num_rows = total_rows / 2;
      }
      // If we have an odd number of rows, rank=0 gets the remainder
      if ( rank == 0 ) my_num_rows += total_rows % 2;
    } else {
      my_num_rows = total_rows;
    }
    const Tpetra::Map<LO,GO,Node> redist_map(total_rows, my_num_rows, 0, comm);

    // Get first the global data copy
    get_1d_copy_helper<ADAPT,SCALAR>::do_get(ptrInArg(*adapter),
                                             global_copy(),
                                             total_rows,
                                             GLOBALLY_REPLICATED);

    // Now get a copy using the map
    Array<SCALAR> my_copy(numVectors * my_num_rows);
    get_1d_copy_helper<ADAPT,SCALAR>::do_get(ptrInArg(*adapter),
                                             my_copy(),
                                             my_num_rows,
                                             ptrInArg(redist_map));

    // Check that you have the data you wanted.  Data is stored in
    // column-major order.
    if ( numprocs > 1 ){
      if ( rank == 0 ){
        // Should get the first ceil(total_rows/2) rows
        size_t vec_num = 0;
        size_t vec_ind = 0;
        for( ; vec_ind < numVectors; ++vec_ind ){
          for( size_t i = 0; i < my_num_rows; ++i ){
            SCALAR mv_value = global_copy[total_rows*vec_num + i];
            SCALAR my_value = my_copy[my_num_rows*vec_num + i];
            TEST_EQUALITY( mv_value, my_value );
          }
        }
      } else if ( rank == 1 ){
        // Should get the last floor(total_rows/2) rows
        size_t vec_num = 0;
        size_t vec_ind = 0;
        for( ; vec_ind < numVectors; ++vec_ind ){
          // Iterate backwards through rows
          for( size_t i = total_rows-1; i >= total_rows - my_num_rows; --i ){
            SCALAR mv_value = global_copy[total_rows*vec_num + i];
            SCALAR my_value = my_copy[my_num_rows*vec_num + i];
            TEST_EQUALITY( mv_value, my_value );
          }
        }
      }
    } else {
      // Otherwise, rank=0 should have gotten the whole thing
      TEST_COMPARE_ARRAYS( global_copy, my_copy );
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( MultiVecAdapter, Globalize, SCALAR, LO, GO )
  {
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef MultiVecAdapter<MV> ADAPT;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();

    RCP<const Comm<int> > comm = getDefaultComm();
    const size_t numprocs = comm->getSize();
    const size_t rank     = comm->getRank();
    // create a Map
    const size_t numVectors = 7;
    const size_t numLocal = 13;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );

    RCP<MV> mv = rcp(new MV(map,numVectors));
    mv->randomize();

    RCP<ADAPT> adapter = createMultiVecAdapter(mv);
    Array<SCALAR> original(numVectors*numLocal*numprocs);
    Array<SCALAR> copy(numVectors*numLocal*numprocs);

    if( rank == 0 ){
      std::fill(original.begin(), original.end(), (SCALAR) 1.9);
    }

    // distribute rank 0's data
    put_1d_data_helper<ADAPT,SCALAR>::do_put(outArg(*adapter), original(),
                                             numLocal*numprocs,
                                             ROOTED);

    // Send rank 0's array to everyone else
    Teuchos::broadcast(*comm, 0, original());

    // Now have everyone get a copy from the multivector adapter
    get_1d_copy_helper<ADAPT,SCALAR>::do_get(ptrInArg(*adapter), copy(),
                                             numLocal*numprocs,
                                             GLOBALLY_REPLICATED);

    TEST_EQUALITY( original, copy );
  }


  /*
   * Instantiations
   */

#ifdef HAVE_TEUCHOS_COMPLEX
#  ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO) \
  typedef std::complex<float> ComplexFloat;             \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexFloat)
#  else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  endif
#  ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)        \
  typedef std::complex<double> ComplexDouble;                   \
  UNIT_TEST_GROUP_ORDINAL_SCALAR(LO, GO, ComplexDouble)
#  else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#  endif
#else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#endif

#ifdef HAVE_TPETRA_INST_FLOAT
#  define UNIT_TEST_GROUP_ORDINAL_FLOAT( LO, GO )       \
  UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, float )
#else
#  define UNIT_TEST_GROUP_ORDINAL_FLOAT( LO, GO )
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
#  define UNIT_TEST_GROUP_ORDINAL_DOUBLE( LO, GO )      \
  UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, double )
#else
#  define UNIT_TEST_GROUP_ORDINAL_DOUBLE( LO, GO )
#endif

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, SCALAR )                \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Initialization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Dimensions, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Copy, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Copy_Map, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( MultiVecAdapter, Globalize, SCALAR, LO, GO ) \

#  define UNIT_TEST_GROUP_ORDINAL( ORDINAL )            \
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

#ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#  define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO )     \
  UNIT_TEST_GROUP_ORDINAL_DOUBLE( LO, GO)               \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT( LO, GO )

  UNIT_TEST_GROUP_ORDINAL(int)

#else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#  define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO )     \
  UNIT_TEST_GROUP_ORDINAL_FLOAT(LO, GO)                 \
  UNIT_TEST_GROUP_ORDINAL_DOUBLE(LO, GO)                \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)         \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO,GO)

  //Modified by JDB (Oct-19-2015)
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
  UNIT_TEST_GROUP_ORDINAL(int)
  typedef long int LongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongInt )
  #ifdef HAVE_TPETRA_INT_LONG_LONG
  typedef long long int LongLongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongLongInt )
  #endif
#else  //ETI
  #ifdef HAVE_TPETRA_INST_INT_INT
  UNIT_TEST_GROUP_ORDINAL(int)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG
  typedef long int LongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL(int,LongInt)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
  typedef long long int LongLongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL(int,LongLongInt)
  #endif
#endif  // EXPL-INST

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

} // end anonymous namespace
