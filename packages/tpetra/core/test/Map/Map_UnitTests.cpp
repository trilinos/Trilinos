// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Map.hpp"
#include "Tpetra_TestingUtilities.hpp"
#include <type_traits> // std::is_same

namespace {

  using Tpetra::createUniformContigMapWithNode;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::outArg;
  using Teuchos::tuple;
  using Teuchos::OrdinalTraits;
  using std::endl;

  using GST = Tpetra::global_size_t;

#define TEST_IS_COMPATIBLE(m1,m2,is_compat)               \
{                                                         \
  Teuchos::OSTab tabCompat0 (out);                        \
  out << "Expect " << (is_compat ? "" : "NOT ") << "compatible" << std::endl; \
  Teuchos::OSTab tabCompat1 (out); \
  out << "Is m1 compatible with itself?" << std::endl;  \
  TEST_EQUALITY_CONST(m1.isCompatible(m1), true);       \
  out << "Is m2 compatible with itself?" << std::endl;  \
  TEST_EQUALITY_CONST(m2.isCompatible(m2), true);       \
  out << "Is m1 compatible with m2?" << std::endl;      \
  TEST_EQUALITY_CONST(m1.isCompatible(m2), is_compat);  \
  out << "Is m2 compatible with m1?" << std::endl;      \
  TEST_EQUALITY_CONST(m2.isCompatible(m1), is_compat);  \
}

#define TEST_IS_SAME_AS(m1,m2,is_sameas)               \
{                                                      \
  out << "Expect " << (is_sameas ? "" : "NOT ") << "same" << std::endl; \
    TEST_EQUALITY_CONST(m1.isSameAs(m1), true);        \
    TEST_EQUALITY_CONST(m2.isSameAs(m2), true);        \
    TEST_EQUALITY_CONST(m1.isSameAs(m2), is_sameas);   \
    TEST_EQUALITY_CONST(m2.isSameAs(m1), is_sameas);   \
}

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

  //
  // UNIT TESTS
  //

  // This test exercises unsigned GlobalOrdinal.  We don't intend to
  // support unsigned GlobalOrdinal going forward.  However, the test
  // should work regardless of the GlobalOrdinal type.
  TEUCHOS_UNIT_TEST( Map, RogersUnsignedGOBugVerification )
  {
    using LO = Tpetra::Map<>::local_ordinal_type;
#if defined(HAVE_TPETRA_INST_INT_UNSIGNED)
    using GO = unsigned int;
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED_LONG)    
    using GO = unsigned long;
#else
    using GO = Tpetra::Map<>::global_ordinal_type;
#endif
    using map_type = Tpetra::Map<LO, GO>;
    
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    if (numImages < 2) return;
    const int myImageID = comm->getRank();
    const GST GSTI = OrdinalTraits<GST>::invalid();
    RCP<map_type> m;
    TEST_NOTHROW( m = rcp(new map_type(GSTI, tuple<GO>(myImageID), 0, comm)) );
    TEST_ASSERT( m.get () != nullptr );
    if (m.get () != nullptr) {
      TEST_EQUALITY( m->getMinAllGlobalIndex(),
		     static_cast<GO> (0) );
      TEST_EQUALITY( m->getMaxAllGlobalIndex(),
		     static_cast<GO> (numImages-1) );
    }

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, compatibilityTests, LO, GO )
  {
    using std::endl;
    typedef Tpetra::Map<LO,GO> M;

    out << "Test: Map, compatibilityTests" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const GST GSTI = OrdinalTraits<GST>::invalid();

    // test isCompatible()
    // m1.isCompatible(m2) should be true if m1 and m2 have the same number of global entries and the same number of local entries on
    // corresponding nodes
    // test the following scenarios:
    // * same number of global and local entries on all nodes
    // * same number of global entries, but different number of local entries on every node
    // * same number of global entries, but different number of local entries on some nodes
    // * different number of global entries, different number of local entries
    //
    // for each, also:
    // test symmetry   : m1.isCompatible(m2) <=> m2.isCompatible(m1)
    // test reflexivity: m1.isCompatible(m1), m2.isCompatible(m2)
    {
      out << "Contiguous nonuniform ctor, same local GID count globally" << endl;
      M m1(GSTI,myImageID,0,comm),
        m2(GSTI,myImageID,0,comm);
      TEST_IS_COMPATIBLE( m1, m2, true );
    }
    {
      out << "Contiguous nonuniform ctor, different local and global GID counts" << endl;
      M m1(GSTI,myImageID+1,0,comm),
        m2(GSTI,myImageID,0,comm);
      TEST_IS_COMPATIBLE( m1, m2, false);
    }
    if (numImages > 1) {
      // want different num local on every proc; map1:numLocal==[0,...,numImages-1], map2:numLocal==[1,...,numImages-1,0]
      {
        out << "Contiguous nonuniform ctor, same global GID count, "
            << "different local GID counts" << endl;
        M m1(GSTI,myImageID,0,comm),
          m2(GSTI,(myImageID+1)%numImages,0,comm);
        out << "myImageID = " << myImageID
            << ", (myImageID+1) % numImages = "
            << ((myImageID+1) % numImages) << endl;
        Teuchos::OSTab tab1 (out);
        out << "m1.getGlobalNumElements() = " << m1.getGlobalNumElements () << endl
            << "m2.getGlobalNumElements() = " << m2.getGlobalNumElements () << endl
            << "m1.getLocalNumElements() = " << m1.getLocalNumElements () << endl
            << "m2.getLocalNumElements() = " << m2.getLocalNumElements () << endl;
        TEST_IS_COMPATIBLE( m1, m2, false);
      }
      if (numImages > 2) {
        // want different num local on a subset of procs
        // image 0 and numImages-1 get map1:numLocal==[0,numImages-1] and map2:numLocal==[numImages-1,0], the others get numLocal==myImageID
        LO mynl1, mynl2;
        if (myImageID == 0) {
          mynl1 = 0;
          mynl2 = numImages-1;
        }
        else if (myImageID == numImages-1) {
          mynl1 = numImages-1;
          mynl2 = 0;
        }
        else {
          mynl1 = mynl2 = myImageID;
        }
        {
          out << "Contiguous nonuniform ctor, same global GID count, "
              << "different local GID counts on subset of processes" << endl;
          M m1(GSTI,mynl1,0,comm),
            m2(GSTI,mynl2,0,comm);
          TEST_IS_COMPATIBLE( m1, m2, false);
        }
      }
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, sameasTests, LO, GO )
  {
    typedef Tpetra::Map<LO,GO> M;

    out << "Test: Map, sameasTests" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const GST GSTI = OrdinalTraits<GST>::invalid();
    {
      M m1(GSTI,0,0,comm),
        m2(GSTI,0,0,comm);
      TEST_IS_SAME_AS(m1, m2, true);
    }
    {
      M m1(GSTI,myImageID,0,comm),
        m2(GSTI,myImageID,0,comm);
      TEST_IS_SAME_AS(m1, m2, true);
    }
    {
      M m1(GSTI,myImageID,0,comm),
        m2(GSTI,myImageID+1,0,comm);
      TEST_IS_SAME_AS(m1, m2, false);
    }
    if (numImages > 1) {
      // FINISH: test all multi-node scenarios, esp. divergent paths
      {
        M m1(GSTI,myImageID,0,comm),
          m2(GSTI,myImageID+(myImageID==1?1:0),0,comm);
        TEST_IS_SAME_AS(m1, m2, false);
      }
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, ContigUniformMap, LO, GO )
  {
    typedef Tpetra::Map<LO,GO> M;

    out << "Test: Map, ContigUniformMap" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with two entries per node
    // this map will have the following entries:
    Array<GO> myGlobal( tuple<GO>(myImageID*2, myImageID*2+1) );
    Array<LO>  myLocal( tuple<LO>(0,1) );

    const size_t numGlobalEntries = numImages*2;
    const GO indexBase = 0;
    const LO localIndexBase = 0;
    M map(numGlobalEntries,indexBase,comm);

    TEST_EQUALITY_CONST(map.isContiguous(), true);
    TEST_EQUALITY_CONST(map.isDistributed(), numImages > 1);
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobalEntries);
    TEST_EQUALITY_CONST(map.getLocalNumElements(), 2);
    TEST_EQUALITY_CONST(map.getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), localIndexBase);
    TEST_EQUALITY_CONST(map.getMaxLocalIndex(), 1);
    TEST_EQUALITY_CONST(map.getMinGlobalIndex(), myGlobal[0]);
    TEST_EQUALITY_CONST(map.getMaxGlobalIndex(), myGlobal[1]);
    TEST_EQUALITY_CONST(map.getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY_CONST(map.getMaxAllGlobalIndex(), as<GO>(numGlobalEntries-1));
    TEST_EQUALITY( map.getLocalElement(myGlobal[0]), myLocal[0] );
    TEST_EQUALITY( map.getLocalElement(myGlobal[1]), myLocal[1] );
    TEST_EQUALITY( map.getGlobalElement(myLocal[0]), myGlobal[0] );
    TEST_EQUALITY( map.getGlobalElement(myLocal[1]), myGlobal[1] );
    TEST_EQUALITY( map.getLocalElement(numGlobalEntries), OrdinalTraits<LO>::invalid() );
    TEST_EQUALITY( map.getGlobalElement(2),               OrdinalTraits<GO>::invalid() );
    TEST_EQUALITY( map.getLocalElement(numGlobalEntries-1), myImageID == numImages-1 ? 1 : OrdinalTraits<LO>::invalid() );
    TEST_COMPARE_ARRAYS( map.getLocalElementList(), myGlobal);
    TEST_EQUALITY_CONST( map.isNodeLocalElement(0), true );
    TEST_EQUALITY_CONST( map.isNodeLocalElement(1), true );
    TEST_EQUALITY_CONST( map.isNodeLocalElement(2), false ); // just try a couple
    TEST_EQUALITY_CONST( map.isNodeLocalElement(3), false );
    for (GO i=0; i < as<GO>(numGlobalEntries); ++i) {
      if (std::find(myGlobal.begin(),myGlobal.end(),i) == myGlobal.end()) {
        TEST_EQUALITY_CONST( map.isNodeGlobalElement(i), false );
      }
      else {
        TEST_EQUALITY_CONST( map.isNodeGlobalElement(i), true );
      }
    }

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, nonTrivialIndexBase, LO, GO )
  {
    typedef Tpetra::Map<LO,GO> Map;

    out << "Test: Map, nonTrivialIndexBase" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per node
    const size_t        numLocal  = 5;
    const GST numGlobal = numImages*numLocal;
    const GO indexBase = 10;

    Map map (numGlobal, indexBase, comm);
    //
    Array<GO> expectedGIDs(numLocal);
    expectedGIDs[0] = indexBase + myImageID*numLocal;
    for (size_t i=1; i<numLocal; ++i) {
      expectedGIDs[i] = expectedGIDs[i-1]+1;
    }
    //
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map.getLocalNumElements(), numLocal);
    TEST_EQUALITY(map.getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), 0);
    TEST_EQUALITY(map.getMaxLocalIndex(), numLocal-1);
    TEST_EQUALITY(map.getMinGlobalIndex(), expectedGIDs[0]);
    TEST_EQUALITY(map.getMaxGlobalIndex(), as<GO>(expectedGIDs[0]+numLocal-1) );
    TEST_EQUALITY(map.getMinAllGlobalIndex(), indexBase);
    TEST_EQUALITY(map.getGlobalElement(0), expectedGIDs[0]);
    TEST_EQUALITY_CONST((GST)map.getMaxAllGlobalIndex(), indexBase+numGlobal-1);
    ArrayView<const GO> glist = map.getLocalElementList();
    TEST_COMPARE_ARRAYS( map.getLocalElementList(), expectedGIDs);

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, indexBaseAndAllMin, LO, GO )
  {
    typedef Tpetra::Map<LO,GO> Map;

    out << "Test: Map, indexBaseAndAllMin" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per node
    const size_t        numLocal  = 5;
    const GST numGlobal = numImages*numLocal;
    const GO indexBase = 0;
    const GO actualBase = 1;
    //
    Array<GO> GIDs(numLocal);
    GIDs[0] = actualBase + myImageID*numLocal;
    for (size_t i=1; i<numLocal; ++i) {
      GIDs[i] = GIDs[i-1]+1;
    }

    Map map (numGlobal, GIDs (), indexBase, comm);
    //
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map.getLocalNumElements(), numLocal);
    TEST_EQUALITY(map.getIndexBase(), indexBase);
    TEST_EQUALITY_CONST(map.getMinLocalIndex(), 0);
    TEST_EQUALITY(map.getMaxLocalIndex(), numLocal-1);
    TEST_EQUALITY(map.getMinGlobalIndex(), GIDs[0]);
    TEST_EQUALITY(map.getMaxGlobalIndex(), as<GO>(GIDs[0]+numLocal-1) );
    TEST_EQUALITY(map.getMinAllGlobalIndex(), actualBase);
    TEST_EQUALITY(map.getGlobalElement(0), GIDs[0]);
    TEST_EQUALITY_CONST((GST)map.getMaxAllGlobalIndex(), actualBase+numGlobal-1);
    ArrayView<const GO> glist = map.getLocalElementList();
    TEST_COMPARE_ARRAYS( map.getLocalElementList(), GIDs);

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, ZeroLocalElements, LO, GO )
  {
    typedef Tpetra::Map<LO,GO> M;

    out << "Test: Map, ZeroLocalElements" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int rank = comm->getRank();

    // Create maps with zero elements on all but the first processor
    Array<GO>  elem_list;
    if (rank == 0)
      elem_list.push_back(0);
    M contig_uniform(1, 0, comm);
    M contig_non_uniform(1, elem_list.size(), 0, comm);
    M non_contig(1, elem_list, 0, comm);

    // Check LID
    LO lid_expected = rank == 0 ? 0 : OrdinalTraits<LO>::invalid();
    TEST_EQUALITY( contig_uniform.getLocalElement(0), lid_expected );
    TEST_EQUALITY( contig_non_uniform.getLocalElement(0), lid_expected );
    TEST_EQUALITY( non_contig.getLocalElement(0), lid_expected );

    // Make sure that the test passed on all MPI processes.
    const int lclSuccess = success ? 1 : 0;
    int gblSuccess = 0;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
  }


  // Test that Map can be declared with no template parameters, so
  // that every template parameter has its default value.
  TEUCHOS_UNIT_TEST( Map, AllDefaultTemplateParameters )
  {
    // If you are letting all template parameters take their default
    // values, you must follow the class name Map with <>.
    typedef Tpetra::Map<> map_type;
    typedef map_type::local_ordinal_type local_ordinal_type;
    typedef map_type::global_ordinal_type global_ordinal_type;

    out << "Test: Map, AllDefaultTemplateParameters" << std::endl;
    Teuchos::OSTab tab0 (out);

    // Verify that the default LocalOrdinal type is int.  We can't put
    // the is_same expression in the macro, since it has a comma
    // (commas separate arguments in a macro).
    static_assert (std::is_same<local_ordinal_type, int>::value,
		   "The default local_ordinal_type should be int, "
		   "currently at least.  That may change in the future.");

    // Verify that the default GlobalOrdinal type has size no less
    // than the default LocalOrdinal type.
    static_assert (sizeof (global_ordinal_type) >= sizeof (local_ordinal_type),
		   "The default global_ordinal_type must have size "
		   "no less than that of the default local_ordinal_type.");
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, KokkosViewConstructor, LO, GO )
 {
   typedef Tpetra::Map<LO,GO> M;
   using Teuchos::ArrayView;
   GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
   using MemSpace = typename Tpetra::KokkosClassic::DefaultNode::DefaultNodeType::memory_space;

   out << "Test: Map, KokkosViewConstructor" << std::endl;

   // create a comm                                                                                                                                                                                                                                                          
   auto comm = Tpetra::getDefaultComm();
   const int rank = comm->getRank();

   // Create a dummy map to force Kokkos to initialize
   {
     M dummy(10,0,comm);
   }


   // View in default space
   int N = 10;
   Kokkos::View<GO*,MemSpace> myview("DeviceView",N);
   auto myview_h = Kokkos::create_mirror_view(myview);

   int RANK_BASE = rank*100;
   // Contiguous part
   int nstop = 5;
   for(int i=0; i<nstop; i++)
     myview_h(i) = RANK_BASE + i;

   // Noncontiguous goop
   myview_h(nstop) = RANK_BASE + 5000;
   
   // The rest, is contiguous, except for the last guy
   for(int i=nstop+1; i<N-1; i++)
     myview_h(i) = RANK_BASE + i;
   
   // The last guy
   myview_h(N-1) = RANK_BASE + 5001;
   Kokkos::deep_copy(myview,myview_h);

   // Now do an Arrayview
   Teuchos::ArrayView<GO> myview_av(myview_h.data(),myview_h.extent(0));

   // Create maps
   M map_kokkos(INVALID,myview,0,comm);
   M map_teuchos(INVALID,myview_av,0,comm);

   // Compare the easy stuff
   TEST_EQUALITY(map_kokkos.getGlobalNumElements(), map_teuchos.getGlobalNumElements());
   TEST_EQUALITY(map_kokkos.getLocalNumElements(), map_teuchos.getLocalNumElements());
   TEST_EQUALITY(map_kokkos.getIndexBase(), map_teuchos.getIndexBase());
   TEST_EQUALITY(map_kokkos.getMinLocalIndex(), map_teuchos.getMinLocalIndex());
   TEST_EQUALITY(map_kokkos.getMaxLocalIndex(), map_teuchos.getMaxLocalIndex());
   TEST_EQUALITY(map_kokkos.getMaxGlobalIndex(), map_teuchos.getMaxGlobalIndex());
   TEST_EQUALITY(map_kokkos.getMinGlobalIndex(), map_teuchos.getMinGlobalIndex());
   TEST_EQUALITY(map_kokkos.getMaxAllGlobalIndex(), map_teuchos.getMaxAllGlobalIndex());
   TEST_EQUALITY(map_kokkos.getMinAllGlobalIndex(), map_teuchos.getMinAllGlobalIndex());

   ArrayView<const GO> glist_kokkos = map_kokkos.getLocalElementList();
   ArrayView<const GO> glist_teuchos = map_teuchos.getLocalElementList();
   TEST_COMPARE_ARRAYS( glist_kokkos, glist_teuchos);

   // Compare the harder stuff by tickling the FHT  
   for(LO i=0; i<N; i++) {
     LO lo_k = map_kokkos.getLocalElement(map_kokkos.getGlobalElement(i));
     LO lo_t = map_teuchos.getLocalElement(map_teuchos.getGlobalElement(i));
     TEST_EQUALITY(lo_k,lo_t);
   }     

 }
 TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, KokkosViewConstructor2, LO, GO )
 {
   typedef Tpetra::Map<LO,GO> M;
   using Teuchos::ArrayView;
   GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
   using MemSpace = typename Tpetra::KokkosClassic::DefaultNode::DefaultNodeType::memory_space;

   out << "Test: Map, KokkosViewConstructor2" << std::endl;

   // create a comm
   auto comm = Tpetra::getDefaultComm();

   // Create a dummy map to force Kokkos to initialize
   {
     M dummy(10,0,comm);
   }

   // View in default space.  Here we make sure we have a value *below* the 
   // value in entry zero.  
   int N = 3;
   Kokkos::View<GO*,MemSpace> myview("DeviceView",N);
   auto myview_h = Kokkos::create_mirror_view(myview);
   myview_h(0) = 101;
   myview_h(1) = 99;
   myview_h(2) = 100;
   Kokkos::deep_copy(myview,myview_h);

   // Now do an Arrayview
   Teuchos::ArrayView<GO> myview_av(myview_h.data(),myview_h.extent(0));

   // Create maps
   M map_kokkos(INVALID,myview,0,comm);
   M map_teuchos(INVALID,myview_av,0,comm);

   // Compare the easy stuff
   TEST_EQUALITY(map_kokkos.getGlobalNumElements(), map_teuchos.getGlobalNumElements());
   TEST_EQUALITY(map_kokkos.getLocalNumElements(), map_teuchos.getLocalNumElements());
   TEST_EQUALITY(map_kokkos.getIndexBase(), map_teuchos.getIndexBase());
   TEST_EQUALITY(map_kokkos.getMinLocalIndex(), map_teuchos.getMinLocalIndex());
   TEST_EQUALITY(map_kokkos.getMaxLocalIndex(), map_teuchos.getMaxLocalIndex());
   TEST_EQUALITY(map_kokkos.getMaxGlobalIndex(), map_teuchos.getMaxGlobalIndex());
   TEST_EQUALITY(map_kokkos.getMinGlobalIndex(), map_teuchos.getMinGlobalIndex());
   TEST_EQUALITY(map_kokkos.getMaxAllGlobalIndex(), map_teuchos.getMaxAllGlobalIndex());
   TEST_EQUALITY(map_kokkos.getMinAllGlobalIndex(), map_teuchos.getMinAllGlobalIndex());

   ArrayView<const GO> glist_kokkos = map_kokkos.getLocalElementList();
   ArrayView<const GO> glist_teuchos = map_teuchos.getLocalElementList();
   TEST_COMPARE_ARRAYS( glist_kokkos, glist_teuchos);

   // Compare the harder stuff by tickling the FHT  
   for(LO i=0; i<N; i++) {
     LO lo_k = map_kokkos.getLocalElement(map_kokkos.getGlobalElement(i));
     LO lo_t = map_teuchos.getLocalElement(map_teuchos.getGlobalElement(i));
     TEST_EQUALITY(lo_k,lo_t);
   }     

 }



  //
  // INSTANTIATIONS
  //

#ifdef HAVE_TPETRA_DEBUG
  // all ordinals, default node
#  define UNIT_TEST_GROUP( LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, compatibilityTests, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, sameasTests, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, nonTrivialIndexBase, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, indexBaseAndAllMin, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, ContigUniformMap, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, ZeroLocalElements, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, KokkosViewConstructor, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, KokkosViewConstructor2, LO, GO )
#else
  // all ordinals, default node
#  define UNIT_TEST_GROUP( LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, compatibilityTests, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, sameasTests, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, nonTrivialIndexBase, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, indexBaseAndAllMin, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, ContigUniformMap, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, ZeroLocalElements, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, KokkosViewConstructor, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, KokkosViewConstructor2, LO, GO )
#endif // HAVE_TPETRA_DEBUG

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LG(UNIT_TEST_GROUP)

}


