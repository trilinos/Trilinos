// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_as.hpp"

// FINISH: test for createFromRecvs(), that negatives in remoteNodeIDs are met by negatives in exportNodeIDs, and that the placement is
//         is preserved. need to understand the semantics of negatives in the node list.

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Tpetra::Distributor;
  using Teuchos::Array;
  //using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::broadcast;
  using Teuchos::Comm;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  using Teuchos::tuple;
  using std::endl;

  bool testMpi = true;
  double errorTolSlack = 1e+1;
  int numRuns = 1;

  int generateValue(int x, int y) {
    // formula for z(x,y) = 0.5(x^2 + y^2 + 3x + y) + xy
    return(((x*x + y*y + x+x+x + y) / 2) + (x*y));
  }

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
    clp.setOption(
        "num-runs", &numRuns,
        "Number of runs to use for timings" );
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return Tpetra::getDefaultComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  //
  // UNIT TESTS
  //


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsContig)
  {
    RCP<const Comm<int> > comm = getDefaultComm();

    // Set debug = true if you want immediate debug output to stderr.
    const bool debug = false;
    Teuchos::RCP<Teuchos::FancyOStream> outPtr =
      debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Distributor createFromSendsContig" << endl;
    Teuchos::OSTab tab1 (myOut);

    const int numImages = comm->getSize();
    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    size_t numImports = 0;
    // fill exportImageIDs with {0,0, -1,-1, 1,1, 2,2, ... numImages-1,numImages-1}
    // two sends to each image, contiguous, in order
    // note the -1s after 0,0; these should not hurt contiguity or be reflected in numImports
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(numImages*2+2);
    exportImageIDs.push_back(0);
    exportImageIDs.push_back(0);
    exportImageIDs.push_back(-1);
    exportImageIDs.push_back(-1);
    for(int i=1; i < numImages; ++i) {
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
    }

    myOut << "Create Distributor from (contiguous) sends" << endl;
    Distributor distributor(comm);
    numImports = distributor.createFromSends(exportImageIDs);

    myOut << "Test the resulting Distributor" << endl;

    TEST_EQUALITY(numImports, as<size_t>(2*numImages));
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<size_t>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), (numImages > 1 ? 2 : 0))
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<size_t>(2*numImages));

    myOut << "Compare getProcsFrom() and getProcsTo()" << endl;
    {
      ArrayView<const int> imgFrom(distributor.getProcsFrom());
      ArrayView<const int> imgTo(distributor.getProcsTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }

    myOut << "Compare getLengthsFrom() and getLengthsTo()" << endl;
    {
      ArrayView<const size_t> lenFrom = distributor.getLengthsFrom();
      ArrayView<const size_t> lenTo   = distributor.getLengthsTo();
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], 2);
        TEST_EQUALITY_CONST( lenTo[i],   2);
      }
    }

    myOut << "Make sure that Distributor output doesn't cause a hang" << endl;
    distributor.describe (out, Teuchos::VERB_EXTREME);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }

  TEUCHOS_UNIT_TEST( Distributor, badArgsFromSends)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // each node i sends to node i+1
    // for the last node, this results in an invalid node id, which should throw an exception on
    // every node
    size_t numImports = 0;

    const bool debug = Tpetra::Details::Behavior::debug("Distributor");

    // create from sends with bad node IDs
    if (debug) {
      {
        Distributor distributor(comm);
        TEST_THROW( numImports = distributor.createFromSends( tuple<int>(myImageID+1)), std::runtime_error );
        // Printing numImports prevents a compiler warning (set but unused).
        out << "numImports result: " << numImports << std::endl;
      }
      {
        Distributor distributor(comm);
        TEST_THROW( numImports = distributor.createFromSends( tuple<int>(0,myImageID+1,0)), std::runtime_error );
        // Printing numImports prevents a compiler warning (set but unused).
        out << "numImports result: " << numImports << std::endl;
      }
    }
    else {
      out << "Debug mode not enabled; set TPETRA_DEBUG=Distributor "
        "to test." << std::endl;
    }

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }

  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsMixedContig)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    // Set debug = true if you want immediate debug output to stderr.
    const bool debug = false;
    Teuchos::RCP<Teuchos::FancyOStream> outPtr =
      debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Distributor createFromSendsMixedContig" << endl;
    Teuchos::OSTab tab1 (myOut);

    if (numImages < 2) {
      myOut << "comm->getSize() = " << numImages << " < 2.  The test makes no "
        "sense to run in this case, so I'll skip it." << endl;
      return;
    }

    // even is black, odd is red
    bool even = ((myImageID % 2) == 0);
    // two exports to each image, including myself
    // on even imageIDs, send data contig
    // on odd imageIDs, send data non-contig
    size_t numImports = 0;
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(numImages*2);
    if (even) {
      // fill exportImageIDs with {0,0, 1,1, 2,2, ... numImages-1,numImages-1}
      for(int i = 0; i < numImages; ++i) {
        exportImageIDs.push_back(i);
        exportImageIDs.push_back(i);
      }
    }
    else {
      // fill exportImageIDs with {0,1,2,...,numImages-1, 0,1,2,...,numImages-1}
      for (int i = 0; i < numImages; ++i) {
        exportImageIDs.push_back(i);
      }
      for (int i = 0; i < numImages; ++i) {
        exportImageIDs.push_back(i);
      }
    }

    myOut << "Create Distributor from sends (mix of contiguous and noncontiguous)" << endl;

    Distributor distributor(comm);
    TEST_NOTHROW( numImports = distributor.createFromSends(exportImageIDs) );

    myOut << "Test the resulting Distributor" << endl;

    TEST_EQUALITY(numImports, as<size_t>(2*numImages));
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<size_t>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 2);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<size_t>(2*numImages));

    myOut << "Compare getProcsFrom() and getProcsTo()" << endl;
    {
      ArrayView<const int> imgFrom(distributor.getProcsFrom());
      ArrayView<const int> imgTo(distributor.getProcsTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }

    myOut << "Compare getLengthsFrom() and getLengthsTo()" << endl;
    {
      ArrayView<const size_t> lenFrom = distributor.getLengthsFrom();
      ArrayView<const size_t> lenTo   = distributor.getLengthsTo();
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], 2);
        TEST_EQUALITY_CONST( lenTo[i],   2);
      }
    }

    myOut << "Make sure that Distributor output doesn't cause a hang" << endl;
    distributor.describe (out, Teuchos::VERB_EXTREME);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsRedBlack)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    // Set debug = true if you want immediate debug output to stderr.
    const bool debug = false;
    Teuchos::RCP<Teuchos::FancyOStream> outPtr =
      debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Distributor createFromSendsRedBlack" << endl;
    Teuchos::OSTab tab1 (myOut);

    if (numImages < 3) {
      myOut << "comm->getSize() = " << numImages << " < 3.  The test makes no "
        "sense to run in this case, so I'll skip it." << endl;
      return;
    }

    // partition world into red/black (according to imageID even/odd)
    // even is black, odd is red
    bool black = ((myImageID % 2) == 0);
    size_t numInMyPartition = 0;
    size_t numImports = 0;
    // fill exportImageIDs with all images from partition
    Array<int> exportImageIDs(0);
    if (black) {
      // evens
      for(int i=0; i < numImages; i+=2) {
        exportImageIDs.push_back(i);
        numInMyPartition++;
      }
    }
    else {
      // odds
      for(int i=1; i < numImages; i+=2) {
        exportImageIDs.push_back(i);
        numInMyPartition++;
      }
    }

    out << "Create Distributor from sends" << endl;

    // create from contiguous sends
    Distributor distributor(comm);
    numImports = distributor.createFromSends(exportImageIDs);

    out << "Test the resulting Distributor" << endl;

    TEST_EQUALITY(numImports, numInMyPartition);
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), numInMyPartition-1);
    TEST_EQUALITY(distributor.getNumReceives(), numInMyPartition-1);
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), (numInMyPartition > 1 ? 1 : 0));
    TEST_EQUALITY(distributor.getTotalReceiveLength(), numInMyPartition);

    myOut << "Compare getProcsFrom() and getProcsTo()" << endl;
    {
      ArrayView<const int> imgFrom(distributor.getProcsFrom());
      ArrayView<const int> imgTo(distributor.getProcsTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }

    myOut << "Compare getLengthsFrom() and getLengthsTo()" << endl;
    {
      ArrayView<const size_t> lenFrom = distributor.getLengthsFrom();
      ArrayView<const size_t> lenTo   = distributor.getLengthsTo();
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numInMyPartition));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numInMyPartition));
      for (size_t i=0; i<numInMyPartition; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], 1);
        TEST_EQUALITY_CONST( lenTo[i],   1);
      }
    }

    myOut << "Make sure that Distributor output doesn't cause a hang" << endl;
    distributor.describe (out, Teuchos::VERB_EXTREME);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsContigNoself)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    // Set debug = true if you want immediate debug output to stderr.
    const bool debug = false;
    Teuchos::RCP<Teuchos::FancyOStream> outPtr =
      debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Distributor createFromSendsContigNoself" << endl;
    Teuchos::OSTab tab1 (myOut);

    if (numImages < 2) {
      myOut << "comm->getSize() = " << numImages << " < 2.  The test makes no "
        "sense to run in this case, so I'll skip it." << endl;
      return;
    }

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    size_t numImports = 0;
    // fill exportImageIDs with {0,1,...,myImageID-1,myImageID+1,...,numImages-1}
    // one send to each image, contiguous, in order, but not to myself
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(numImages-1);
    for (int i=0; i < myImageID; ++i) {
      exportImageIDs.push_back(i);
    }
    for (int i = myImageID+1; i < numImages; ++i) {
      exportImageIDs.push_back(i);
    }

    myOut << "Create Distributor from sends" << endl;

    // create from contiguous sends
    Distributor distributor(comm);
    numImports = distributor.createFromSends(exportImageIDs);

    myOut << "Test the resulting Distributor" << endl;

    TEST_EQUALITY(numImports, as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), false);
    TEST_EQUALITY(distributor.getNumSends(), as<size_t>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 1);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<size_t>(numImages-1));

    myOut << "Compare getProcsFrom() and getProcsTo()" << endl;
    {
      ArrayView<const int> imgFrom(distributor.getProcsFrom());
      ArrayView<const int> imgTo(distributor.getProcsTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }

    myOut << "Compare getLengthsFrom() and getLengthsTo()" << endl;
    {
      ArrayView<const size_t> lenFrom(distributor.getLengthsFrom());
      ArrayView<const size_t> lenTo(distributor.getLengthsTo());
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages-1));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages-1));
      for (int i=0; i<numImages-1; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], 1);
        TEST_EQUALITY_CONST( lenTo[i],   1);
      }
    }

    myOut << "Make sure that Distributor output doesn't cause a hang" << endl;
    distributor.describe (out, Teuchos::VERB_EXTREME);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsContigUnordered)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // Set debug = true if you want immediate debug output to stderr.
    const bool debug = false;
    Teuchos::RCP<Teuchos::FancyOStream> outPtr =
      debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Distributor createFromSendsContigUnordered" << endl;
    Teuchos::OSTab tab1 (myOut);

    if (numImages < 3) {
      myOut << "comm->getSize() = " << numImages << " < 3.  The test makes no "
        "sense to run in this case, so I'll skip it." << endl;
      return;
    }

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    size_t numImports = 0;
    // fill exportImageIDs with {0,0,0, 1,1,1, 2,2,2, ... numImages-1,numImages-1,numImages-1}
    // three sends to each image, out of order (even first, then odd)
    // only test if numImages > 2
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(numImages*3);
    // even first: {0,0,0, 2,2,2, 4,4,4, ...}
    for(int i=0; i < numImages; i+=2) {
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
    }
    // then odd: {1,1,1, 3,3,3, 5,5,5, ...}
    for(int i=1; i < numImages; i+=2) {
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
      exportImageIDs.push_back(i);
    }

    myOut << "Create Distributor from contiguous sends" << endl;

    Distributor distributor(comm);
    numImports = distributor.createFromSends(exportImageIDs);

    myOut << "Test the resulting Distributor" << endl;

    TEST_EQUALITY(numImports, as<size_t>(3*numImages));
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<size_t>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 3);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<size_t>(3*numImages));

    myOut << "Compare getProcsFrom() and getProcsTo()" << endl;
    {
      ArrayView<const int> imgFrom(distributor.getProcsFrom());
      ArrayView<const int> imgTo(distributor.getProcsTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }

    myOut << "Compare getLengthsFrom() and getLengthsTo()" << endl;
    {
      ArrayView<const size_t> lenFrom = distributor.getLengthsFrom();
      ArrayView<const size_t> lenTo   = distributor.getLengthsTo();
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST(lenFrom[i], 3);
        TEST_EQUALITY_CONST(lenTo[i],   3);
      }
    }

    myOut << "Make sure that Distributor output doesn't cause a hang" << endl;
    distributor.describe (out, Teuchos::VERB_EXTREME);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsNonContig)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

    // Set debug = true if you want immediate debug output to stderr.
    const bool debug = false;
    Teuchos::RCP<Teuchos::FancyOStream> outPtr =
      debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Distributor createFromSendsNonContig" << endl;
    Teuchos::OSTab tab1 (myOut);

    // send data to each image, including myself
    // the consequence is that each image will send to every other images
    size_t numImports = 0;
    // put some -1s in there
    // fill exportImageIDs with {0, 1, 2, ... numImages-1,
    //                           -1,
    //                           0, 1, 2, ... numImages-1}
    Array<int> exportImageIDs(0);
    exportImageIDs.reserve(2*numImages);
    for(int i=0; i < numImages; ++i) {
      exportImageIDs.push_back(i);
    }
    exportImageIDs.push_back(-1);
    for(int i=0; i < numImages; ++i) {
      exportImageIDs.push_back(i);
    }

    myOut << "Create Distributor from noncontiguous sends" << endl;

    Distributor distributor(comm);
    numImports = distributor.createFromSends(exportImageIDs);

    myOut << "Test the resulting Distributor" << endl;

    TEST_EQUALITY(numImports, as<size_t>(2*numImages));
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<size_t>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<size_t>(numImages-1));
    if (numImages == 1) {
      TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 0);
    }
    else {
      TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 2);
    }
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<size_t>(2*numImages));

    myOut << "Compare getProcsFrom() and getProcsTo()" << endl;
    {
      ArrayView<const int> imgFrom(distributor.getProcsFrom());
      ArrayView<const int> imgTo(distributor.getProcsTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }

    myOut << "Compare getLengthsFrom() and getLengthsTo()" << endl;
    {
      ArrayView<const size_t> lenFrom = distributor.getLengthsFrom();
      ArrayView<const size_t> lenTo   = distributor.getLengthsTo();
      TEST_EQUALITY(lenFrom.size(),as<Teuchos_Ordinal>(numImages));
      TEST_EQUALITY(lenTo.size()  ,as<Teuchos_Ordinal>(numImages));
      for (int i=0; i<numImages; ++i) {
        TEST_EQUALITY_CONST( lenFrom[i], 2);
        TEST_EQUALITY_CONST( lenTo[i],   2);
      }
    }

    myOut << "Make sure that Distributor output doesn't cause a hang" << endl;
    distributor.describe (out, Teuchos::VERB_EXTREME);

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, doPostsContig, Packet )
  {
    typedef Teuchos::ScalarTraits<Packet>   PT;
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    // Set debug = true if you want immediate debug output to stderr.
    const bool debug = false;
    Teuchos::RCP<Teuchos::FancyOStream> outPtr =
      debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Distributor doPostsContig" << endl;
    Teuchos::OSTab tab1 (myOut);

    // send data to each image, including myself
    const size_t numExportIDs = numImages;
    size_t numRemoteIDs = 0;
    // fill exportImageIDs with {0, -1, 1, -1, 2, -1, ... numImages-1}
    // on root node only, interlace node IDs with invalid nodes, corresponding to untouched data in import/export buffers
    Array<int> exportImageIDs;
    if (myImageID == 0) {
      exportImageIDs.reserve(2*numExportIDs-1);
      exportImageIDs.push_back(0);
      for(int i=1; i < as<int>(numExportIDs); ++i) {
        exportImageIDs.push_back(-1);
        exportImageIDs.push_back(i);
      }
    }
    else {
      exportImageIDs.reserve(numExportIDs);
      for(int i=0; i < as<int>(numExportIDs); ++i) {
        exportImageIDs.push_back(i);
      }
    }

    myOut << "Create Distributor from sends" << endl;
    Distributor distributor(comm);
    numRemoteIDs = distributor.createFromSends(exportImageIDs);

    myOut << "Make sure that Distributor output doesn't cause a hang" << endl;
    distributor.describe (out, Teuchos::VERB_EXTREME);

    myOut << "Check getReverse(create=false)" << endl;
    RCP<Distributor> revDistor = distributor.getReverse(false);
    TEUCHOS_ASSERT(revDistor.is_null());

    myOut << "Check getReverse(create=true)" << endl;
    revDistor = distributor.getReverse();
    TEUCHOS_ASSERT(!revDistor.is_null());

    myOut << "Check return value of createFromSends" << endl;
    TEST_EQUALITY(numRemoteIDs, as<size_t>(numImages));

    myOut << "Generate data set for doPosts on Process 0" << endl;
    // generate global random data set: each image sends 1 packet to each image
    // we need numImages*numImages "unique" values (we don't want redundant data allowing false positives)
    // root node generates all values, sends them to the others.
    Array<Packet> exports(numImages*numImages);
    if (myImageID == 0) {
      for (int i=0; i<numImages*numImages; i++) {
        exports[i] = PT::random();
      }
    }
    myOut << "Broadcast data set from Process 0" << endl;
    broadcast(*comm,0,exports());

    // pick a subset of entries to post
    Array<Packet> myExports(0);
    if (myImageID == 0) {
      myExports.resize(2*numImages-1,0);
      for (int i=0; i<numImages; ++i) {
        myExports[2*i] = exports[i];
      }
    }
    else {
      myExports.resize(numImages);
      std::copy(exports.begin()+myImageID*numImages, exports.begin()+(myImageID+1)*numImages, myExports.begin() );
    }
    // do posts, one Packet to each image
    Kokkos::View<Packet*, Kokkos::HostSpace> imports("imports", 1*distributor.getTotalReceiveLength());
    Kokkos::View<const Packet*, Kokkos::HostSpace> myExportsConst(myExports.data(), myExports.size());
    myOut << "Call doPostsAndWaits" << endl;
    distributor.doPostsAndWaits(myExportsConst, 1, imports);

    myOut << "Test results" << endl;
    // imports[i] came from image i. it was element "myImageID" in his "myExports" vector.
    // it corresponds to element i*numImages+myImageID in the global export vector
    // make a copy of the corresponding entries in the global vector, then compare these against the
    // entries that I received
    Array<Packet> expectedImports(numImages);
    {
      typename Array<Packet>::iterator eI = expectedImports.begin(),
                                        E = exports.begin()+myImageID;
      int left = numImages;
      while (true) {
        *eI = *E;
        if (--left > 0) {
          eI++;
          E += numImages;
        }
        else {
          break;
        }
      }
    }
    // check the values
    TEST_COMPARE_ARRAYS(expectedImports,imports);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, doPostsNonContig, Packet )
  {
    typedef Teuchos::ScalarTraits<Packet>   PT;
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();

    for (int run=0; run<numRuns; run++) {
      // send data to each image, including myself
      size_t numRemoteIDs = 0;
      // exportImageIDs = {0, 1, 2, ..., numImages-1, 0, 1, 2, ..., numImages-1}
      //
      // on root node only, put some invalid nodes in the middle, corresponding to untouched data in import/export buffers
      // like so:
      // exportImageIDs = {0, 1, 2, ..., numImages-1, -1, -1, 0, 1, 2, ..., numImages-1}
      Array<int> exportImageIDs;
      if (myImageID == 0) {
        exportImageIDs.reserve(2*numImages+2);
        for(int i=0; i < numImages; ++i) {
          exportImageIDs.push_back(i);
        }
        exportImageIDs.push_back(-1);
        exportImageIDs.push_back(-1);
        for(int i=0; i < numImages; ++i) {
          exportImageIDs.push_back(i);
        }
      }
      else {
        exportImageIDs.reserve(2*numImages);
        for(int i=0; i < numImages; ++i) {
          exportImageIDs.push_back(i);
        }
        for(int i=0; i < numImages; ++i) {
          exportImageIDs.push_back(i);
        }
      }
      Distributor distributor(comm);
      numRemoteIDs = distributor.createFromSends(exportImageIDs);

      // Make sure that Distributor output doesn't cause a hang.
      distributor.describe (out, Teuchos::VERB_EXTREME);

      TEST_EQUALITY(numRemoteIDs, as<size_t>(2*numImages));
      // generate global random data set: each image sends 2 packets to each image
      // we need 2*numImages*numImages "unique" values (we don't want redundant data allowing false positives)
      // root node generates all values, sends them to the others.
      Array<Packet> exports(numImages*2*numImages);
      if (myImageID == 0) {
        for (int i=0; i<2*numImages*numImages; ++i) {
          exports[i] = PT::random();
        }
      }
      // broadcast
      broadcast(*comm,0,exports());
      // pick a subset of entries to post
      Array<Packet> myExports(0);
      if (myImageID == 0) {
        myExports.resize(2*numImages+2,PT::zero());
        for (int i=0; i<numImages; ++i) {
          myExports[i] = exports[i];
        }
        for (int i=0; i<numImages; ++i) {
          myExports[numImages+2+i] = exports[numImages+i];
        }
      }
      else {
        myExports.resize(2*numImages,PT::zero());
        std::copy(exports.begin()+myImageID*2*numImages, exports.begin()+(myImageID+1)*2*numImages, myExports.begin() );
      }
      // do posts, one Packet to each image
      Kokkos::View<Packet*, Kokkos::HostSpace> imports("imports", 1*distributor.getTotalReceiveLength());
      Kokkos::View<const Packet*, Kokkos::HostSpace> myExportsConst(myExports.data(), myExports.size());
      distributor.doPostsAndWaits(myExportsConst, 1, imports);
      // imports[i] came from image i. it was element "myImageID" in his "myExports" vector.
      // it corresponds to element i*numImages+myImageID in the global export vector
      // make a copy of the corresponding entries in the global vector, then compare these against the
      // entries that I received
      Array<Packet> expectedImports(2*numImages,PT::zero());
      {
        typename Array<Packet>::iterator eI = expectedImports.begin(),
                 E = exports.begin()+myImageID;
        for (int i=0; i<numImages-1; ++i) {
          (*eI++) = *E;
          E += numImages;
          (*eI++) = *E;
          E += numImages;
        }
        (*eI++) = *E;
        E += numImages;
        (*eI++) = *E;
      }
      // check the values
      TEST_COMPARE_ARRAYS(expectedImports,imports);
      // All procs fail if any proc fails
      int globalSuccess_int = -1;
      reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
      TEST_EQUALITY_CONST( globalSuccess_int, 0 );
    }
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, badArgsFromRecvs, Ordinal )
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    const bool debug = Tpetra::Details::Behavior::debug("Distributor");

    // Each (MPI) process i sends to process i+1.  For the last
    // process, calling createFromRecvs with these data result in an
    // invalid process id.  In debug mode, this should throw an
    // exception on every process.

    if (debug) {
      {
        Distributor distributor(comm);
        Array<Ordinal> exportIDs;
        Array<int> exportNodeIDs;
        TEST_THROW( distributor.createFromRecvs<Ordinal>( tuple<Ordinal>(0), tuple<int>(myImageID+1), exportIDs, exportNodeIDs), std::runtime_error );
      }
      {
        Distributor distributor(comm);
        Array<Ordinal> exportIDs;
        Array<int> exportNodeIDs;
        TEST_THROW( distributor.createFromRecvs<Ordinal>( tuple<Ordinal>(0,0,0), tuple<int>(0,myImageID+1,0), exportIDs, exportNodeIDs), std::runtime_error );
      }
      // create from recvs with conflicting sizes, but otherwise valid entries
      {
        Distributor distributor(comm);
        Array<Ordinal> exportIDs;
        Array<int> exportNodeIDs;
        TEST_THROW( distributor.createFromRecvs<Ordinal>( tuple<Ordinal>(0), tuple<int>(0,0), exportIDs, exportNodeIDs), std::runtime_error );
      }
      {
        Distributor distributor(comm);
        Array<Ordinal> exportIDs;
        Array<int> exportNodeIDs;
        TEST_THROW( distributor.createFromRecvs<Ordinal>( tuple<Ordinal>(0,0), tuple<int>(0), exportIDs, exportNodeIDs), std::runtime_error );
      }
    }
    else {
      out << "Debug mode not enabled; set TPETRA_DEBUG=Distributor "
        "to test." << std::endl;
    }

    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromRecvs, Ordinal )
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    const int length = numImages;
    // fill remoteImageIDs with {0, 1, 2, ... length-1}
    // we'll receive one GID from every image
    //
    // fill remoteGIDs with row from generator
    // we'll receive generateValue(i,myImageID) from proc "i"
    // "i" sends us generateValue(i,myImageID)
    // similarly, we send generateValue(myImageID,i) to proc "i"
    Array<int> importImageIDs;
    Array<Ordinal> importGIDs;
    importImageIDs.reserve(length);
    importGIDs.reserve(length);
    for(int i=0; i < length; ++i) {
      importImageIDs.push_back(i);
      importGIDs.push_back( as<Ordinal>(generateValue(i, myImageID)) );
    }
    Distributor distributor(comm);
    Array<int> exportImageIDs;
    Array<Ordinal> exportGIDs;
    distributor.createFromRecvs<Ordinal>(importGIDs, importImageIDs, exportGIDs, exportImageIDs);

    // Make sure that Distributor output doesn't cause a hang.
    distributor.describe (out, Teuchos::VERB_EXTREME);

    TEST_EQUALITY(exportGIDs.size(), exportImageIDs.size());  // should *always* be the case
    Array<Ordinal> expectedGIDs;
    for (int i=0; i < length; ++i) {
      expectedGIDs.push_back( as<Ordinal>(generateValue(myImageID,i)) );
    }
    TEST_COMPARE_ARRAYS(importImageIDs, exportImageIDs);
    TEST_COMPARE_ARRAYS(expectedGIDs, exportGIDs);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
#ifdef HAVE_TPETRA_DISTRIBUTOR_TIMINGS
    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();
#endif
  }


  //
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#ifdef HAVE_TPETRA_DEBUG
#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromRecvs, ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, badArgsFromRecvs, ORDINAL )
#else
#   define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, createFromRecvs, ORDINAL )
#endif // HAVE_TPETRA_DEBUG

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD

    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, doPostsContig,    double )
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Distributor, doPostsNonContig,    double )
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_FLOAT( Distributor, doPostsContig )
    UNIT_TEST_GROUP_ORDINAL(int)

# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( Distributor, doPostsContig )
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( Distributor, doPostsNonContig )

    typedef short int ShortInt;
    UNIT_TEST_GROUP_ORDINAL(ShortInt)
    UNIT_TEST_GROUP_ORDINAL(int)
    typedef long int LongInt;
    UNIT_TEST_GROUP_ORDINAL(LongInt)
#   ifdef HAVE_TPETRA_INT_LONG_LONG
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard scopeGuard(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
