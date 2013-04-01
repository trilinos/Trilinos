/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Teuchos_as.hpp"
#include "Tpetra_Distributor.hpp"
#include <Teuchos_Array.hpp>

// FINISH: test for createFromRecvs(), that negatives in remoteNodeIDs are met by negatives in exportNodeIDs, and that the placement is
//         is preserved. need to understand the semantics of negatives in the node list.

namespace {

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Tpetra::Distributor;
  using Tpetra::DefaultPlatform;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::tuple;

  bool testMpi = true;
  double errorTolSlack = 1e+1;

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
  }

  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
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
    // create from contiguous sends
    Distributor distributor(comm);
    numImports = distributor.createFromSends(exportImageIDs);
    // tests
    TEST_EQUALITY(numImports, as<size_t>(2*numImages));
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<size_t>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), (numImages > 1 ? 2 : 0))
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<size_t>(2*numImages));
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
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
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


// mfh 01 Apr 2013: Distributor only checks input arguments in a
// debug build, so this test is only enabled in a debug build.
#ifdef HAVE_TPETRA_DEBUG
  ////
  TEUCHOS_UNIT_TEST( Distributor, badArgsFromSends)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // each node i sends to node i+1
    // for the last node, this results in an invalid node id, which should throw an exception on
    // every node
    size_t numImports = 0; (void)numImports;
    // create from sends with bad node IDs
    {
      Distributor distributor(comm);
      TEST_THROW( numImports = distributor.createFromSends( tuple<int>(myImageID+1)), std::runtime_error );
    }
    {
      Distributor distributor(comm);
      TEST_THROW( numImports = distributor.createFromSends( tuple<int>(0,myImageID+1,0)), std::runtime_error );
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }
#endif // HAVE_TPETRA_DEBUG

  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsMixedContig)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
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
    // create from sends, contiguous and non-contiguous
    Distributor distributor(comm);
#ifdef HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
    TEST_THROW( distributor.createFromSends(exportImageIDs), std::runtime_error );
#else
    TEST_NOTHROW( numImports = distributor.createFromSends(exportImageIDs) );
    // tests
    TEST_EQUALITY(numImports, as<size_t>(2*numImages));
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<size_t>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 2);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<size_t>(2*numImages));
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
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
#endif
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsRedBlack)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 3) return;
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
    // create from contiguous sends
    Distributor distributor(comm);
    numImports = distributor.createFromSends(exportImageIDs);
    // tests
    TEST_EQUALITY(numImports, numInMyPartition);
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), numInMyPartition-1);
    TEST_EQUALITY(distributor.getNumReceives(), numInMyPartition-1);
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), (numInMyPartition > 1 ? 1 : 0));
    TEST_EQUALITY(distributor.getTotalReceiveLength(), numInMyPartition);
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
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
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsContigNoself)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
    if (numImages < 2) return;
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
    // create from contiguous sends
    Distributor distributor(comm);
    numImports = distributor.createFromSends(exportImageIDs);
    // tests
    TEST_EQUALITY(numImports, as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), false);
    TEST_EQUALITY(distributor.getNumSends(), as<size_t>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 1);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<size_t>(numImages-1));
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
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
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsContigUnordered)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    if (numImages < 3) return;
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
    // create from contiguous sends
    Distributor distributor(comm);
    numImports = distributor.createFromSends(exportImageIDs);
    // tests
    TEST_EQUALITY(numImports, as<size_t>(3*numImages));
    TEST_EQUALITY_CONST(distributor.hasSelfMessage(), true);
    TEST_EQUALITY(distributor.getNumSends(), as<size_t>(numImages-1));
    TEST_EQUALITY(distributor.getNumReceives(), as<size_t>(numImages-1));
    TEST_EQUALITY_CONST(distributor.getMaxSendLength(), 3);
    TEST_EQUALITY(distributor.getTotalReceiveLength(), as<size_t>(3*numImages));
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
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
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST( Distributor, createFromSendsNonContig)
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();

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
    // create from non-contiguous sends
    Distributor distributor(comm);
#ifdef HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
    TEST_THROW( distributor.createFromSends(exportImageIDs), std::runtime_error );
#else
    numImports = distributor.createFromSends(exportImageIDs);
    // tests
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
    {
      ArrayView<const int> imgFrom(distributor.getImagesFrom());
      ArrayView<const int> imgTo(distributor.getImagesTo());
      TEST_COMPARE_ARRAYS(imgFrom, imgTo);
    }
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
#endif
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, doPostsContig, Packet )
  {
    typedef Teuchos::ScalarTraits<Packet>   PT;
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
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
    Distributor distributor(comm);
    numRemoteIDs = distributor.createFromSends(exportImageIDs);
    TEST_EQUALITY(numRemoteIDs, as<size_t>(numImages));
    // generate global random data set: each image sends 1 packet to each image
    // we need numImages*numImages "unique" values (we don't want redundant data allowing false positives)
    // root node generates all values, sends them to the others.
    Array<Packet> exports(numImages*numImages);
    if (myImageID == 0) {
      for (int i=0; i<numImages*numImages; i++) {
        exports[i] = PT::random();
      }
    }
    // broadcast
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
    Array<Packet> imports(1*distributor.getTotalReceiveLength());
    distributor.doPostsAndWaits(myExports().getConst(), 1, imports());
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
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, doPostsNonContig, Packet )
  {
    typedef Teuchos::ScalarTraits<Packet>   PT;
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
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
#ifdef HAVE_TPETRA_THROW_EFFICIENCY_WARNINGS
    TEST_THROW( numRemoteIDs = distributor.createFromSends(exportImageIDs), std::runtime_error );
#else
    numRemoteIDs = distributor.createFromSends(exportImageIDs);
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
    Array<Packet> imports(1*distributor.getTotalReceiveLength());
    distributor.doPostsAndWaits(myExports().getConst(), 1, imports());
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
#endif
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }


// mfh 01 Apr 2013: Distributor only checks input arguments in a
// debug build, so this test is only enabled in a debug build.
#ifdef HAVE_TPETRA_DEBUG
  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, badArgsFromRecvs, Ordinal )
  {
    RCP<const Comm<int> > comm = getDefaultComm();
    const int myImageID = comm->getRank();
    // each node i sends to node i+1
    // for the last node, this results in an invalid node id, which should throw an exception on
    // every node
    // create from recvs with bad node IDs
    {
      Distributor distributor(comm);
      ArrayRCP<Ordinal> exportIDs;
      ArrayRCP<int> exportNodeIDs;
      TEST_THROW( distributor.createFromRecvs<Ordinal>( tuple<Ordinal>(0), tuple<int>(myImageID+1), exportIDs, exportNodeIDs), std::runtime_error );
    }
    {
      Distributor distributor(comm);
      ArrayRCP<Ordinal> exportIDs;
      ArrayRCP<int> exportNodeIDs;
      TEST_THROW( distributor.createFromRecvs<Ordinal>( tuple<Ordinal>(0,0,0), tuple<int>(0,myImageID+1,0), exportIDs, exportNodeIDs), std::runtime_error );
    }
    // create from recvs with conflicting sizes, but otherwise valid entries
    {
      Distributor distributor(comm);
      ArrayRCP<Ordinal> exportIDs;
      ArrayRCP<int> exportNodeIDs;
      TEST_THROW( distributor.createFromRecvs<Ordinal>( tuple<Ordinal>(0), tuple<int>(0,0), exportIDs, exportNodeIDs), std::runtime_error );
    }
    {
      Distributor distributor(comm);
      ArrayRCP<Ordinal> exportIDs;
      ArrayRCP<int> exportNodeIDs;
      TEST_THROW( distributor.createFromRecvs<Ordinal>( tuple<Ordinal>(0,0), tuple<int>(0), exportIDs, exportNodeIDs), std::runtime_error );
    }
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }
#endif // HAVE_TPETRA_DEBUG

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Distributor, createFromRecvs, Ordinal )
  {
    typedef Teuchos::OrdinalTraits<Ordinal> OT;
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
    ArrayRCP<int> exportImageIDs;
    ArrayRCP<Ordinal> exportGIDs;
    distributor.createFromRecvs<Ordinal>(importGIDs, importImageIDs, exportGIDs, exportImageIDs);
    TEST_EQUALITY(exportGIDs.size(), exportImageIDs.size());  // should *always* be the case
    Array<Ordinal> expectedGIDs;
    for (int i=0; i < length; ++i) {
      expectedGIDs.push_back( as<Ordinal>(generateValue(myImageID,i)) );
    }
    TEST_COMPARE_ARRAYS(importImageIDs, exportImageIDs);
    TEST_COMPARE_ARRAYS(expectedGIDs, exportGIDs);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
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
#   ifdef HAVE_TEUCHOS_LONG_LONG_INT
      typedef long long int LongLongInt;
      UNIT_TEST_GROUP_ORDINAL(LongLongInt)
#   endif

# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}
