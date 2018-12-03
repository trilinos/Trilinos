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

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <type_traits>      // std::is_same

namespace {      // (anonymous)

using std::endl;
using std::max;
using std::min;
using Teuchos::arcp;
using Teuchos::arcpClone;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::null;
using Teuchos::outArg;
using Teuchos::ParameterList;
using Teuchos::parameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::tuple;
using Teuchos::VERB_EXTREME;
using Teuchos::VERB_HIGH;
using Teuchos::VERB_LOW;
using Teuchos::VERB_MEDIUM;
using Teuchos::VERB_NONE;
using Tpetra::DynamicProfile;
using Tpetra::ProfileType;
using Tpetra::StaticProfile;
using Tpetra::TestingUtilities::getDefaultComm;
typedef Tpetra::global_size_t GST;

double      errorTolSlack = 1e+1;
std::string filedir;

#define STD_TESTS(graph)                                                                                      \
    {                                                                                                         \
        auto   STCOMM   = graph.getComm();                                                                    \
        auto   STMYGIDS = graph.getRowMap()->getNodeElementList();                                            \
        size_t STMAX    = 0;                                                                                  \
                                                                                                              \
        for(size_t STR = 0; STR < graph.getNodeNumRows(); ++STR)                                              \
        {                                                                                                     \
            TEST_EQUALITY(graph.getNumEntriesInLocalRow(STR), graph.getNumEntriesInGlobalRow(STMYGIDS[STR])); \
            STMAX = std::max(STMAX, graph.getNumEntriesInLocalRow(STR));                                      \
        }                                                                                                     \
        TEST_EQUALITY(graph.getNodeMaxNumRowEntries(), STMAX);                                                \
        GST STGMAX;                                                                                           \
        Teuchos::reduceAll<int, GST>(*STCOMM, Teuchos::REDUCE_MAX, STMAX, Teuchos::outArg(STGMAX));           \
        TEST_EQUALITY(graph.getGlobalMaxNumRowEntries(), STGMAX);                                             \
    }


TEUCHOS_STATIC_SETUP()
{
    Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption("filedir", &filedir, "Directory of expected input files.");
    clp.addOutputSetupOptions(true);
    clp.setOption("test-mpi",
                  "test-serial",
                  &Tpetra::TestingUtilities::testMpi,
                  "Test MPI (if available) or force test of serial.  In a serial build,"
                  " this option is ignored and a serial comm is always used.");
    clp.setOption("error-tol-slack", &errorTolSlack, "Slack off of machine epsilon used to check test results");
}

//
// UNIT TESTS
//
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, Swap, LO, GO, Node)
{
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;

    typedef Tpetra::CrsGraph<LO, GO, Node> graph_t;
    typedef Tpetra::Map<LO, GO, Node>      map_t;

//    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();

    auto initialComm = getDefaultComm();
    TEUCHOS_TEST_FOR_EXCEPTION(initialComm->getSize() < 2, std::runtime_error, "This test requires at least two processors.");

    // Set up a communicator that has exactly two processors in it for the actual test.
    const int color = (initialComm->getRank() < 2) ? 0 : 1;
    const int key   = 0;
    auto      comm  = initialComm->split(color, key);

    if(0 == color)
    {
        const int numProcs = comm->getSize();
        const int myRank   = comm->getRank();
        assert(numProcs == 2);

        out << ">>> CrsGraph::swap() Unit Test" << std::endl;
        out << ">>> numProcs: " << numProcs << std::endl;

        success               = true;
//        const size_t numLocal = 1;

        //        RCP<const map_t> rmap = rcp(new map_t(INVALID, numLocal, 0, comm));
        //        RCP<const map_t> cmap = rcp(new map_t(INVALID, numLocal, 0, comm));

        // Set up Row Map
        const GO gblNumInds = 5;
        const LO lclNumInds = (myRank == 0) ? 3 : 2;      // 3 on p0, 2 on p1 (works b/c we have exactly 2 processes)

        std::vector<GO> myGIDs;
        if(0 == myRank)
        {
            myGIDs = {0, 1, 3};
        }
        else if(1 == myRank)
        {
            myGIDs = {7, 10};
        }

        RCP<const map_t> rowMap(new map_t(gblNumInds, myGIDs.data(), lclNumInds, 0, comm));

        Teuchos::ArrayRCP<size_t> numEntPerRow(lclNumInds);
        if(0 == myRank)
        {
            numEntPerRow[0] = 2;      // (0,0), (0,11)
            numEntPerRow[1] = 2;      // (1,3), (1,4)
            numEntPerRow[2] = 1;      // (3,2)
        }
        else if(1 == myRank)
        {
            numEntPerRow[0] = 2;      // (7,5), (7,7)
            numEntPerRow[1] = 1;      // (10,6)
        }
        RCP<graph_t> G1(new graph_t(rowMap, numEntPerRow, Tpetra::StaticProfile));

        std::vector<GO> myGblInds(2);
        if(0 == myRank)
        {
            myGblInds[0] = 0;
            myGblInds[1] = 11;
            G1->insertGlobalIndices(0, 2, myGblInds.data());      // (0,0), (0,11)
            myGblInds[0] = 3;
            myGblInds[1] = 4;
            G1->insertGlobalIndices(1, 2, myGblInds.data());      // (1,3), (1,4)
            myGblInds[0] = 2;
            G1->insertGlobalIndices(3, 1, myGblInds.data());      // (3,2)
        }
        else if(1 == myRank)
        {
            myGblInds[0] = 5;
            myGblInds[1] = 7;
            G1->insertGlobalIndices(7, 2, myGblInds.data());      // (7,5), (7,7)
            myGblInds[0] = 6;
            G1->insertGlobalIndices(10, 1, myGblInds.data());      // (10,6)
        }


        RCP<const map_t> domainMap = rowMap;

        const GO         indexBase = 0;
        RCP<const map_t> rangeMap(new map_t(8, indexBase, comm));

        G1->fillComplete(domainMap, rangeMap);

        // RCP<graph_t> graph = rcp(new graph_t(rmap, cmap, 2, StaticProfile));
    }
    else
    {
        throw std::runtime_error("FAILED on comm split!");
    }
}


////
#if 0
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(CrsGraph, WithColMap, LO, GO, Node)
{
    using Teuchos::Comm;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef Tpetra::CrsGraph<LO, GO, Node> GRAPH;
    typedef Tpetra::Map<LO, GO, Node>      map_type;

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid();

    RCP<const Comm<int>> comm     = getDefaultComm();
    const int            numProcs = comm->getSize();
    if(numProcs > 1)
    {
        const size_t        numLocal = 1;
        RCP<const map_type> rmap     = rcp(new map_type(INVALID, numLocal, 0, comm));
        RCP<const map_type> cmap     = rcp(new map_type(INVALID, numLocal, 0, comm));

        // must allocate enough for all submitted indices.
        RCP<GRAPH> G = rcp(new GRAPH(rmap, cmap, 2, StaticProfile));
        TEST_EQUALITY_CONST(G->hasColMap(), true);
        const GO myrowind = rmap->getGlobalElement(0);

        // mfh 16 May 2013: insertGlobalIndices doesn't do column Map
        // filtering anymore, so we have to test whether each of the
        // column indices to insert is invalid.
        if(cmap->isNodeGlobalElement(myrowind))
        {
            if(cmap->isNodeGlobalElement(myrowind + 1))
            {
                TEST_NOTHROW(G->insertGlobalIndices(myrowind, tuple<GO>(myrowind, myrowind + 1)));
            }
            else
            {
                TEST_NOTHROW(G->insertGlobalIndices(myrowind, tuple<GO>(myrowind)));
            }
        }
        TEST_NOTHROW(G->fillComplete());
        TEST_EQUALITY(G->getRowMap(), rmap);
        TEST_EQUALITY(G->getColMap(), cmap);
        TEST_EQUALITY(G->getNumEntriesInGlobalRow(myrowind), 1);
    }
    // All procs fail if any node fails
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int>(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));

    if(gblSuccess == 1)
    {
        out << "Succeeded on all processes!" << endl;
    }
    else
    {
        out << "FAILED on at least one process!" << endl;
    }
    TEST_EQUALITY_CONST(gblSuccess, 1);
}
#endif


  //
  // INSTANTIATIONS
  //



  // Tests to build and run in both debug and release modes.  We will
  // instantiate them over all enabled local ordinal (LO), global
  // ordinal (GO), and Kokkos Node (NODE) types.
#define UNIT_TEST_GROUP_DEBUG_AND_RELEASE(LO, GO, NODE) \
      TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(CrsGraph, Swap, LO, GO, NODE)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP_DEBUG_AND_RELEASE)

  }      // namespace
