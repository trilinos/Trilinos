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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ParameterList.hpp"

namespace { // (anonymous)
  TEUCHOS_UNIT_TEST(CrsGraph, DetectNonlocals)
  {
    using Tpetra::TestingUtilities::getDefaultComm;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using std::endl;
    using crs_graph_type = Tpetra::CrsGraph<>;
    using map_type = Tpetra::Map<>;
    using LO = crs_graph_type::local_ordinal_type;
    using GO = crs_graph_type::global_ordinal_type;

    const bool debug = Tpetra::Details::Behavior::debug("CrsGraph");
    if (! debug) {
      out << "This test only does anything in debug mode.  "
        "Try rerunning with the environment variable "
        "TPETRA_DEBUG=CrsGraph set." << endl;
      return;
    }
    out << "Test that Tpetra::CrsGraph::fillComplete reports an "
      "error if the \"No Nonlocal Changes\" parameter is true, but "
      "there were nonlocal inserts." << endl;

    auto comm = getDefaultComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    TEST_EQUALITY( numProcs, 2 );

    const LO lclNumRows = 1;
    const GO gblNumRows = 2;
    const GO indexBase = 0;
    RCP<const map_type> rowMap =
      rcp(new map_type(gblNumRows, lclNumRows, indexBase, comm));

    for (const bool assertNoNonlocalChanges : {false, true}) {
      crs_graph_type G(rowMap, 1);

      // Attempt an insert into a nonowned row, but only on Process 1.
      // In debug mode, the processes should do an all-reduce to check
      // whether any process had a nonowned insert.
      if (myRank == 1) {
        const GO gblRow(0);
        const GO gblColInd(0);
        TEST_ASSERT( ! rowMap->isNodeGlobalElement(gblRow) );
        G.insertGlobalIndices(gblRow, LO(1), &gblColInd);
      }

      auto params = Teuchos::parameterList("CrsGraph::fillComplete");
      // If set to true, this is a lie; the point is that in debug
      // mode, fillComplete should check whether it is a lie.  Don't
      // use TEST_THROW because that requires a specific type of
      // exception; we don't care about the exception type here, as
      // long as it's a subclass of std::exception.
      if (assertNoNonlocalChanges) {
        params->set("No Nonlocal Changes", true);
      }
      bool threw = false;
      try {
        G.fillComplete(params);
      }
      catch (std::exception& e) {
        threw = true;
        out << "G.fillComplete() threw: " << e.what() << endl;
      }
      if (assertNoNonlocalChanges) {
        TEST_ASSERT( threw );
      }
      else {
        TEST_ASSERT( ! threw );
      }

      int lclSuccess = success ? 1 : 0;
      int gblSuccess = 0;
      using Teuchos::outArg;
      using Teuchos::REDUCE_MIN;
      using Teuchos::reduceAll;
      reduceAll(*comm, REDUCE_MIN, lclSuccess, outArg(gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
    }
  }
} // namespace (anonymous)
