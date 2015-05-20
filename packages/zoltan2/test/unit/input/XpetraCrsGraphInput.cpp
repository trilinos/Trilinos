// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Basic testing of Zoltan2::XpetraCrsGraphAdapter
/*!  \file XpetraCrsGraphInput.cpp
 *   \brief Test of Zoltan2::XpetraCrsGraphAdapter class.
 *  \todo add weights and coordinates
 */

#include <string>

#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::Array;
using Teuchos::ArrayView;

typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tgraph_t;
typedef Xpetra::CrsGraph<zlno_t, zgno_t, znode_t> xgraph_t;
typedef Epetra_CrsGraph egraph_t;

void printGraph(RCP<const Comm<int> > &comm, zlno_t nvtx,
    const zgno_t *vtxIds, const zlno_t *offsets, const zgno_t *edgeIds)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (zlno_t i=0; i < nvtx; i++){
        std::cout << " vertex " << vtxIds[i] << ": ";
        for (zlno_t j=offsets[i]; j < offsets[i+1]; j++){
          std::cout << edgeIds[j] << " ";
        }
        std::cout << std::endl;
      }
      std::cout.flush();
    }
    comm->barrier();
  }
  comm->barrier();
}

template <typename User>
int verifyInputAdapter(
  Zoltan2::XpetraCrsGraphAdapter<User> &ia,
  tgraph_t &graph
)
{
  RCP<const Comm<int> > comm = graph.getComm();
  int fail = 0, gfail=0;

  if (!fail && 
      ia.getLocalNumVertices() != graph.getNodeNumRows())
    fail = 4;

  if (!fail && 
      ia.getLocalNumEdges() != graph.getNodeNumEntries())
      fail = 6;

  gfail = globalFail(comm, fail);

  const zgno_t *vtxIds=NULL, *edgeIds=NULL;
  const zlno_t *offsets=NULL;
  size_t nvtx=0;

  if (!gfail){

    nvtx = ia.getLocalNumVertices();
    ia.getVertexIDsView(vtxIds);
    ia.getEdgesView(offsets, edgeIds);

    if (nvtx != graph.getNodeNumRows())
      fail = 8;

    gfail = globalFail(comm, fail);

    if (gfail == 0){
      printGraph(comm, nvtx, vtxIds, offsets, edgeIds);
    }
    else{
      if (!fail) fail = 10;
    }
  }
  return fail;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int fail = 0, gfail=0;

  // Create an object that can give us test Tpetra, Xpetra
  // and Epetra graphs for testing.

  RCP<UserInputForTests> uinput;

  try{
    uinput =
      rcp(new UserInputForTests(testDataFilePath,std::string("simple"), comm, true));
  }
  catch(std::exception &e){
    TEST_FAIL_AND_EXIT(*comm, 0, string("input ")+e.what(), 1);
  }

  RCP<tgraph_t> tG;     // original graph (for checking)
  RCP<tgraph_t> newG;   // migrated graph

  tG = uinput->getUITpetraCrsGraph();
  size_t nvtx = tG->getNodeNumRows();
  ArrayView<const zgno_t> rowGids = tG->getRowMap()->getNodeElementList();

  // To test migration in the input adapter we need a Solution
  // object.  The Solution needs an IdentifierMap.
  // Our solution just assigns all objects to part zero.

  typedef Zoltan2::IdentifierMap<tgraph_t> idmap_t;

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  ArrayRCP<const zgno_t> gidArray = arcpFromArrayView(rowGids);
  RCP<const idmap_t> idMap = rcp(new idmap_t(env, comm, gidArray));

  int nWeights = 1;

  typedef Zoltan2::XpetraCrsGraphAdapter<tgraph_t>  adapter_t;
  typedef Zoltan2::PartitioningSolution<adapter_t> soln_t;
  typedef adapter_t::part_t part_t;

  part_t *p = new part_t [nvtx];
  memset(p, 0, sizeof(part_t) * nvtx);
  ArrayRCP<part_t> solnParts(p, 0, nvtx, true);

  soln_t solution(env, comm, idMap, nWeights);
  solution.setParts(solnParts);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::CrsGraph
  if (!gfail){
    RCP<const tgraph_t> ctG = rcp_const_cast<const tgraph_t>(tG);
    RCP<Zoltan2::XpetraCrsGraphAdapter<tgraph_t> > tGInput;

    try {
      tGInput =
        rcp(new Zoltan2::XpetraCrsGraphAdapter<tgraph_t>(ctG));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0,
        string("XpetraCrsGraphAdapter ")+e.what(), 1);
    }

    if (rank==0)
      std::cout << "Input adapter for Tpetra::CrsGraph" << std::endl;

    fail = verifyInputAdapter<tgraph_t>(*tGInput, *tG);

    gfail = globalFail(comm, fail);

    if (!gfail){
      tgraph_t *mMigrate = NULL;
      try{
        tGInput->applyPartitioningSolution( *tG, mMigrate, solution);
        newG = rcp(mMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);

      if (!gfail){
        RCP<const tgraph_t> cnewG = rcp_const_cast<const tgraph_t>(newG);
        RCP<Zoltan2::XpetraCrsGraphAdapter<tgraph_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraCrsGraphAdapter<tgraph_t>(cnewG));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0,
            string("XpetraCrsGraphAdapter 2 ")+e.what(), 1);
        }

        if (rank==0){
          std::cout <<
           "Input adapter for Tpetra::CrsGraph migrated to proc 0" <<
           std::endl;
        }
        fail = verifyInputAdapter<tgraph_t>(*newInput, *newG);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::CrsGraph
  if (!gfail){
    RCP<xgraph_t> xG = uinput->getUIXpetraCrsGraph();
    RCP<const xgraph_t> cxG = rcp_const_cast<const xgraph_t>(xG);
    RCP<Zoltan2::XpetraCrsGraphAdapter<xgraph_t> > xGInput;

    try {
      xGInput =
        rcp(new Zoltan2::XpetraCrsGraphAdapter<xgraph_t>(cxG));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0,
        string("XpetraCrsGraphAdapter 3 ")+e.what(), 1);
    }

    if (rank==0){
      std::cout << "Input adapter for Xpetra::CrsGraph" << std::endl;
    }
    fail = verifyInputAdapter<xgraph_t>(*xGInput, *tG);

    gfail = globalFail(comm, fail);

    if (!gfail){
      xgraph_t *mMigrate =NULL;
      try{
        xGInput->applyPartitioningSolution( *xG, mMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);

      if (!gfail){
        RCP<const xgraph_t> cnewG(mMigrate);
        RCP<Zoltan2::XpetraCrsGraphAdapter<xgraph_t> > newInput;
        try{
          newInput =
            rcp(new Zoltan2::XpetraCrsGraphAdapter<xgraph_t>(cnewG));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0,
            string("XpetraCrsGraphAdapter 4 ")+e.what(), 1);
        }

        if (rank==0){
          std::cout <<
           "Input adapter for Xpetra::CrsGraph migrated to proc 0" <<
           std::endl;
        }
        fail = verifyInputAdapter<xgraph_t>(*newInput, *newG);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

#ifdef HAVE_EPETRA_DATA_TYPES
  /////////////////////////////////////////////////////////////
  // User object is Epetra_CrsGraph
  if (!gfail){
    RCP<egraph_t> eG = uinput->getUIEpetraCrsGraph();
    RCP<const egraph_t> ceG = rcp_const_cast<const egraph_t>(eG);
    RCP<Zoltan2::XpetraCrsGraphAdapter<egraph_t> > eGInput;

    try {
      eGInput =
        rcp(new Zoltan2::XpetraCrsGraphAdapter<egraph_t>(ceG));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0,
        string("XpetraCrsGraphAdapter 5 ")+e.what(), 1);
    }

    if (rank==0){
      std::cout << "Input adapter for Epetra_CrsGraph" << std::endl;
    }
    fail = verifyInputAdapter<egraph_t>(*eGInput, *tG);

    gfail = globalFail(comm, fail);

    if (!gfail){
      egraph_t *mMigrate =NULL;
      try{
        eGInput->applyPartitioningSolution( *eG, mMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);

      if (!gfail){
        RCP<const egraph_t> cnewG(mMigrate, true);
        RCP<Zoltan2::XpetraCrsGraphAdapter<egraph_t> > newInput;
        try{
          newInput =
            rcp(new Zoltan2::XpetraCrsGraphAdapter<egraph_t>(cnewG));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0,
            string("XpetraCrsGraphAdapter 6 ")+e.what(), 1);
        }

        if (rank==0){
          std::cout <<
           "Input adapter for Epetra_CrsGraph migrated to proc 0" <<
           std::endl;
        }
        fail = verifyInputAdapter<egraph_t>(*newInput, *newG);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }
#endif

  /////////////////////////////////////////////////////////////
  // DONE

  if (rank==0)
    std::cout << "PASS" << std::endl;
}
