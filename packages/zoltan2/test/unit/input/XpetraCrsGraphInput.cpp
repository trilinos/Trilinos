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

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayView;

typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> tgraph_t;
typedef Xpetra::CrsGraph<zlno_t, zgno_t, znode_t> xgraph_t;

template<typename offset_t>
void printGraph(RCP<const Comm<int> > &comm, zlno_t nvtx,
    const zgno_t *vtxIds, const offset_t *offsets, const zgno_t *edgeIds)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (zlno_t i=0; i < nvtx; i++){
        std::cout << " vertex " << vtxIds[i] << ": ";
        for (offset_t j=offsets[i]; j < offsets[i+1]; j++){
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
  typedef typename Zoltan2::InputTraits<User>::offset_t offset_t;

  RCP<const Comm<int> > comm = graph.getComm();
  int fail = 0, gfail=0;

  if (!fail && 
      ia.getLocalNumVertices() != graph.getNodeNumRows())
    fail = 4;

  if (!fail && 
      ia.getLocalNumEdges() != graph.getNodeNumEntries())
      fail = 6;

  gfail = globalFail(*comm, fail);

  const zgno_t *vtxIds=NULL, *edgeIds=NULL;
  const offset_t *offsets=NULL;
  size_t nvtx=0;

  if (!gfail){

    nvtx = ia.getLocalNumVertices();
    ia.getVertexIDsView(vtxIds);
    ia.getEdgesView(offsets, edgeIds);

    if (nvtx != graph.getNodeNumRows())
      fail = 8;

    gfail = globalFail(*comm, fail);

    if (gfail == 0){
      printGraph<offset_t>(comm, nvtx, vtxIds, offsets, edgeIds);
    }
    else{
      if (!fail) fail = 10;
    }
  }
  return fail;
}

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int fail = 0, gfail=0;
  bool aok = true;

  // Create an object that can give us test Tpetra, Xpetra
  // and Epetra graphs for testing.

  RCP<UserInputForTests> uinput;
  Teuchos::ParameterList params;
  params.set("input file", "simple");
  params.set("file type", "Chaco");

  try{
    uinput = rcp(new UserInputForTests(params, comm));
  }
  catch(std::exception &e){
    aok = false;
    std::cout << e.what() << std::endl;
  }
  TEST_FAIL_AND_EXIT(*comm, aok, "input ", 1);

  RCP<tgraph_t> tG;     // original graph (for checking)
  RCP<tgraph_t> newG;   // migrated graph

  tG = uinput->getUITpetraCrsGraph();
  size_t nvtx = tG->getNodeNumRows();

  // To test migration in the input adapter we need a Solution object.  
  // Our solution just assigns all objects to part zero.

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));

  int nWeights = 1;

  typedef Zoltan2::XpetraCrsGraphAdapter<tgraph_t>  adapter_t;
  typedef Zoltan2::PartitioningSolution<adapter_t> soln_t;
  typedef adapter_t::part_t part_t;

  part_t *p = new part_t [nvtx];
  memset(p, 0, sizeof(part_t) * nvtx);
  ArrayRCP<part_t> solnParts(p, 0, nvtx, true);

  soln_t solution(env, comm, nWeights);
  solution.setParts(solnParts);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::CrsGraph
  if (!gfail){
    if (rank==0)
      std::cout << "Input adapter for Tpetra::CrsGraph" << std::endl;

    RCP<const tgraph_t> ctG = rcp_const_cast<const tgraph_t>(tG);
    RCP<Zoltan2::XpetraCrsGraphAdapter<tgraph_t> > tGInput;

    try {
      tGInput =
        rcp(new Zoltan2::XpetraCrsGraphAdapter<tgraph_t>(ctG));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsGraphAdapter ", 1);

    fail = verifyInputAdapter<tgraph_t>(*tGInput, *tG);

    gfail = globalFail(*comm, fail);

    if (!gfail){
      tgraph_t *mMigrate = NULL;
      try{
        tGInput->applyPartitioningSolution( *tG, mMigrate, solution);
        newG = rcp(mMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(*comm, fail);

      if (!gfail){
        RCP<const tgraph_t> cnewG = rcp_const_cast<const tgraph_t>(newG);
        RCP<Zoltan2::XpetraCrsGraphAdapter<tgraph_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraCrsGraphAdapter<tgraph_t>(cnewG));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsGraphAdapter 2 ", 1);

        if (rank==0){
          std::cout <<
           "Input adapter for Tpetra::CrsGraph migrated to proc 0" <<
           std::endl;
        }
        fail = verifyInputAdapter<tgraph_t>(*newInput, *newG);
        if (fail) fail += 100;
        gfail = globalFail(*comm, fail);
      }
    }
    if (gfail){
      printFailureCode(*comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::CrsGraph
  if (!gfail){
    if (rank==0)
      std::cout << "Input adapter for Xpetra::CrsGraph" << std::endl;
     
    RCP<xgraph_t> xG = uinput->getUIXpetraCrsGraph();
    RCP<const xgraph_t> cxG = rcp_const_cast<const xgraph_t>(xG);
    RCP<Zoltan2::XpetraCrsGraphAdapter<xgraph_t> > xGInput;

    try {
      xGInput =
        rcp(new Zoltan2::XpetraCrsGraphAdapter<xgraph_t>(cxG));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsGraphAdapter 3 ", 1);

    fail = verifyInputAdapter<xgraph_t>(*xGInput, *tG);

    gfail = globalFail(*comm, fail);

    if (!gfail){
      xgraph_t *mMigrate =NULL;
      try{
        xGInput->applyPartitioningSolution( *xG, mMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(*comm, fail);

      if (!gfail){
        RCP<const xgraph_t> cnewG(mMigrate);
        RCP<Zoltan2::XpetraCrsGraphAdapter<xgraph_t> > newInput;
        try{
          newInput =
            rcp(new Zoltan2::XpetraCrsGraphAdapter<xgraph_t>(cnewG));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsGraphAdapter 4 ", 1);

        if (rank==0){
          std::cout <<
           "Input adapter for Xpetra::CrsGraph migrated to proc 0" <<
           std::endl;
        }
        fail = verifyInputAdapter<xgraph_t>(*newInput, *newG);
        if (fail) fail += 100;
        gfail = globalFail(*comm, fail);
      }
    }
    if (gfail){
      printFailureCode(*comm, fail);
    }
  }

#ifdef HAVE_EPETRA_DATA_TYPES
  /////////////////////////////////////////////////////////////
  // User object is Epetra_CrsGraph
  typedef Epetra_CrsGraph egraph_t;
  if (!gfail){
    if (rank==0)
      std::cout << "Input adapter for Epetra_CrsGraph" << std::endl;

    RCP<egraph_t> eG = uinput->getUIEpetraCrsGraph();
    RCP<const egraph_t> ceG = rcp_const_cast<const egraph_t>(eG);
    RCP<Zoltan2::XpetraCrsGraphAdapter<egraph_t> > eGInput;

    bool goodAdapter = true;
    try {
      eGInput =
        rcp(new Zoltan2::XpetraCrsGraphAdapter<egraph_t>(ceG));
    }
    catch (std::exception &e){
      if (std::is_same<znode_t, Xpetra::EpetraNode>::value) {
        aok = false;
        goodAdapter = false;
        std::cout << e.what() << std::endl;
      }
      else {
        // We expect an error from Xpetra when znode_t != Xpetra::EpetraNode
        // Ignore it, but skip tests using this matrix adapter.
        std::cout << "Node type is not supported by Xpetra's Epetra interface;"
                  << " Skipping this test." << std::endl;
        std::cout << "FYI:  Here's the exception message: " << std::endl
                  << e.what() << std::endl;
        goodAdapter = false;
      }
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsGraphAdapter 5 ", 1);

    if (goodAdapter) {
      fail = verifyInputAdapter<egraph_t>(*eGInput, *tG);

      gfail = globalFail(*comm, fail);

      if (!gfail){
        egraph_t *mMigrate =NULL;
        try{
          eGInput->applyPartitioningSolution( *eG, mMigrate, solution);
        }
        catch (std::exception &e){
          fail = 11;
        }

        gfail = globalFail(*comm, fail);

        if (!gfail){
          RCP<const egraph_t> cnewG(mMigrate, true);
          RCP<Zoltan2::XpetraCrsGraphAdapter<egraph_t> > newInput;
          try{
            newInput =
              rcp(new Zoltan2::XpetraCrsGraphAdapter<egraph_t>(cnewG));
          }
          catch (std::exception &e){
            aok = false;
            std::cout << e.what() << std::endl;
          }
          TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsGraphAdapter 6 ", 1);

          if (rank==0){
            std::cout <<
             "Input adapter for Epetra_CrsGraph migrated to proc 0" <<
             std::endl;
          }
          fail = verifyInputAdapter<egraph_t>(*newInput, *newG);
          if (fail) fail += 100;
          gfail = globalFail(*comm, fail);
        }
      }
      if (gfail){
        printFailureCode(*comm, fail);
      }
    }
  }
#endif

  /////////////////////////////////////////////////////////////
  // DONE

  if (rank==0)
    std::cout << "PASS" << std::endl;
}
