// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Basic testing of Zoltan2::TpetraRowGraphAdapter
/*!  \file TpetraRowGraphAdapter.cpp
 *   \brief Test of Zoltan2::TpetraRowGraphAdapter class.
 *  \todo add weights and coordinates
 */

#include <string>

#include <Zoltan2_TpetraRowGraphAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::Comm;

typedef Tpetra::CrsGraph<zlno_t, zgno_t, znode_t> ztcrsgraph_t;
typedef Tpetra::RowGraph<zlno_t, zgno_t, znode_t> ztrowgraph_t;

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
  Zoltan2::TpetraRowGraphAdapter<User> &ia, ztrowgraph_t &graph)
{
  typedef typename Zoltan2::InputTraits<User>::offset_t offset_t;

  auto comm = graph.getComm();
  int fail = 0, gfail=0;

  if (!fail &&
      ia.getLocalNumVertices() != graph.getLocalNumRows())
    fail = 4;

  if (!fail &&
      ia.getLocalNumEdges() != graph.getLocalNumEntries())
      fail = 6;

  gfail = globalFail(*comm, fail);

  const zgno_t *vtxIds=NULL, *edgeIds=NULL;
  const offset_t *offsets=NULL;
  size_t nvtx=0;

  if (!gfail){

    nvtx = ia.getLocalNumVertices();
    ia.getVertexIDsView(vtxIds);
    ia.getEdgesView(offsets, edgeIds);

    if (nvtx != graph.getLocalNumRows())
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

  // Create an object that can give us test Tpetra graphs for testing

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

  // Input crs graph and row graph cast from it.
  auto tG = uinput->getUITpetraCrsGraph();
  auto trG = rcp_dynamic_cast<ztrowgraph_t>(tG);

  RCP<ztrowgraph_t> newG;   // migrated graph

  size_t nvtx = tG->getLocalNumRows();

  // To test migration in the input adapter we need a Solution object.

  const auto env = rcp(new Zoltan2::Environment(comm));

  int nWeights = 1;

  typedef Zoltan2::TpetraRowGraphAdapter<ztrowgraph_t>  adapter_t;
  typedef Zoltan2::PartitioningSolution<adapter_t> soln_t;
  typedef adapter_t::part_t part_t;

  part_t *p = new part_t [nvtx];
  memset(p, 0, sizeof(part_t) * nvtx);
  ArrayRCP<part_t> solnParts(p, 0, nvtx, true);

  soln_t solution(env, comm, nWeights);
  solution.setParts(solnParts);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::RowGraph
  if (!gfail){
    if (rank==0)
      std::cout << "Input adapter for Tpetra::RowGraph" << std::endl;

    RCP<const ztrowgraph_t> ctrG = rcp_const_cast<const ztrowgraph_t>(
                                   rcp_dynamic_cast<ztrowgraph_t>(tG));

    RCP<adapter_t> trGInput;

    try {
      trGInput = rcp(new adapter_t(ctrG));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "TpetraRowGraphAdapter ", 1);

    fail = verifyInputAdapter<ztrowgraph_t>(*trGInput, *trG);

    gfail = globalFail(*comm, fail);

    if (!gfail){
      ztrowgraph_t *mMigrate = NULL;
      try{
        trGInput->applyPartitioningSolution( *trG, mMigrate, solution);
        newG = rcp(mMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(*comm, fail);

      if (!gfail){
        auto cnewG = rcp_const_cast<const ztrowgraph_t>(newG);
        RCP<adapter_t> newInput;
        try{
          newInput = rcp(new adapter_t(cnewG));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "TpetraRowGraphAdapter 2 ", 1);

        if (rank==0){
          std::cout <<
           "Input adapter for Tpetra::RowGraph migrated to proc 0" <<
           std::endl;
        }
        fail = verifyInputAdapter<ztrowgraph_t>(*newInput, *newG);
        if (fail) fail += 100;
        gfail = globalFail(*comm, fail);
      }
    }
    if (gfail){
      printFailureCode(*comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // DONE

  if (rank==0)
    std::cout << "PASS" << std::endl;
}
