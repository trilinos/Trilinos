// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Basic testing of Zoltan2::XpetraCrsGraphInput

#include <string>

#include <UserInputForTests.hpp>

#include <Zoltan2_XpetraCrsGraphInput.hpp>
#include <Zoltan2_InputTraits.hpp>

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

typedef double scalar_t;
typedef int lno_t;
typedef int gno_t;
typedef Zoltan2::default_node_t node_t;

typedef UserInputForTests<scalar_t, lno_t, gno_t> uinput_t;
typedef Tpetra::CrsGraph<lno_t, gno_t, node_t> tgraph_t;
typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> xgraph_t;
typedef Epetra_CrsGraph egraph_t;

template <typename L, typename G>
  void printGraph(RCP<const Comm<int> > &comm, L nvtx,
    const G *vtxIds, const L *offsets, const G *edgeIds)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (L i=0; i < nvtx; i++){
        std::cout << " vertex " << vtxIds[i] << ": ";
        for (L j=offsets[i]; j < offsets[i+1]; j++){
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
  Zoltan2::XpetraCrsGraphInput<User> &ia, tgraph_t &graph)
{
  typedef typename Zoltan2::InputTraits<User>::scalar_t S;
  typedef typename Zoltan2::InputTraits<User>::lno_t L;
  typedef typename Zoltan2::InputTraits<User>::gno_t G;

  RCP<const Comm<int> > comm = graph.getComm();
  int fail = 0, gfail=0;

  if (!ia.haveLocalIds())
    fail = 1;

  size_t base;
  if (!fail && !ia.haveConsecutiveLocalIds(base))
    fail = 2;

  if (!fail && base != 0)
    fail = 3;

  if (!fail && ia.getLocalNumVertices() != graph.getNodeNumRows())
    fail = 4;

  if (!fail && ia.getGlobalNumVertices() != graph.getGlobalNumRows())
    fail = 5;

  if (!fail && ia.getLocalNumEdges() != graph.getNodeNumEntries())
      fail = 6;

  if (!fail && ia.getGlobalNumEdges() != graph.getGlobalNumEntries())
    fail = 7;

  gfail = globalFail(comm, fail);

  const G *vtxIds=NULL, *edgeIds=NULL;
  const L *lids=NULL, *offsets=NULL;
  size_t nvtx=0;

  if (!gfail){

    nvtx = ia.getVertexListView(vtxIds, lids, offsets, edgeIds);

    if (nvtx != graph.getNodeNumRows())
      fail = 8;
    if (!fail && lids != NULL)
      fail = 9;

    gfail = globalFail(comm, fail);

    if (gfail == 0){
      printGraph<L, G>(comm, nvtx, vtxIds, offsets, edgeIds);
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
  int nprocs = comm->getSize();
  int fail = 0, gfail=0;

  // Create object that can give us test Tpetra, Xpetra
  // and Epetra graphs for testing.

  RCP<uinput_t> uinput;

  try{
    uinput =
      rcp(new uinput_t(std::string("../data/simple.mtx"), comm));
  }
  catch(std::exception &e){
    TEST_FAIL_AND_EXIT(*comm, 0, string("input ")+e.what(), 1);
  }

  RCP<tgraph_t> tG;     // original graph (for checking)
  RCP<tgraph_t> newG;   // migrated graph

  tG = uinput->getTpetraCrsGraph();
  size_t nvtx = tG->getNodeNumRows();
  Teuchos::ArrayView<const gno_t> rowGids =
    tG->getRowMap()->getNodeElementList();

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::CrsGraph
  if (!gfail){
    RCP<const tgraph_t> ctG = rcp_const_cast<const tgraph_t>(tG);
    RCP<Zoltan2::XpetraCrsGraphInput<tgraph_t> > tGInput;

    try {
      tGInput =
        rcp(new Zoltan2::XpetraCrsGraphInput<tgraph_t>(ctG));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0,
        string("XpetraCrsGraphInput ")+e.what(), 1);
    }

    if (rank==0)
      std::cout << "Input adapter for Tpetra::CrsGraph" << std::endl;

    fail = verifyInputAdapter<tgraph_t>(*tGInput, *tG);

    gfail = globalFail(comm, fail);

    if (!gfail){
      Zoltan2::PartitioningSolution<gno_t, lno_t, lno_t>
               solution(nprocs, nvtx, 0);
      ArrayRCP<gno_t> &solnGids = solution.getGidsRCP();
      ArrayRCP<size_t> &solnParts = solution.getPartsRCP();
      for (size_t i = 0; i < nvtx; i++) solnGids[i] = rowGids[i];
      memset(solnParts.getRawPtr(), 0, sizeof(size_t) * nvtx);

      tgraph_t *mMigrate = NULL;
      try{
        tGInput->applyPartitioningSolution(*tG, mMigrate, solution);
        newG = rcp(mMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);

      if (!gfail){
        RCP<const tgraph_t> cnewG = rcp_const_cast<const tgraph_t>(newG);
        RCP<Zoltan2::XpetraCrsGraphInput<tgraph_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraCrsGraphInput<tgraph_t>(cnewG));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0,
            string("XpetraCrsGraphInput 2 ")+e.what(), 1);
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
    RCP<xgraph_t> xG = uinput->getXpetraCrsGraph();
    RCP<const xgraph_t> cxG = rcp_const_cast<const xgraph_t>(xG);
    RCP<Zoltan2::XpetraCrsGraphInput<xgraph_t> > xGInput;

    try {
      xGInput =
        rcp(new Zoltan2::XpetraCrsGraphInput<xgraph_t>(cxG));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0,
        string("XpetraCrsGraphInput 3 ")+e.what(), 1);
    }

    if (rank==0){
      std::cout << "Input adapter for Xpetra::CrsGraph" << std::endl;
    }
    fail = verifyInputAdapter<xgraph_t>(*xGInput, *tG);

    gfail = globalFail(comm, fail);

    if (!gfail){
      Zoltan2::PartitioningSolution<gno_t, lno_t, lno_t>
               solution(nprocs, nvtx, 0);
      ArrayRCP<gno_t> &solnGids = solution.getGidsRCP();
      ArrayRCP<size_t> &solnParts = solution.getPartsRCP();
      for (size_t i = 0; i < nvtx; i++) solnGids[i] = rowGids[i];
      memset(solnParts.getRawPtr(), 0, sizeof(size_t) * nvtx);

      xgraph_t *mMigrate =NULL;
      try{
        xGInput->applyPartitioningSolution(*xG, mMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);

      if (!gfail){
        RCP<const xgraph_t> cnewG(mMigrate);
        RCP<Zoltan2::XpetraCrsGraphInput<xgraph_t> > newInput;
        try{
          newInput =
            rcp(new Zoltan2::XpetraCrsGraphInput<xgraph_t>(cnewG));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0,
            string("XpetraCrsGraphInput 4 ")+e.what(), 1);
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

  /////////////////////////////////////////////////////////////
  // User object is Epetra_CrsGraph
  if (!gfail){
    RCP<egraph_t> eG = uinput->getEpetraCrsGraph();
    RCP<const egraph_t> ceG = rcp_const_cast<const egraph_t>(eG);
    RCP<Zoltan2::XpetraCrsGraphInput<egraph_t> > eGInput;

    try {
      eGInput =
        rcp(new Zoltan2::XpetraCrsGraphInput<egraph_t>(ceG));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0,
        string("XpetraCrsGraphInput 5 ")+e.what(), 1);
    }

    if (rank==0){
      std::cout << "Input adapter for Epetra_CrsGraph" << std::endl;
    }
    fail = verifyInputAdapter<egraph_t>(*eGInput, *tG);

    gfail = globalFail(comm, fail);

    if (!gfail){
      Zoltan2::PartitioningSolution<gno_t, lno_t, lno_t>
               solution(nprocs, nvtx, 0);
      ArrayRCP<gno_t> &solnGids = solution.getGidsRCP();
      ArrayRCP<size_t> &solnParts = solution.getPartsRCP();
      for (size_t i = 0; i < nvtx; i++) solnGids[i] = rowGids[i];
      memset(solnParts.getRawPtr(), 0, sizeof(size_t) * nvtx);

      egraph_t *mMigrate =NULL;
      try{
        eGInput->applyPartitioningSolution(*eG, mMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);

      if (!gfail){
        RCP<const egraph_t> cnewG(mMigrate, true);
        RCP<Zoltan2::XpetraCrsGraphInput<egraph_t> > newInput;
        try{
          newInput =
            rcp(new Zoltan2::XpetraCrsGraphInput<egraph_t>(cnewG));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0,
            string("XpetraCrsGraphInput 6 ")+e.what(), 1);
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

  /////////////////////////////////////////////////////////////
  // DONE

  if (rank==0)
    std::cout << "PASS" << std::endl;
}
