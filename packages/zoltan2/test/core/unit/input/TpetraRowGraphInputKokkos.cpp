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

#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_TpetraCrsGraphAdapter.hpp>
#include <Zoltan2_TpetraRowGraphAdapter.hpp>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <cstdlib>
#include <stdexcept>

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

using ztcrsgraph_t = Tpetra::CrsGraph<zlno_t, zgno_t, znode_t>;
using ztrowgraph_t = Tpetra::RowGraph<zlno_t, zgno_t, znode_t>;
using node_t = typename Zoltan2::InputTraits<ztrowgraph_t>::node_t;
using device_t = typename node_t::device_type;
using rowAdapter_t = Zoltan2::TpetraRowGraphAdapter<ztrowgraph_t>;
using crsAdapter_t = Zoltan2::TpetraCrsGraphAdapter<ztcrsgraph_t>;
using execspace_t =
    typename rowAdapter_t::ConstWeightsHostView1D::execution_space;

template <typename offset_t>
void printGraph(RCP<const Comm<int>> &comm, zlno_t nvtx, const zgno_t *vtxIds,
                const offset_t *offsets, const zgno_t *edgeIds) {
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p = 0; p < nprocs; p++) {
    if (p == rank) {
      std::cout << rank << ":" << std::endl;
      for (zlno_t i = 0; i < nvtx; i++) {
        std::cout << " vertex " << vtxIds[i] << ": ";
        for (offset_t j = offsets[i]; j < offsets[i + 1]; j++) {
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

template <typename adapter_t, typename graph_t>
void TestGraphIds(adapter_t &ia, graph_t &graph) {

  using idsHost_t = typename adapter_t::ConstIdsHostView;
  using offsetsHost_t = typename adapter_t::ConstOffsetsHostView;
  using localInds_t =
      typename adapter_t::user_t::nonconst_local_inds_host_view_type;

  const auto nvtx = graph.getLocalNumRows();
  const auto nedges = graph.getLocalNumEntries();
  const auto maxNumEntries = graph.getLocalMaxNumRowEntries();

  typename adapter_t::Base::ConstIdsHostView adjIdsHost_("adjIdsHost_", nedges);
  typename adapter_t::Base::ConstOffsetsHostView offsHost_("offsHost_",
                                                           nvtx + 1);

  localInds_t nbors("nbors", maxNumEntries);

  for (size_t v = 0; v < nvtx; v++) {
    size_t numColInds = 0;
    graph.getLocalRowCopy(v, nbors, numColInds);

    offsHost_(v + 1) = offsHost_(v) + numColInds;
    for (size_t e = offsHost_(v), i = 0; e < offsHost_(v + 1); e++) {
      adjIdsHost_(e) = graph.getColMap()->getGlobalElement(nbors(i++));
    }
  }

  idsHost_t vtxIdsHost;
  ia.getVertexIDsHostView(vtxIdsHost);

  const auto graphIDS = graph.getRowMap()->getLocalElementList();

  Z2_TEST_COMPARE_ARRAYS(graphIDS, vtxIdsHost);

  idsHost_t adjIdsHost;
  offsetsHost_t offsetsHost;
  ia.getEdgesHostView(offsetsHost, adjIdsHost);

  Z2_TEST_COMPARE_ARRAYS(adjIdsHost_, adjIdsHost);
  Z2_TEST_COMPARE_ARRAYS(offsHost_, offsetsHost);
}

template <typename adapter_t, typename graph_t>
void verifyInputAdapter(adapter_t &ia, graph_t &graph) {
  using idsDevice_t = typename adapter_t::ConstIdsDeviceView;
  using idsHost_t = typename adapter_t::ConstIdsHostView;
  using offsetsDevice_t = typename adapter_t::ConstOffsetsDeviceView;
  using offsetsHost_t = typename adapter_t::ConstOffsetsHostView;
  using weightsDevice_t = typename adapter_t::WeightsDeviceView1D;
  using weightsHost_t = typename adapter_t::WeightsHostView1D;

  const auto nVtx = ia.getLocalNumIDs();

  Z2_TEST_EQUALITY(ia.getLocalNumVertices(), graph.getLocalNumRows());
  Z2_TEST_EQUALITY(ia.getLocalNumEdges(), graph.getLocalNumEntries());

  /////////////////////////////////
  //// getVertexIdsView
  /////////////////////////////////

  idsDevice_t vtxIdsDevice;
  ia.getVertexIDsDeviceView(vtxIdsDevice);
  idsHost_t vtxIdsHost;
  ia.getVertexIDsHostView(vtxIdsHost);

  Z2_TEST_DEVICE_HOST_VIEWS(vtxIdsDevice, vtxIdsHost);

  /////////////////////////////////
  //// getEdgesView
  /////////////////////////////////

  idsDevice_t adjIdsDevice;
  offsetsDevice_t offsetsDevice;

  ia.getEdgesDeviceView(offsetsDevice, adjIdsDevice);

  idsHost_t adjIdsHost;
  offsetsHost_t offsetsHost;
  ia.getEdgesHostView(offsetsHost, adjIdsHost);

  Z2_TEST_DEVICE_HOST_VIEWS(adjIdsDevice, adjIdsHost);
  Z2_TEST_DEVICE_HOST_VIEWS(offsetsDevice, offsetsHost);

  /////////////////////////////////
  //// setVertexWeightsDevice
  /////////////////////////////////
  Z2_TEST_THROW(ia.setVertexWeightsDevice(
                    typename adapter_t::ConstWeightsDeviceView1D{}, 50),
                std::runtime_error);

  weightsDevice_t wgts0("wgts0", nVtx);
  Kokkos::parallel_for(
      nVtx, KOKKOS_LAMBDA(const int idx) { wgts0(idx) = idx * 2; });
  Kokkos::fence();

  Z2_TEST_NOTHROW(ia.setVertexWeightsDevice(wgts0, 0));

  // Don't reuse the same View, since we don't copy the values,
  // we just assign the View (increase use count)
  weightsDevice_t wgts1("wgts1", nVtx);
  Kokkos::parallel_for(
      nVtx, KOKKOS_LAMBDA(const int idx) { wgts1(idx) = idx * 3; });

  Z2_TEST_NOTHROW(ia.setVertexWeightsDevice(wgts1, 1));

  /////////////////////////////////
  //// getVertexWeightsDevice
  /////////////////////////////////
  {
    weightsDevice_t weightsDevice;
    Z2_TEST_NOTHROW(ia.getVertexWeightsDeviceView(weightsDevice, 0));

    weightsHost_t weightsHost;
    Z2_TEST_NOTHROW(ia.getVertexWeightsHostView(weightsHost, 0));

    Z2_TEST_DEVICE_HOST_VIEWS(weightsDevice, weightsHost);

    Z2_TEST_DEVICE_HOST_VIEWS(wgts0, weightsHost);
  }
  {
    weightsDevice_t weightsDevice;
    Z2_TEST_NOTHROW(ia.getVertexWeightsDeviceView(weightsDevice, 1));

    weightsHost_t weightsHost;
    Z2_TEST_NOTHROW(ia.getVertexWeightsHostView(weightsHost, 1));

    Z2_TEST_DEVICE_HOST_VIEWS(weightsDevice, weightsHost);

    Z2_TEST_DEVICE_HOST_VIEWS(wgts1, weightsHost);
  }
  {
    weightsDevice_t wgtsDevice;
    Z2_TEST_THROW(ia.getVertexWeightsDeviceView(wgtsDevice, 2),
                  std::runtime_error);

    weightsHost_t wgtsHost;
    Z2_TEST_THROW(ia.getVertexWeightsHostView(wgtsHost, 2), std::runtime_error);
  }

  TestGraphIds(ia, graph);
}

int main(int narg, char *arg[]) {
  using rowSoln_t = Zoltan2::PartitioningSolution<rowAdapter_t>;
  using rowPart_t = rowAdapter_t::part_t;

  using crsSoln_t = Zoltan2::PartitioningSolution<crsAdapter_t>;
  using crsPart_t = crsAdapter_t::part_t;

  Tpetra::ScopeGuard tscope(&narg, &arg);
  const auto comm = Tpetra::getDefaultComm();

  try {
    Teuchos::ParameterList params;
    params.set("input file", "simple");
    params.set("file type", "Chaco");

    auto uinput = rcp(new UserInputForTests(params, comm));

    // Input crs graph and row graph cast from it.
    const auto crsGraph = uinput->getUITpetraCrsGraph();
    const auto rowGraph = rcp_dynamic_cast<ztrowgraph_t>(crsGraph);

    const auto nvtx = rowGraph->getLocalNumRows();

    // To test migration in the input adapter we need a Solution object.
    const auto env = rcp(new Zoltan2::Environment(comm));

    const int nWeights = 2;

    /////////////////////////////////////////////////////////////
    // User object is Tpetra::CrsGraph
    /////////////////////////////////////////////////////////////
    {
      PrintFromRoot("Input adapter for Tpetra::CrsGraph");

      auto tpetraCrsGraphInput = rcp(new crsAdapter_t(crsGraph, nWeights));

      verifyInputAdapter(*tpetraCrsGraphInput, *crsGraph);

      ztcrsgraph_t *mMigrate = NULL;
      crsPart_t *p = new crsPart_t[nvtx];
      memset(p, 0, sizeof(crsPart_t) * nvtx);
      ArrayRCP<crsPart_t> solnParts(p, 0, nvtx, true);

      crsSoln_t solution(env, comm, nWeights);
      solution.setParts(solnParts);
      tpetraCrsGraphInput->applyPartitioningSolution(*crsGraph, mMigrate,
                                                     solution);
      const auto newG = rcp(mMigrate);

      auto cnewG = rcp_const_cast<const ztcrsgraph_t>(newG);
      auto newInput = rcp(new crsAdapter_t(cnewG, nWeights));

      PrintFromRoot("Input adapter for Tpetra::RowGraph migrated to proc 0");

      verifyInputAdapter(*newInput, *newG);
    }

    /////////////////////////////////////////////////////////////
    // User object is Tpetra::RowGraph
    /////////////////////////////////////////////////////////////
    {
      PrintFromRoot("Input adapter for Tpetra::RowGraph");

      auto tpetraRowGraphInput = rcp(new rowAdapter_t(rowGraph, nWeights));

      verifyInputAdapter(*tpetraRowGraphInput, *crsGraph);

      rowPart_t *p = new rowPart_t[nvtx];
      memset(p, 0, sizeof(rowPart_t) * nvtx);
      ArrayRCP<rowPart_t> solnParts(p, 0, nvtx, true);

      rowSoln_t solution(env, comm, nWeights);
      solution.setParts(solnParts);

      ztrowgraph_t *mMigrate = NULL;
      tpetraRowGraphInput->applyPartitioningSolution(*crsGraph, mMigrate,
                                                     solution);
      const auto newG = rcp(mMigrate);

      auto cnewG = rcp_const_cast<const ztrowgraph_t>(newG);
      auto newInput = rcp(new rowAdapter_t(cnewG, nWeights));

      PrintFromRoot("Input adapter for Tpetra::RowGraph migrated to proc 0");

      verifyInputAdapter(*newInput, *newG);
    }
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  PrintFromRoot("PASS");
}
