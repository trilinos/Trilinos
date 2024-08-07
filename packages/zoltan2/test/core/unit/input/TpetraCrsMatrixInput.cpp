// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Basic testing of Zoltan2::TpetraCrsMatrixAdapter

/*! \file TpetraCrsMatrixInput.cpp
 *  \brief Test of Zoltan2::TpetraCrsMatrixAdapter class.
 *  \todo test with geometric row coordinates.
 */

#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_TpetraRowMatrixAdapter.hpp>
#include <Zoltan2_TpetraCrsMatrixAdapter.hpp>

#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <cstdlib>
#include <stdexcept>

using Teuchos::Comm;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

using ztcrsmatrix_t = Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t>;
using node_t = typename Zoltan2::InputTraits<ztcrsmatrix_t>::node_t;
using device_t = typename node_t::device_type;
using crsAdapter_t = Zoltan2::TpetraCrsMatrixAdapter<ztcrsmatrix_t>;
using execspace_t =
    typename crsAdapter_t::ConstWeightsHostView1D::execution_space;

//////////////////////////////////////////////////////////////////////////

template <typename adapter_t, typename matrix_t>
void TestMatrixIds(adapter_t &ia, matrix_t &matrix) {

  using idsHost_t = typename adapter_t::ConstIdsHostView;
  using offsetsHost_t = typename adapter_t::ConstOffsetsHostView;
  using localInds_t =
      typename adapter_t::user_t::nonconst_local_inds_host_view_type;
  using localVals_t =
      typename adapter_t::user_t::nonconst_values_host_view_type;


  const auto nrows = matrix.getLocalNumRows();
  const auto ncols = matrix.getLocalNumEntries();
  const auto maxNumEntries = matrix.getLocalMaxNumRowEntries();

  typename adapter_t::Base::ConstIdsHostView colIdsHost_("colIdsHost_", ncols);
  typename adapter_t::Base::ConstOffsetsHostView offsHost_("offsHost_",
                                                           nrows + 1);

  localInds_t localColInds("localColInds", maxNumEntries);
  localVals_t localVals("localVals", maxNumEntries);

  for (size_t r = 0; r < nrows; r++) {
    size_t numEntries = 0;
    matrix.getLocalRowCopy(r, localColInds, localVals, numEntries);;

    offsHost_(r + 1) = offsHost_(r) + numEntries;
    for (size_t e = offsHost_(r), i = 0; e < offsHost_(r + 1); e++) {
      colIdsHost_(e) = matrix.getColMap()->getGlobalElement(localColInds(i++));
    }
  }

  idsHost_t rowIdsHost;
  ia.getRowIDsHostView(rowIdsHost);

  const auto matrixIDS = matrix.getRowMap()->getLocalElementList();

  Z2_TEST_COMPARE_ARRAYS(matrixIDS, rowIdsHost);

  idsHost_t colIdsHost;
  offsetsHost_t offsetsHost;
  ia.getCRSHostView(offsetsHost, colIdsHost);

  Z2_TEST_COMPARE_ARRAYS(colIdsHost_, colIdsHost);
  Z2_TEST_COMPARE_ARRAYS(offsHost_, offsetsHost);
}

template <typename adapter_t, typename matrix_t>
void verifyInputAdapter(adapter_t &ia, matrix_t &matrix) {
  using idsDevice_t = typename adapter_t::ConstIdsDeviceView;
  using idsHost_t = typename adapter_t::ConstIdsHostView;
  using weightsDevice_t = typename adapter_t::WeightsDeviceView1D;
  using weightsHost_t = typename adapter_t::WeightsHostView1D;

  const auto nrows = ia.getLocalNumIDs();

  Z2_TEST_EQUALITY(ia.getLocalNumRows(), matrix.getLocalNumRows());
  Z2_TEST_EQUALITY(ia.getLocalNumColumns(), matrix.getLocalNumCols());
  Z2_TEST_EQUALITY(ia.getLocalNumEntries(), matrix.getLocalNumEntries());

  /////////////////////////////////
  //// getRowIdsView
  /////////////////////////////////

  idsDevice_t rowIdsDevice;
  ia.getRowIDsDeviceView(rowIdsDevice);
  idsHost_t rowIdsHost;
  ia.getRowIDsHostView(rowIdsHost);

  Z2_TEST_DEVICE_HOST_VIEWS(rowIdsDevice, rowIdsHost);

  /////////////////////////////////
  //// setRowWeightsDevice
  /////////////////////////////////
  Z2_TEST_THROW(ia.setRowWeightsDevice(
                    typename adapter_t::WeightsDeviceView1D{}, 50),
                std::runtime_error);

  weightsDevice_t wgts0("wgts0", nrows);
  Kokkos::parallel_for(
      nrows, KOKKOS_LAMBDA(const int idx) { wgts0(idx) = idx * 2; });

  Z2_TEST_NOTHROW(ia.setRowWeightsDevice(wgts0, 0));

  // Don't reuse the same View, since we don't copy the values,
  // we just assign the View (increase use count)
  weightsDevice_t wgts1("wgts1", nrows);
  Kokkos::parallel_for(
      nrows, KOKKOS_LAMBDA(const int idx) { wgts1(idx) = idx * 3; });

  Z2_TEST_NOTHROW(ia.setRowWeightsDevice(wgts1, 1));

  /////////////////////////////////
  //// getRowWeightsDevice
  /////////////////////////////////
  {
    weightsDevice_t weightsDevice;
    Z2_TEST_NOTHROW(ia.getRowWeightsDeviceView(weightsDevice, 0));

    weightsHost_t weightsHost;
    Z2_TEST_NOTHROW(ia.getRowWeightsHostView(weightsHost, 0));

    Z2_TEST_DEVICE_HOST_VIEWS(weightsDevice, weightsHost);

    Z2_TEST_DEVICE_HOST_VIEWS(wgts0, weightsHost);
  }
  {
    weightsDevice_t weightsDevice;
    Z2_TEST_NOTHROW(ia.getRowWeightsDeviceView(weightsDevice, 1));

    weightsHost_t weightsHost;
    Z2_TEST_NOTHROW(ia.getRowWeightsHostView(weightsHost, 1));

    Z2_TEST_DEVICE_HOST_VIEWS(weightsDevice, weightsHost);

    Z2_TEST_DEVICE_HOST_VIEWS(wgts1, weightsHost);
  }
  {
    weightsDevice_t wgtsDevice;
    Z2_TEST_THROW(ia.getRowWeightsDeviceView(wgtsDevice, 2),
                  std::runtime_error);

    weightsHost_t wgtsHost;
    Z2_TEST_THROW(ia.getRowWeightsHostView(wgtsHost, 2), std::runtime_error);
  }

  TestMatrixIds(ia, matrix);
}

//////////////////////////////////////////////////////////////////////////

int main(int narg, char *arg[]) {
  using crsSoln_t = Zoltan2::PartitioningSolution<crsAdapter_t>;
  using crsPart_t = crsAdapter_t::part_t;

  Tpetra::ScopeGuard tscope(&narg, &arg);
  const auto comm = Tpetra::getDefaultComm();

  try {
    Teuchos::ParameterList params;
    params.set("input file", "simple");
    params.set("file type", "Chaco");

    auto uinput = rcp(new UserInputForTests(params, comm));

    // Input crs matrix
    const auto crsMatrix = uinput->getUITpetraCrsMatrix();

    const auto nrows = crsMatrix->getLocalNumRows();

    // To test migration in the input adapter we need a Solution object.
    const auto env = rcp(new Zoltan2::Environment(comm));

    const int nWeights = 2;

    /////////////////////////////////////////////////////////////
    // User object is Tpetra::CrsMatrix
    /////////////////////////////////////////////////////////////

    PrintFromRoot("Input adapter for Tpetra::CrsMatrix");

    // Graph Adapters use crsGraph, original TpetraInput uses trM (=CrsMatrix)
    auto tpetraCrsMatrixInput = rcp(new crsAdapter_t(crsMatrix, nWeights));

    verifyInputAdapter(*tpetraCrsMatrixInput, *crsMatrix);

    crsPart_t *p = new crsPart_t[nrows];
    memset(p, 0, sizeof(crsPart_t) * nrows);
    ArrayRCP<crsPart_t> solnParts(p, 0, nrows, true);

    crsSoln_t solution(env, comm, nWeights);
    solution.setParts(solnParts);

    ztcrsmatrix_t *mMigrate = NULL;
    tpetraCrsMatrixInput->applyPartitioningSolution(*crsMatrix, mMigrate,
                                                    solution);
    const auto newM = rcp(mMigrate);
    auto cnewM = rcp_const_cast<const ztcrsmatrix_t>(newM);
    auto newInput = rcp(new crsAdapter_t(cnewM, nWeights));

    PrintFromRoot("Input adapter for Tpetra::CrsMatrix migrated to proc 0");

    verifyInputAdapter(*newInput, *newM);

  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  PrintFromRoot("PASS");

}
