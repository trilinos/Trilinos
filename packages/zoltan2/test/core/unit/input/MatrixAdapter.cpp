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
// Testing of GraphAdapter built from Xpetra matrix input adapters.
//

/*! \brief Test of GraphAdapter interface.
 *
 */

#include "Kokkos_StaticCrsGraph.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_TpetraRowMatrixAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>

#include <bitset>
#include <iostream>
#include <string>

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_TestingHelpers.hpp>

using Teuchos::ArrayView;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

using simpleUser_t = Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t>;

using tcrsMatrix_t = Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t>;
using trowMatrix_t = Tpetra::RowMatrix<zscalar_t, zlno_t, zgno_t, znode_t>;
using tmap_t = Tpetra::Map<zlno_t, zgno_t, znode_t>;

using simpleVAdapter_t = Zoltan2::BasicVectorAdapter<simpleUser_t>;
using baseMAdapter_t = Zoltan2::MatrixAdapter<tcrsMatrix_t, simpleUser_t>;
using tRowMAdapter_t = typename Zoltan2::TpetraRowMatrixAdapter<trowMatrix_t, simpleUser_t>;
using zoffset_t = typename baseMAdapter_t::offset_t;

/////////////////////////////////////////////////////////////////////////////
int main(int narg, char *arg[]) {
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();

  int nVtxWeights = 5;
  int nnzWgtIdx = -1;
  std::string fname("simple");

  if (rank == 0)
    std::cout << "TESTING base case (global)" << std::endl;

  // Input generator
  UserInputForTests *uinput;

  uinput = new UserInputForTests(testDataFilePath, fname, comm, true);

  RCP<tcrsMatrix_t> M = uinput->getUITpetraCrsMatrix();
  RCP<trowMatrix_t> trM = rcp_dynamic_cast<trowMatrix_t>(M);
  RCP<const trowMatrix_t> ctrM = rcp_const_cast<const trowMatrix_t>(trM);
  zlno_t nLocalRows = M->getLocalNumRows();

  // Weights:
  zscalar_t **rowWeights = NULL;
  Zoltan2::BaseAdapter<trowMatrix_t>::WeightsDeviceView kRowWeights(
      Kokkos::ViewAllocateWithoutInitializing("kRowWeights"), nVtxWeights,
      nLocalRows);

  auto kRowWeightsHost = Kokkos::create_mirror_view(kRowWeights);

  if (nVtxWeights > 0) {
    rowWeights = new zscalar_t *[nVtxWeights];
    for (int i = 0; i < nVtxWeights; i++) {
      if (nnzWgtIdx == i) {
        rowWeights[i] = NULL;
        kRowWeightsHost(i, 0) = 0;
      } else {
        rowWeights[i] = new zscalar_t[nLocalRows];
        for (zlno_t j = 0; j < nLocalRows; j++) {
          rowWeights[i][j] = 200000 * i + j;
          kRowWeightsHost(i, j) = 200000 * i + j;
        }
      }
    }
  }

  Kokkos::deep_copy(kRowWeights, kRowWeightsHost);

  tRowMAdapter_t tmi(ctrM, nVtxWeights);
  for (int i = 0; i < nVtxWeights; i++) {
    tmi.setWeights(rowWeights[i], 1, i);
  }
  tmi.setRowWeights(kRowWeights);

  simpleVAdapter_t *via = NULL;

  // Set up some fake input
  zscalar_t **coords = NULL;
  int coordDim = 3;

  if (coordDim > 0) {
    coords = new zscalar_t *[coordDim];
    for (int i = 0; i < coordDim; i++) {
      coords[i] = new zscalar_t[nLocalRows];
      for (zlno_t j = 0; j < nLocalRows; j++) {
        coords[i][j] = 100000 * i + j;
      }
    }
  }

  zgno_t *gids = NULL;
  if (coordDim > 0) {
    gids = new zgno_t[nLocalRows];
    for (zlno_t i = 0; i < nLocalRows; i++)
      gids[i] = M->getRowMap()->getGlobalElement(i);
    via = new simpleVAdapter_t(nLocalRows, gids, coords[0],
                               (coordDim > 1 ? coords[1] : NULL),
                               (coordDim > 2 ? coords[2] : NULL), 1, 1, 1);
    tmi.setCoordinateInput(via);
  }

  // TEST of getIDsView, getIDsKokkosView, getRowIDsView, getRowIDsHostView and
  // getRowIDsDeviceView

  const zgno_t *ids;
  const zgno_t *rowIds;
  Zoltan2::BaseAdapter<simpleUser_t>::ConstIdsDeviceView kIds;
  Zoltan2::BaseAdapter<simpleUser_t>::ConstIdsHostView kHostIds;
  Zoltan2::BaseAdapter<simpleUser_t>::ConstIdsDeviceView kDeviceIds;

  tmi.getIDsView(ids);
  tmi.getRowIDsView(rowIds);
  tmi.getIDsKokkosView(kIds);
  auto kIdsHost = Kokkos::create_mirror_view(kIds);
  Kokkos::deep_copy(kIdsHost, kIds);

  tmi.getRowIDsHostView(kHostIds);
  tmi.getRowIDsDeviceView(kDeviceIds);

  auto kDeviceIdsHost = Kokkos::create_mirror_view(kDeviceIds);
  Kokkos::deep_copy(kDeviceIdsHost, kDeviceIds);

  bool success = true;
  for (size_t i = 0; i < tmi.getLocalNumIDs(); ++i) {
    TEUCHOS_TEST_EQUALITY(ids[i], kIdsHost(i), std::cout, success);
    TEUCHOS_TEST_EQUALITY(ids[i], rowIds[i], std::cout, success);
    TEUCHOS_TEST_EQUALITY(kIdsHost(i), kHostIds(i), std::cout, success);
    TEUCHOS_TEST_EQUALITY(kIdsHost(i), kDeviceIdsHost(i), std::cout, success);
  }
  TEST_FAIL_AND_EXIT(*comm, success, "ids != vertexIds != kIds", 1)

  // TEST of getCRSView, getCRSHostView and getCRSDeviceView
  //  const zoffset_t *offsets;
  ArrayRCP<const zoffset_t> offsets;
  ArrayRCP<const zgno_t> colIds;
  Zoltan2::BaseAdapter<trowMatrix_t>::ConstOffsetsHostView kHostOffsets;
  Zoltan2::BaseAdapter<trowMatrix_t>::ConstIdsHostView kHostColIds;
  Zoltan2::BaseAdapter<trowMatrix_t>::ConstOffsetsDeviceView kDeviceOffsets;
  Zoltan2::BaseAdapter<trowMatrix_t>::ConstIdsDeviceView kDeviceColIds;
  tmi.getCRSView(offsets, colIds);
  tmi.getCRSHostView(kHostOffsets, kHostColIds);
  tmi.getCRSDeviceView(kDeviceOffsets, kDeviceColIds);

  auto kDeviceColIdsHost = Kokkos::create_mirror_view(kDeviceColIds);
  Kokkos::deep_copy(kDeviceColIdsHost, kDeviceColIds);

  auto kDeviceOffsetsHost = Kokkos::create_mirror_view(kDeviceOffsets);
  Kokkos::deep_copy(kDeviceOffsetsHost, kDeviceOffsets);

  for (int i = 0; success && i < colIds.size(); i++) {
    TEUCHOS_TEST_EQUALITY(colIds[i], kHostColIds(i), std::cout, success);
    TEUCHOS_TEST_EQUALITY(colIds[i], kDeviceColIdsHost(i), std::cout, success);
  }
  TEST_FAIL_AND_EXIT(*comm, success, "colIds != kHostColIds != kDeviceColIds",
                     1)

  for (int i = 0; success && i < offsets.size(); i++) {
    TEUCHOS_TEST_EQUALITY(offsets[i], kHostOffsets(i), std::cout, success);
    TEUCHOS_TEST_EQUALITY(offsets[i], kDeviceOffsetsHost(i), std::cout,
                          success);
  }
  TEST_FAIL_AND_EXIT(*comm, success,
                     "offsets != kHostOffsets != kDeviceOffsets", 1)

  ArrayRCP<const zscalar_t> values;
  Zoltan2::BaseAdapter<trowMatrix_t>::ConstScalarsHostView kHostValues;
  Zoltan2::BaseAdapter<trowMatrix_t>::ConstScalarsDeviceView kDeviceValues;

  tmi.getCRSView(offsets, colIds, values);
  tmi.getCRSHostView(kHostOffsets, kHostColIds, kHostValues);
  tmi.getCRSDeviceView(kDeviceOffsets, kDeviceColIds, kDeviceValues);
  auto kDeviceValuesHost = Kokkos::create_mirror_view(kDeviceValues);
  Kokkos::deep_copy(kDeviceValuesHost, kDeviceValues);

  for (int i = 0; success && i < colIds.size(); i++) {
    TEUCHOS_TEST_EQUALITY(colIds[i], kHostColIds(i), std::cout, success);
    TEUCHOS_TEST_EQUALITY(colIds[i], kDeviceColIdsHost(i), std::cout, success);
  }
  TEST_FAIL_AND_EXIT(*comm, success, "colIds != kHostColIds != kDeviceColIds",
                     1)

  for (int i = 0; success && i < offsets.size(); i++) {
    TEUCHOS_TEST_EQUALITY(offsets[i], kHostOffsets(i), std::cout, success);
    TEUCHOS_TEST_EQUALITY(offsets[i], kDeviceOffsetsHost(i), std::cout,
                          success);
  }
  TEST_FAIL_AND_EXIT(*comm, success,
                     "offsets != kHostOffsets != kDeviceOffsets", 1)

  for (int i = 0; success && i < values.size(); i++) {
    TEUCHOS_TEST_EQUALITY(values[i], kHostValues(i), std::cout, success);
    TEUCHOS_TEST_EQUALITY(values[i], kDeviceValuesHost(i), std::cout, success);
  }
  TEST_FAIL_AND_EXIT(*comm, success, "values != kHostValues != kDeviceValues",
                     1)

  // TEST of getRowWeightsView, getRowWeightsHostView and
  // getRowWeightsDeviceView
  for (int w = 0; success && w < nVtxWeights; w++) {
    const zscalar_t *wgts;
    Zoltan2::BaseAdapter<trowMatrix_t>::WeightsHostView kHostWgts;
    Zoltan2::BaseAdapter<trowMatrix_t>::WeightsDeviceView kDeviceWgts;

    Kokkos::View<zscalar_t **, typename znode_t::device_type> wkgts;
    int stride;
    tmi.getRowWeightsView(wgts, stride, w);
    tmi.getRowWeightsHostView(kHostWgts);
    tmi.getRowWeightsDeviceView(kDeviceWgts);

    auto kDeviceWgtsHost = Kokkos::create_mirror_view(kDeviceWgts);
    Kokkos::deep_copy(kDeviceWgtsHost, kDeviceWgts);

    for (zlno_t i = 0; success && i < nLocalRows; ++i) {
      TEUCHOS_TEST_EQUALITY(wgts[stride * i], kHostWgts(w, i), std::cout,
                            success);
      TEUCHOS_TEST_EQUALITY(wgts[stride * i], kDeviceWgtsHost(w, i), std::cout,
                            success);
    }
  }
  TEST_FAIL_AND_EXIT(*comm, success, "wgts != vwgts != wkgts", 1)

  // TEST of useNumNonzerosAsRowWeight, useNumNonzerosAsRowWeightHost and
  // useNumNonzerosAsRowWeightDevice
  // for (int w = 0; success && w < nVtxWeights; w++) {
  //   TEUCHOS_TEST_EQUALITY(tmi.useNumNonzerosAsRowWeight(w),
  //                         tmi.useNumNonzerosAsRowWeightHost(w), std::cout,
  //                         success);
  //   TEUCHOS_TEST_EQUALITY(tmi.useNumNonzerosAsRowWeight(w),
  //                         tmi.useNumNonzerosAsRowWeightDevice(w), std::cout,
  //                         success);
  // }
  // TEST_FAIL_AND_EXIT(
  //     *comm, success,
  //     "useNumNonzerosAsRowWeight != useNumNonzerosAsRowWeightHost != "
  //     "useNumNonzerosAsRowWeightDevice",
  //     1)

  // clean
  if (nVtxWeights > 0) {
    for (int i = 0; i < nVtxWeights; i++) {
      if (rowWeights[i])
        delete[] rowWeights[i];
    }
    delete[] rowWeights;
  }

  if (coordDim > 0) {
    delete via;
    delete[] gids;
    for (int i = 0; i < coordDim; i++) {
      if (coords[i])
        delete[] coords[i];
    }
    delete[] coords;
  }

  delete uinput;

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  return 0;
}
