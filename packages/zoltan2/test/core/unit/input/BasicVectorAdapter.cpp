// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Test for Zoltan2::BasicVectorAdapter

#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_TestingHelpers.hpp>

#include <iostream>

typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> userTypes_t;

void testBasisVector(Zoltan2::BasicVectorAdapter<userTypes_t> *ia, int *valueStrides, int wdim, int mvdim) {

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    // Ids
    const zgno_t *ids;
    Zoltan2::BaseAdapter<userTypes_t>::ConstIdsHostView kHostIds;
    Zoltan2::BaseAdapter<userTypes_t>::ConstIdsDeviceView kDeviceIds;

    ia->getIDsView(ids);
    ia->getIDsHostView(kHostIds);
    ia->getIDsDeviceView(kDeviceIds);

    auto kDeviceIdsHost = Kokkos::create_mirror_view(kDeviceIds);
    Kokkos::deep_copy(kDeviceIdsHost, kDeviceIds);

    bool success = true;
    for (size_t i = 0; i < ia->getLocalNumIDs(); ++i) {
      TEUCHOS_TEST_EQUALITY(ids[i], kHostIds(i), std::cout, success);
      TEUCHOS_TEST_EQUALITY(kHostIds(i), kDeviceIdsHost(i), std::cout, success);
    }
    TEST_FAIL_AND_EXIT(*comm, success, "ids != hostIds != deviceIds", 1)

    // Weights
    Zoltan2::BaseAdapter<userTypes_t>::WeightsHostView kHostWgts;
    Zoltan2::BaseAdapter<userTypes_t>::WeightsDeviceView kDeviceWgts;
    ia->getWeightsHostView(kHostWgts);
    ia->getWeightsDeviceView(kDeviceWgts);
    auto kDeviceWgtsHost = Kokkos::create_mirror_view(kDeviceWgts);
    Kokkos::deep_copy(kDeviceWgtsHost, kDeviceWgts);

    for (int w = 0; success && w < wdim; w++) {
      const zscalar_t *wgts;

      Kokkos::View<zscalar_t **, typename znode_t::device_type> wkgts;
      int stride;
      ia->getWeightsView(wgts, stride, w);
      for (zlno_t i = 0; success && i < ia->getNumWeightsPerID(); ++i) {
        TEUCHOS_TEST_EQUALITY(wgts[stride * i], kHostWgts(i, w), std::cout,
                              success);
        TEUCHOS_TEST_EQUALITY(wgts[stride * i], kDeviceWgtsHost(i, w), std::cout,
                              success);
      }
    }
    TEST_FAIL_AND_EXIT(*comm, success, "wgts != vwgts != wkgts", 1)

    // Coords
    Zoltan2::AdapterWithCoords<userTypes_t>::CoordsHostView kHostCoords;
    Zoltan2::AdapterWithCoords<userTypes_t>::CoordsDeviceView kDeviceCoords;
    ia->getCoordinatesHostView(kHostCoords);
    ia->getCoordinatesDeviceView(kDeviceCoords);
    auto kHostCoordsMV = Kokkos::create_mirror_view(kHostCoords);
    Kokkos::deep_copy(kHostCoordsMV, kHostCoords);

    auto kDeviceCoordsMV = Kokkos::create_mirror_view(kDeviceCoords);
    Kokkos::deep_copy(kDeviceCoordsMV, kDeviceCoords);
    for (int v=0; v < mvdim; v++){
        const zscalar_t *coords;
        int stride;

        ia->getCoordinatesView(coords, stride, v);

        success = true;
        for (size_t i = 0; i < ia->getLocalNumIDs(); ++i) {
          TEUCHOS_TEST_EQUALITY(coords[i*stride], kHostCoordsMV(i, v), std::cout, success);
          TEUCHOS_TEST_EQUALITY(coords[i*stride], kDeviceCoordsMV(i, v), std::cout, success);
        }
        TEST_FAIL_AND_EXIT(*comm, success, "ids != kHostCoordsMV != kDeviceCoordsMV", 1)
    }
}

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0;

  // Create a single vector and a strided multi-vector with
  // strided multi-weights.

  zlno_t numLocalIds = 10;
  zgno_t *myIds = new zgno_t[numLocalIds];
  zgno_t myFirstId = rank * numLocalIds;

  int wdim = 2;
  zscalar_t *weights = new zscalar_t [numLocalIds*wdim];
  int *weightStrides = new int [wdim];
  const zscalar_t **weightPtrs = new const zscalar_t * [wdim];

  int mvdim = 3;
  zscalar_t *v_values= new zscalar_t [numLocalIds];
  zscalar_t *mv_values= new zscalar_t [mvdim * numLocalIds];
  int *valueStrides = new int [mvdim];
  const zscalar_t **valuePtrs = new const zscalar_t * [mvdim];

  for (zlno_t i=0; i < numLocalIds; i++){
    myIds[i] = myFirstId+i;

    for (int w=0; w < wdim; w++)
      weights[w*numLocalIds + i] = w + 1 + nprocs - rank;

    v_values[i] = numLocalIds-i;

    for (int v=0; v < mvdim; v++)
      mv_values[i*mvdim + v] = (v+1) * (nprocs-rank) / (i+1);
  }

  for (int w=0; w < wdim; w++){
    weightStrides[w] = 1;
    weightPtrs[w] = weights + numLocalIds*w;
  }

  for (int v=0; v < mvdim; v++){
    valueStrides[v] = mvdim;
    valuePtrs[v] = mv_values + v;
  }

  Zoltan2::BasicVectorAdapter<userTypes_t> *ia = NULL;

  {
    // A Zoltan2::BasicVectorAdapter object with one vector and weights

    std::vector<const zscalar_t *> weightValues;
    std::vector<int> strides;

    weightValues.push_back(weightPtrs[0]);
    weightValues.push_back(weightPtrs[1]);
    strides.push_back(1);
    strides.push_back(1);

    try{
     ia = new Zoltan2::BasicVectorAdapter<userTypes_t>(numLocalIds, myIds,
       v_values, 1, true, weightPtrs[0], 1);
    }
    catch (std::exception &e){
      fail = 1;
    }

    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 1", fail);

    testBasisVector(ia, valueStrides, 1, 1);

    delete ia;


  }

  {
    // A Zoltan2::BasicVectorAdapter object with a multivector with weights

    std::vector<const zscalar_t *> weightValues, values;
    std::vector<int> wstrides, vstrides;

    for (int dim=0; dim < wdim; dim++){
      weightValues.push_back(weightPtrs[dim]);
      wstrides.push_back(1);
    }

    for (int dim=0; dim < mvdim; dim++){
      values.push_back(valuePtrs[dim]);
      vstrides.push_back(mvdim);
    }

    try{
     ia = new Zoltan2::BasicVectorAdapter<userTypes_t>(
        numLocalIds, myIds, values, vstrides, weightValues, wstrides);
    }
    catch (std::exception &e){
      fail = 1;
    }

    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 2", fail);

    testBasisVector(ia, valueStrides, wdim, mvdim);

    delete ia;

  }

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  delete [] myIds;
  delete [] weights;
  delete [] weightStrides;
  delete [] weightPtrs;
  delete [] v_values;
  delete [] mv_values;
  delete [] valueStrides;
  delete [] valuePtrs;

  return fail;
}

