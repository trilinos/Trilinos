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
// Test for Zoltan2::BasicVectorAdapter for coordinate-based problems

#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::Array;

typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> userTypes_t;

int checkBasicCoordinate(
  Zoltan2::BasicVectorAdapter<userTypes_t> *ia, 
  int len, int glen, zgno_t *ids,
  zscalar_t *xyz,
  zscalar_t *weights,
  int nCoords, int nWeights)
{
  int fail = 0;

  if (ia->getNumEntriesPerID() != nCoords)
    fail = 100;

  if (!fail && ia->getNumWeightsPerID() != nWeights)
    fail = 101;

  if (!fail && ia->getLocalNumIDs() != size_t(len))
    fail = 102;

  for (int x=0; !fail && x < nCoords; x++){
    const zgno_t *idList;
    const zscalar_t *vals;
    int stride;

    ia->getIDsView(idList);
    ia->getEntriesView(vals, stride, x);

    zscalar_t *coordVal = xyz + x;
    for (int i=0; !fail && i < len; i++, coordVal += 3){

      if (idList[i] != ids[i])
        fail = 105;

      if (!fail && vals[stride*i] != *coordVal)
        fail = 106;
    }
  }

  for (int w=0; !fail && w < nWeights; w++){
    const zscalar_t *wgts;
    int stride;

    ia->getWeightsView(wgts, stride, w);

    zscalar_t *weightVal = weights + len*w;
    for (int i=0; !fail && i < len; i++, weightVal++){
      if (wgts[stride*i] != *weightVal)
        fail = 110;
    }
  }

  return fail;
}
  

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int fail = 0;

  // Get some coordinates

  typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> mv_t;
  RCP<UserInputForTests> uinput;
  Teuchos::ParameterList params;
  params.set("input file", "simple");
  params.set("file type", "Chaco");

  try{
    uinput = rcp(new UserInputForTests(params, comm));
  }
  catch(std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "input constructor", 1);

  RCP<mv_t> coords;

  try{
    coords = uinput->getUICoordinates();
  }
  catch(std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "getting coordinates", 1);

  int numLocalIds = coords->getLocalLength();
  int numGlobalIds = coords->getGlobalLength();
  int coordDim = coords->getNumVectors();
  ArrayView<const zgno_t> idList = coords->getMap()->getNodeElementList();

  // Create global Ids, x-, y- and z-coordinates, and also arrays of weights.

  Array<zgno_t> myIds(numLocalIds);
  zgno_t myFirstId = rank * numLocalIds;
  
  int wdim = 2;
  Array<zscalar_t> weights(numLocalIds*wdim);
  for (int i = 0; i < numLocalIds*wdim; i++) weights[i] = zscalar_t(i);

  zscalar_t *x_values= coords->getDataNonConst(0).getRawPtr();
  zscalar_t *y_values= x_values;  // fake 3 dimensions if needed
  zscalar_t *z_values= x_values;

  if (coordDim > 1){
    y_values= coords->getDataNonConst(1).getRawPtr();
    if (coordDim > 2)
      z_values= coords->getDataNonConst(2).getRawPtr();
  }

  Array<zscalar_t> xyz_values(3*numLocalIds);

  for (zlno_t i=0; i < numLocalIds; i++)   // global Ids
    myIds[i] = myFirstId+i;

  zscalar_t *x = xyz_values.getRawPtr();   // a stride-3 coordinate array
  zscalar_t *y = x+1;
  zscalar_t *z = y+1;

  for (int i=0, ii=0; i < numLocalIds; i++, ii += 3){
    x[ii] = x_values[i];
    y[ii] = y_values[i];
    z[ii] = z_values[i];
  }

  RCP<Zoltan2::BasicVectorAdapter<userTypes_t> > ia;

  {
    ////////////////////////////////////////////////////////////////
    // 3-dimensional coordinates with stride one and no weights,
    //   using simpler constructor
  
    int ncoords = 3;
    int nweights = 0;
  
    try{
     ia = rcp(new Zoltan2::BasicVectorAdapter<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), x_values, y_values, z_values));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 0", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 0", fail);
  }

  {
    ////////////////////////////////////////////////////////////////
    // 3-dimensional coordinates with stride one and one weight
    //   using simpler constructor
  
    int ncoords = 3;
    int nweights = 1;
  
    try{
     ia = rcp(new Zoltan2::BasicVectorAdapter<userTypes_t>(
       numLocalIds, myIds.getRawPtr(),
       x_values, y_values, z_values, 1, 1, 1,
       true, weights.getRawPtr(), 1));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 0a", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 0a", fail);
  }

  {
    ////////////////////////////////////////////////////////////////
    // 3-dimensional coordinates with stride one and no weights
  
    int ncoords = 3;
    int nweights = 0;

    std::vector<const zscalar_t *> values, weightValues;
    std::vector<int> valueStrides, weightStrides;
  
    values.push_back(x_values);
    values.push_back(y_values);
    values.push_back(z_values);
    valueStrides.push_back(1);
    valueStrides.push_back(1);
    valueStrides.push_back(1);
  
    try{
     ia = rcp(new Zoltan2::BasicVectorAdapter<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), values, valueStrides, 
       weightValues, weightStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 1", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 1", fail);
  
    // Try using the default: no strides supplied means strides are one.

    std::vector<int> emptyStrides;
  
    try{
     ia = rcp(new Zoltan2::BasicVectorAdapter<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), values, emptyStrides, 
       weightValues, emptyStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 2", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 2", fail);
  }

  {
    ////////////////////////////////////////////////////////////////
    // 2-dimensional coordinates with stride three and two weights
  
    int ncoords = 2;
    int nweights = 2;

    std::vector<const zscalar_t *> values, weightValues;
    std::vector<int> valueStrides, weightStrides;
  
    values.push_back(xyz_values.getRawPtr());
    values.push_back(xyz_values.getRawPtr() + 1);
    valueStrides.push_back(3);
    valueStrides.push_back(3);
  
    weightValues.push_back(weights.getRawPtr());
    weightValues.push_back(weights.getRawPtr() + numLocalIds);
    weightStrides.push_back(1);
    weightStrides.push_back(1);
  
    try{
     ia = rcp(new Zoltan2::BasicVectorAdapter<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), values, valueStrides, 
       weightValues, weightStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 3", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 3", fail);
  
    // Try using default weight strides

    std::vector<int> emptyStrides;
  
    try{
     ia = rcp(new Zoltan2::BasicVectorAdapter<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), values, valueStrides, 
       weightValues, emptyStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 4", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 4", fail);
  }

  {
    ////////////////////////////////////////////////////////////////
    // 1-dimensional coordinates with stride one and two weights

    int ncoords = 1;
    int nweights = 2;

    std::vector<const zscalar_t *> values, weightValues;
    std::vector<int> valueStrides, weightStrides;
  
    values.push_back(x_values);
    valueStrides.push_back(1);
  
    weightValues.push_back(weights.getRawPtr());
    weightValues.push_back(weights.getRawPtr() + numLocalIds);
    weightStrides.push_back(1);
    weightStrides.push_back(1);
  
    try{
     ia = rcp(new Zoltan2::BasicVectorAdapter<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), values, valueStrides, 
       weightValues, weightStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 4", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 4", fail);
  }

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  return fail;
}

