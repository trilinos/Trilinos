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
// Test for Zoltan2::BasicVectorAdapter 

#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>


typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> userTypes_t;

int checkBasicVector(
  Zoltan2::BasicVectorAdapter<userTypes_t> *ia, int len, int glen,
  zgno_t *ids, int mvdim, const zscalar_t **values, int *valueStrides,
  int wdim, const zscalar_t **weights, int *weightStrides)
{
  int fail = 0;
  bool strideOne = false;

  if (valueStrides == NULL) strideOne = true;

  if (ia->getNumEntriesPerID() != mvdim)
    fail = 100;

  if (!fail && ia->getNumWeightsPerID() != wdim)
    fail = 101;

  if (!fail && ia->getLocalNumIDs() != size_t(len))
    fail = 102;

  const zgno_t *idList;
  ia->getIDsView(idList);
  for (int i=0; !fail && i < len; i++)
    if (!fail && idList[i] != ids[i])
      fail = 107;

  for (int v=0; !fail && v < mvdim; v++){
    const zscalar_t *vals;
    int correctStride = (strideOne ? 1 : valueStrides[v]);
    int stride;

    ia->getEntriesView(vals, stride, v);

    if (!fail && stride != correctStride)
      fail = 105;

    for (int i=0; !fail && i < len; i++){
// TODO fix values check
//      if (vals[stride*i] != values[v][correctStride*i])
//        fail = 106;
    }
  }

  for (int w=0; !fail && w < wdim; w++){
    const zscalar_t *wgts;
    int stride;

    ia->getWeightsView(wgts, stride, w);

    if (!fail && stride != weightStrides[w])
      fail = 109;

    for (int i=0; !fail && i < len; i++){
      if (wgts[stride*i] != weights[w][weightStrides[w]*i])
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
    // A Zoltan2::BasicVectorAdapter object with one vector and no weights

    std::vector<const zscalar_t *> weightValues;
    std::vector<int> strides;
  
    try{
     ia = new Zoltan2::BasicVectorAdapter<userTypes_t>(numLocalIds, myIds,
       v_values, 1);
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 1", fail);
  
    fail = checkBasicVector(ia, numLocalIds, numLocalIds*nprocs,
      myIds, 1, valuePtrs, NULL, 0, NULL, NULL);
  
    delete ia;

    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check vector 1", fail);
  }

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
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 2", fail);
  
    fail = checkBasicVector(ia, numLocalIds, numLocalIds*nprocs,
      myIds, 1, valuePtrs, NULL, 1, weightPtrs, weightStrides);
  
    delete ia;

    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check vector 2", fail);
  }

  {
    // A Zoltan2::BasicVectorAdapter object with a multivector and no weights

    std::vector<const zscalar_t *> weightValues, values;
    std::vector<int> wstrides, vstrides;

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
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 3", fail);
  
    fail = checkBasicVector(ia, numLocalIds, numLocalIds*nprocs,
      myIds, mvdim, valuePtrs, valueStrides, 0, NULL, NULL);
  
    delete ia;

    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check vector 3", fail);
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
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 4", fail);
  
    fail = checkBasicVector(ia, numLocalIds, numLocalIds*nprocs,
      myIds, mvdim, valuePtrs, valueStrides,
      wdim, weightPtrs, weightStrides);
  
    delete ia;

    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check vector 4", fail);
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

