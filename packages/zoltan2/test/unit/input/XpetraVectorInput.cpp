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

/*! \file XpetraVectorInput.cpp
 *  \brief Test of Zoltan2::XpetraVectorInput class.
 *  \todo add test with weights
 */

#include <string>

#include <Zoltan2_XpetraVectorInput.hpp>
#include <Zoltan2_InputTraits.hpp>
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

typedef UserInputForTests uinput_t;
typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> tvector_t;
typedef Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> xvector_t;
typedef Epetra_Vector evector_t;

void printVector(RCP<const Comm<int> > &comm, lno_t vlen,
    const gno_t *vtxIds, const scalar_t *vals)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (lno_t i=0; i < vlen; i++){
        std::cout << " " << vtxIds[i] << ": " << vals[i] << std::endl;
      }
      std::cout.flush();
    }
    comm->barrier();
  }
  comm->barrier();
}

template <typename User>
int verifyInputAdapter(
  Zoltan2::XpetraVectorInput<User> &ia, tvector_t &vector, int wdim, 
    scalar_t **weights, int *strides)
{
  RCP<const Comm<int> > comm = vector.getMap()->getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getNumberOfVectors() !=1) 
    fail = 42;

  if (!fail && ia.getNumberOfWeights() !=wdim) 
    fail = 41;

  if (!fail && ia.getLocalLength() != vector.getLocalLength())
    fail = 4;

  gfail = globalFail(comm, fail);

  if (!gfail){
    const gno_t *vtxIds=NULL;
    const scalar_t *vals=NULL;
    int stride;

    size_t nvals = ia.getVector(vtxIds, vals, stride);

    if (nvals != vector.getLocalLength())
      fail = 8;
    if (!fail && stride != 1)
      fail = 10;

    gfail = globalFail(comm, fail);

    if (gfail == 0){
      printVector(comm, nvals, vtxIds, vals);
    }
  }

  if (!gfail && wdim){
    const scalar_t *wgt =NULL;
    int stride;

    for (int w=0; !fail && w < wdim; w++){
      size_t nvals = ia.getVectorWeights(w, wgt, stride);

      if (nvals != vector.getLocalLength())
        fail = 100;
      if (!fail && stride != strides[w])
        fail = 101;

      for (size_t v=0; !fail && v < nvals; v++){
        if (wgt[v*stride] != weights[w][v*stride])
          fail=102;
      }
    }

    gfail = globalFail(comm, fail);
  }

  return gfail;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int fail = 0, gfail=0;

  // Create object that can give us test Tpetra, Xpetra
  // and Epetra vectors for testing.

  RCP<uinput_t> uinput;

  try{
    uinput = 
      rcp(new uinput_t(testDataFilePath,std::string("simple"), comm, true));
  }
  catch(std::exception &e){
    TEST_FAIL_AND_EXIT(*comm, 0, string("input ")+e.what(), 1);
  }

  RCP<tvector_t> tV;     // original vector (for checking)
  RCP<tvector_t> newV;   // migrated vector

  tV = uinput->getTpetraVector();
  size_t vlen = tV->getLocalLength();
  Teuchos::ArrayView<const gno_t> rowGids = tV->getMap()->getNodeElementList();

  // To test migration in the input adapter we need a Solution
  // object.  The Solution needs an IdentifierMap.

  typedef Zoltan2::IdentifierMap<tvector_t> idmap_t;

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  ArrayRCP<const gno_t> gidArray = arcpFromArrayView(rowGids);
  RCP<const idmap_t> idMap = rcp(new idmap_t(env, comm, gidArray));

  int weightDim = 1;

  zoltan2_partId_t *p = new zoltan2_partId_t [vlen];
  memset(p, 0, sizeof(zoltan2_partId_t) * vlen);
  ArrayRCP<zoltan2_partId_t> solnParts(p, 0, vlen, true);

  std::vector<const scalar_t *> emptyWeights;
  std::vector<int> emptyStrides;

  typedef Zoltan2::XpetraVectorInput<tvector_t> adapter_t;
  Zoltan2::PartitioningSolution<adapter_t> solution(
    env, comm, idMap, weightDim);
  solution.setParts(gidArray, solnParts);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::Vector, no weights
  if (!gfail){ 
    RCP<const tvector_t> ctV = rcp_const_cast<const tvector_t>(tV);
    RCP<adapter_t> tVInput;
  
    try {
      tVInput = 
        rcp(new Zoltan2::XpetraVectorInput<tvector_t>(ctV, 
          emptyWeights, emptyStrides));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraVectorInput ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << tVInput->inputAdapterName() << ", constructed with ";
      std::cout  << "Tpetra::Vector" << std::endl;
    }
    
    fail = verifyInputAdapter<tvector_t>(*tVInput, *tV, 0, NULL, NULL);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      tvector_t *vMigrate = NULL;
      try{
        tVInput->applyPartitioningSolution(*tV, vMigrate, solution);
        newV = rcp(vMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const tvector_t> cnewV = rcp_const_cast<const tvector_t>(newV);
        RCP<Zoltan2::XpetraVectorInput<tvector_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraVectorInput<tvector_t>(cnewV,
            emptyWeights, emptyStrides));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraVectorInput 2 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << tVInput->inputAdapterName() << ", constructed with ";
          std::cout << "Tpetra::Vector migrated to proc 0" << std::endl;
        }
        fail = verifyInputAdapter<tvector_t>(*newInput, *newV, 0, NULL, NULL);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::Vector
  if (!gfail){ 
    RCP<xvector_t> xV = uinput->getXpetraVector();
    RCP<const xvector_t> cxV = rcp_const_cast<const xvector_t>(xV);
    RCP<Zoltan2::XpetraVectorInput<xvector_t> > xVInput;
  
    try {
      xVInput = 
        rcp(new Zoltan2::XpetraVectorInput<xvector_t>(cxV,
          emptyWeights, emptyStrides));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraVectorInput 3 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << xVInput->inputAdapterName() << ", constructed with ";
      std::cout << "Xpetra::Vector" << std::endl;
    }
    fail = verifyInputAdapter<xvector_t>(*xVInput, *tV, 0, NULL, NULL);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      xvector_t *vMigrate =NULL;
      try{
        xVInput->applyPartitioningSolution(*xV, vMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const xvector_t> cnewV(vMigrate);
        RCP<Zoltan2::XpetraVectorInput<xvector_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraVectorInput<xvector_t>(cnewV, 
              emptyWeights, emptyStrides));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraVectorInput 4 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << xVInput->inputAdapterName() << ", constructed with ";
          std::cout << "Xpetra::Vector migrated to proc 0" << std::endl;
        }
        fail = verifyInputAdapter<xvector_t>(*newInput, *newV, 0, NULL, NULL);
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
  // User object is Epetra_Vector
  if (!gfail){ 
    RCP<evector_t> eV = uinput->getEpetraVector();
    RCP<const evector_t> ceV = rcp_const_cast<const evector_t>(eV);
    RCP<Zoltan2::XpetraVectorInput<evector_t> > eVInput;
  
    try {
      eVInput = 
        rcp(new Zoltan2::XpetraVectorInput<evector_t>(ceV,
          emptyWeights, emptyStrides));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraVectorInput 5 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << eVInput->inputAdapterName() << ", constructed with ";
      std::cout << "Epetra_Vector" << std::endl;
    }
    fail = verifyInputAdapter<evector_t>(*eVInput, *tV, 0, NULL, NULL);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      evector_t *vMigrate =NULL;
      try{
        eVInput->applyPartitioningSolution(*eV, vMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const evector_t> cnewV(vMigrate, true);
        RCP<Zoltan2::XpetraVectorInput<evector_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraVectorInput<evector_t>(cnewV, 
              emptyWeights, emptyStrides));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraVectorInput 6 ")+e.what(), 1);
        }
  
        if (rank==0){
           std::cout << eVInput->inputAdapterName() << ", constructed with ";
           std::cout << "Epetra_Vector migrated to proc 0" << std::endl;
        }
        fail = verifyInputAdapter<evector_t>(*newInput, *newV, 0, NULL, NULL);
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
