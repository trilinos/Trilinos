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
 *  \brief Test of Zoltan2::XpetraMultiVectorAdapter class with vector input.
 *  \todo add test with weights
 */

#include <string>

#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::Comm;

typedef Tpetra::Vector<zscalar_t, zlno_t, zgno_t, znode_t> tvector_t;
typedef Xpetra::Vector<zscalar_t, zlno_t, zgno_t, znode_t> xvector_t;

void printVector(RCP<const Comm<int> > &comm, zlno_t vlen,
    const zgno_t *vtxIds, const zscalar_t *vals)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (zlno_t i=0; i < vlen; i++){
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
  Zoltan2::XpetraMultiVectorAdapter<User> &ia, tvector_t &vector, int wdim, 
    zscalar_t **weights, int *strides)
{
  RCP<const Comm<int> > comm = vector.getMap()->getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getNumEntriesPerID() !=1) 
    fail = 42;

  if (!fail && ia.getNumWeightsPerID() !=wdim) 
    fail = 41;

  if (!fail && ia.getLocalNumIDs() != vector.getLocalLength())
    fail = 4;

  gfail = globalFail(*comm, fail);

  if (!gfail){
    const zgno_t *vtxIds=NULL;
    const zscalar_t *vals=NULL;
    int stride;

    size_t nvals = ia.getLocalNumIDs();
    if (nvals != vector.getLocalLength())
      fail = 8;

    ia.getIDsView(vtxIds);
    ia.getEntriesView(vals, stride);

    if (!fail && stride != 1)
      fail = 10;

    gfail = globalFail(*comm, fail);

    if (gfail == 0){
      printVector(comm, nvals, vtxIds, vals);
    }
  }

  if (!gfail && wdim){
    const zscalar_t *wgt =NULL;
    int stride;

    for (int w=0; !fail && w < wdim; w++){
      ia.getWeightsView(wgt, stride, w);

      if (!fail && stride != strides[w])
        fail = 101;

      for (size_t v=0; !fail && v < vector.getLocalLength(); v++){
        if (wgt[v*stride] != weights[w][v*stride])
          fail=102;
      }
    }

    gfail = globalFail(*comm, fail);
  }

  return gfail;
}


int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int fail = 0, gfail=0;
  bool aok = true;

  // Create object that can give us test Tpetra, Xpetra
  // and Epetra vectors for testing.

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

  RCP<tvector_t> tV;     // original vector (for checking)
  RCP<tvector_t> newV;   // migrated vector

  tV = rcp(new tvector_t(uinput->getUITpetraCrsGraph()->getRowMap()));
  tV->randomize();
  size_t vlen = tV->getLocalLength();

  // To test migration in the input adapter we need a Solution object.

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));

  int nWeights = 1;

  typedef Zoltan2::XpetraMultiVectorAdapter<tvector_t> adapter_t;
  typedef adapter_t::part_t part_t;

  part_t *p = new part_t [vlen];
  memset(p, 0, sizeof(part_t) * vlen);
  ArrayRCP<part_t> solnParts(p, 0, vlen, true);

  std::vector<const zscalar_t *> emptyWeights;
  std::vector<int> emptyStrides;

  Zoltan2::PartitioningSolution<adapter_t> solution(env, comm, nWeights);
  solution.setParts(solnParts);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::Vector, no weights
  if (!gfail){ 
    if (rank==0)
      std::cout << "Constructed with Tpetra::Vector" << std::endl;
    
    RCP<const tvector_t> ctV = rcp_const_cast<const tvector_t>(tV);
    RCP<adapter_t> tVInput;
  
    try {
      tVInput = 
        rcp(new Zoltan2::XpetraMultiVectorAdapter<tvector_t>(ctV, 
          emptyWeights, emptyStrides));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter ", 1);
  
    fail = verifyInputAdapter<tvector_t>(*tVInput, *tV, 0, NULL, NULL);
  
    gfail = globalFail(*comm, fail);
  
    if (!gfail){
      tvector_t *vMigrate = NULL;
      try{
        tVInput->applyPartitioningSolution(*tV, vMigrate, solution);
        newV = rcp(vMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(*comm, fail);
  
      if (!gfail){
        RCP<const tvector_t> cnewV = rcp_const_cast<const tvector_t>(newV);
        RCP<Zoltan2::XpetraMultiVectorAdapter<tvector_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraMultiVectorAdapter<tvector_t>(cnewV,
            emptyWeights, emptyStrides));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter 2 ", 1);
  
        if (rank==0){
          std::cout << "Constructed with ";
          std::cout << "Tpetra::Vector migrated to proc 0" << std::endl;
        }
        fail = verifyInputAdapter<tvector_t>(*newInput, *newV, 0, NULL, NULL);
        if (fail) fail += 100;
        gfail = globalFail(*comm, fail);
      }
    }
    if (gfail){
      printFailureCode(*comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::Vector
  if (!gfail){ 
    if (rank==0)
      std::cout << "Constructed with Xpetra::Vector" << std::endl;

    RCP<tvector_t> xtV =
        rcp(new tvector_t(uinput->getUITpetraCrsGraph()->getRowMap()));
    xtV->randomize();
    RCP<xvector_t> xV = Zoltan2::XpetraTraits<tvector_t>::convertToXpetra(xtV);
    RCP<const xvector_t> cxV = rcp_const_cast<const xvector_t>(xV);
    RCP<Zoltan2::XpetraMultiVectorAdapter<xvector_t> > xVInput;
  
    try {
      xVInput = 
        rcp(new Zoltan2::XpetraMultiVectorAdapter<xvector_t>(cxV,
          emptyWeights, emptyStrides));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter 3 ", 1);
  
    fail = verifyInputAdapter<xvector_t>(*xVInput, *tV, 0, NULL, NULL);
  
    gfail = globalFail(*comm, fail);
  
    if (!gfail){
      xvector_t *vMigrate =NULL;
      try{
        xVInput->applyPartitioningSolution(*xV, vMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }
  
      gfail = globalFail(*comm, fail);
  
      if (!gfail){
        RCP<const xvector_t> cnewV(vMigrate);
        RCP<Zoltan2::XpetraMultiVectorAdapter<xvector_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraMultiVectorAdapter<xvector_t>(cnewV, 
              emptyWeights, emptyStrides));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter 4 ", 1);
  
        if (rank==0){
          std::cout << "Constructed with ";
          std::cout << "Xpetra::Vector migrated to proc 0" << std::endl;
        }
        fail = verifyInputAdapter<xvector_t>(*newInput, *newV, 0, NULL, NULL);
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
  // User object is Epetra_Vector
  typedef Epetra_Vector evector_t;
  if (!gfail){ 
    if (rank==0)
      std::cout << "Constructed with Epetra_Vector" << std::endl;

    RCP<evector_t> eV = 
        rcp(new Epetra_Vector(uinput->getUIEpetraCrsGraph()->RowMap()));
    eV->Random();
    RCP<const evector_t> ceV = rcp_const_cast<const evector_t>(eV);
    RCP<Zoltan2::XpetraMultiVectorAdapter<evector_t> > eVInput;
  
    bool goodAdapter = true;
    try {
      eVInput = 
        rcp(new Zoltan2::XpetraMultiVectorAdapter<evector_t>(ceV,
          emptyWeights, emptyStrides));
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
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter 5 ", 1);
  
    if (goodAdapter) {
      fail = verifyInputAdapter<evector_t>(*eVInput, *tV, 0, NULL, NULL);
  
      gfail = globalFail(*comm, fail);
  
      if (!gfail){
        evector_t *vMigrate =NULL;
        try{
          eVInput->applyPartitioningSolution(*eV, vMigrate, solution);
        }
        catch (std::exception &e){
          fail = 11;
        }
  
        gfail = globalFail(*comm, fail);
  
        if (!gfail){
          RCP<const evector_t> cnewV(vMigrate, true);
          RCP<Zoltan2::XpetraMultiVectorAdapter<evector_t> > newInput;
          try{
            newInput = 
              rcp(new Zoltan2::XpetraMultiVectorAdapter<evector_t>(cnewV, 
                emptyWeights, emptyStrides));
          }
          catch (std::exception &e){
            aok = false;
            std::cout << e.what() << std::endl;
          }
          TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter 6 ", 1);
    
          if (rank==0){
             std::cout << "Constructed with ";
             std::cout << "Epetra_Vector migrated to proc 0" << std::endl;
          }
          fail = verifyInputAdapter<evector_t>(*newInput, *newV, 0, NULL, NULL);
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
