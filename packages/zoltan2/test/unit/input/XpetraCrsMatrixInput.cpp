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
// Basic testing of Zoltan2::XpetraCrsMatrixAdapter 

/*! \file XpetraCrsMatrixInput.cpp
 *  \brief Test of Zoltan2::XpetraCrsMatrixAdapter class.
 *  \todo test with geometric row coordinates.
 */

#include <string>

#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
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
typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> tmatrix_t;
typedef Xpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> xmatrix_t;
typedef Epetra_CrsMatrix ematrix_t;

void printMatrix(RCP<const Comm<int> > &comm, zlno_t nrows,
    const zgno_t *rowIds, const zlno_t *offsets, const zgno_t *colIds)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (zlno_t i=0; i < nrows; i++){
        std::cout << " row " << rowIds[i] << ": ";
        for (zlno_t j=offsets[i]; j < offsets[i+1]; j++){
          std::cout << colIds[j] << " ";
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
  Zoltan2::XpetraCrsMatrixAdapter<User> &ia, tmatrix_t &M)
{
  RCP<const Comm<int> > comm = M.getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getLocalNumRows() != M.getNodeNumRows())
    fail = 4;

  if (M.getNodeNumRows()){
    if (!fail && ia.getLocalNumColumns() != M.getNodeNumCols())
      fail = 6;
  }

  gfail = globalFail(comm, fail);

  const zgno_t *rowIds=NULL, *colIds=NULL;
  const zlno_t *offsets=NULL;
  size_t nrows=0;

  if (!gfail){

    nrows = ia.getLocalNumRows();
    ia.getRowIDsView(rowIds);
    ia.getCRSView(offsets, colIds);

    if (nrows != M.getNodeNumRows())
      fail = 8;

    gfail = globalFail(comm, fail);

    if (gfail == 0){
      printMatrix(comm, nrows, rowIds, offsets, colIds);
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
  int fail = 0, gfail=0;

  // Create object that can give us test Tpetra, Xpetra
  // and Epetra matrices for testing.

  RCP<uinput_t> uinput;

  try{
    uinput = 
      rcp(new uinput_t(
        testDataFilePath,std::string("simple"), comm, true));
  }
  catch(std::exception &e){
    TEST_FAIL_AND_EXIT(*comm, 0, string("input ")+e.what(), 1);
  }

  RCP<tmatrix_t> tM;     // original matrix (for checking)
  RCP<tmatrix_t> newM;   // migrated matrix

  tM = uinput->getUITpetraCrsMatrix();
  size_t nrows = tM->getNodeNumRows();
  Teuchos::ArrayView<const zgno_t> rowGids = 
    tM->getRowMap()->getNodeElementList();

  // To test migration in the input adapter we need a Solution
  // object.  The Solution needs an IdentifierMap.

  typedef Zoltan2::IdentifierMap<tmatrix_t> idmap_t;

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  ArrayRCP<const zgno_t> gidArray = arcpFromArrayView(rowGids);
  RCP<const idmap_t> idMap = rcp(new idmap_t(env, comm, gidArray));

  int nWeights = 1;


  typedef Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t> adapter_t;
  typedef Zoltan2::PartitioningSolution<adapter_t> soln_t;
  typedef adapter_t::part_t part_t;

  part_t *p = new part_t [nrows];
  memset(p, 0, sizeof(part_t) * nrows);
  ArrayRCP<part_t> solnParts(p, 0, nrows, true);

  soln_t solution(env, comm, idMap, nWeights);
  solution.setParts(gidArray, solnParts);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::CrsMatrix
  if (!gfail){ 
    RCP<const tmatrix_t> ctM = rcp_const_cast<const tmatrix_t>(tM);
    RCP<Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t> > tMInput;
  
    try {
      tMInput = 
        rcp(new Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t>(ctM));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraCrsMatrixAdapter ")+e.what(), 1);
    }
  
    if (rank==0)
      std::cout << "Input adapter for Tpetra::CrsMatrix" << std::endl;
    
    fail = verifyInputAdapter<tmatrix_t>(*tMInput, *tM);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      tmatrix_t *mMigrate = NULL;
      try{
        tMInput->applyPartitioningSolution(*tM, mMigrate, solution);
        newM = rcp(mMigrate);
      }
      catch (std::exception &e){
        fail = 11;
        cout << "Error caught:  " << e.what() << endl;
      }

      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const tmatrix_t> cnewM = rcp_const_cast<const tmatrix_t>(newM);
        RCP<Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t>(cnewM));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraCrsMatrixAdapter 2 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Tpetra::CrsMatrix migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<tmatrix_t>(*newInput, *newM);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::CrsMatrix
  if (!gfail){ 
    RCP<xmatrix_t> xM = uinput->getUIXpetraCrsMatrix();
    RCP<const xmatrix_t> cxM = rcp_const_cast<const xmatrix_t>(xM);
    RCP<Zoltan2::XpetraCrsMatrixAdapter<xmatrix_t> > xMInput;
  
    try {
      xMInput = 
        rcp(new Zoltan2::XpetraCrsMatrixAdapter<xmatrix_t>(cxM));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraCrsMatrixAdapter 3 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << "Input adapter for Xpetra::CrsMatrix" << std::endl;
    }
    fail = verifyInputAdapter<xmatrix_t>(*xMInput, *tM);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      xmatrix_t *mMigrate =NULL;
      try{
        xMInput->applyPartitioningSolution(*xM, mMigrate, solution);
      }
      catch (std::exception &e){
        cout << "Error caught:  " << e.what() << endl;
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const xmatrix_t> cnewM(mMigrate);
        RCP<Zoltan2::XpetraCrsMatrixAdapter<xmatrix_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraCrsMatrixAdapter<xmatrix_t>(cnewM));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraCrsMatrixAdapter 4 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Xpetra::CrsMatrix migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<xmatrix_t>(*newInput, *newM);
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
  // User object is Epetra_CrsMatrix
  if (!gfail){ 
    RCP<ematrix_t> eM = uinput->getUIEpetraCrsMatrix();
    RCP<const ematrix_t> ceM = rcp_const_cast<const ematrix_t>(eM);
    RCP<Zoltan2::XpetraCrsMatrixAdapter<ematrix_t> > eMInput;
  
    try {
      eMInput = 
        rcp(new Zoltan2::XpetraCrsMatrixAdapter<ematrix_t>(ceM));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraCrsMatrixAdapter 5 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << "Input adapter for Epetra_CrsMatrix" << std::endl;
    }
    fail = verifyInputAdapter<ematrix_t>(*eMInput, *tM);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      ematrix_t *mMigrate =NULL;
      try{
        eMInput->applyPartitioningSolution(*eM, mMigrate, solution);
      }
      catch (std::exception &e){
        cout << "Error caught:  " << e.what() << endl;
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const ematrix_t> cnewM(mMigrate, true);
        RCP<Zoltan2::XpetraCrsMatrixAdapter<ematrix_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraCrsMatrixAdapter<ematrix_t>(cnewM));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraCrsMatrixAdapter 6 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Epetra_CrsMatrix migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<ematrix_t>(*newInput, *newM);
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
