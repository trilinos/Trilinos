// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::Comm;

typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> tmatrix_t;
typedef Xpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> xmatrix_t;

template<typename offset_t>
void printMatrix(RCP<const Comm<int> > &comm, zlno_t nrows,
    const zgno_t *rowIds, const offset_t *offsets, const zgno_t *colIds)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (zlno_t i=0; i < nrows; i++){
        std::cout << " row " << rowIds[i] << ": ";
        for (offset_t j=offsets[i]; j < offsets[i+1]; j++){
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
  typedef typename Zoltan2::InputTraits<User>::offset_t offset_t;

  RCP<const Comm<int> > comm = M.getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getLocalNumRows() != M.getLocalNumRows())
    fail = 4;

  if (M.getLocalNumRows()){
    if (!fail && ia.getLocalNumColumns() != M.getLocalNumCols())
      fail = 6;
  }

  gfail = globalFail(*comm, fail);

  const zgno_t *rowIds=NULL;
  ArrayRCP<const zgno_t> colIds;
  ArrayRCP<const offset_t> offsets;
  size_t nrows=0;

  if (!gfail){

    nrows = ia.getLocalNumRows();
    ia.getRowIDsView(rowIds);
    ia.getCRSView(offsets, colIds);

    if (nrows != M.getLocalNumRows())
      fail = 8;

    gfail = globalFail(*comm, fail);

    if (gfail == 0){
      printMatrix<offset_t>(comm, nrows, rowIds, offsets.getRawPtr(), colIds.getRawPtr());
    }
    else{
      if (!fail) fail = 10;
    }
  }
  return fail;
}

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int fail = 0, gfail=0;
  bool aok = true;

  // Create object that can give us test Tpetra, Xpetra
  // and Epetra matrices for testing.

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

  RCP<tmatrix_t> tM;     // original matrix (for checking)
  RCP<tmatrix_t> newM;   // migrated matrix

  tM = uinput->getUITpetraCrsMatrix();
  size_t nrows = tM->getLocalNumRows();

  // To test migration in the input adapter we need a Solution object. 

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));

  int nWeights = 1;

  typedef Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t> adapter_t;
  typedef Zoltan2::PartitioningSolution<adapter_t> soln_t;
  typedef adapter_t::part_t part_t;

  part_t *p = new part_t [nrows];
  memset(p, 0, sizeof(part_t) * nrows);
  ArrayRCP<part_t> solnParts(p, 0, nrows, true);

  soln_t solution(env, comm, nWeights);
  solution.setParts(solnParts);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::CrsMatrix
  if (!gfail){ 
    if (rank==0)
      std::cout << "Input adapter for Tpetra::CrsMatrix" << std::endl;
    
    RCP<const tmatrix_t> ctM = rcp_const_cast<const tmatrix_t>(tM);
    RCP<Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t> > tMInput;
  
    try {
      tMInput = 
        rcp(new Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t>(ctM));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsMatrixAdapter ", 1);
  
    fail = verifyInputAdapter<tmatrix_t>(*tMInput, *tM);
  
    gfail = globalFail(*comm, fail);
  
    if (!gfail){
      tmatrix_t *mMigrate = NULL;
      try{
        tMInput->applyPartitioningSolution(*tM, mMigrate, solution);
        newM = rcp(mMigrate);
      }
      catch (std::exception &e){
        fail = 11;
        std::cout << "Error caught:  " << e.what() << std::endl;
      }

      gfail = globalFail(*comm, fail);
  
      if (!gfail){
        RCP<const tmatrix_t> cnewM = rcp_const_cast<const tmatrix_t>(newM);
        RCP<Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraCrsMatrixAdapter<tmatrix_t>(cnewM));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsMatrixAdapter 2 ", 1);
  
        if (rank==0){
          std::cout << 
           "Input adapter for Tpetra::CrsMatrix migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<tmatrix_t>(*newInput, *newM);
        if (fail) fail += 100;
        gfail = globalFail(*comm, fail);
      }
    }
    if (gfail){
      printFailureCode(*comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::CrsMatrix
  if (!gfail){ 
    if (rank==0)
      std::cout << "Input adapter for Xpetra::CrsMatrix" << std::endl;

    RCP<xmatrix_t> xM = uinput->getUIXpetraCrsMatrix();
    RCP<const xmatrix_t> cxM = rcp_const_cast<const xmatrix_t>(xM);
    RCP<Zoltan2::XpetraCrsMatrixAdapter<xmatrix_t> > xMInput;
  
    try {
      xMInput = 
        rcp(new Zoltan2::XpetraCrsMatrixAdapter<xmatrix_t>(cxM));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsMatrixAdapter 3 ", 1);
  
    fail = verifyInputAdapter<xmatrix_t>(*xMInput, *tM);
  
    gfail = globalFail(*comm, fail);
  
    if (!gfail){
      xmatrix_t *mMigrate =NULL;
      try{
        xMInput->applyPartitioningSolution(*xM, mMigrate, solution);
      }
      catch (std::exception &e){
        std::cout << "Error caught:  " << e.what() << std::endl;
        fail = 11;
      }
  
      gfail = globalFail(*comm, fail);
  
      if (!gfail){
        RCP<const xmatrix_t> cnewM(mMigrate);
        RCP<Zoltan2::XpetraCrsMatrixAdapter<xmatrix_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraCrsMatrixAdapter<xmatrix_t>(cnewM));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "XpetraCrsMatrixAdapter 4 ", 1);
  
        if (rank==0){
          std::cout << 
           "Input adapter for Xpetra::CrsMatrix migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<xmatrix_t>(*newInput, *newM);
        if (fail) fail += 100;
        gfail = globalFail(*comm, fail);
      }
    }
    if (gfail){
      printFailureCode(*comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // DONE

  if (rank==0)
    std::cout << "PASS" << std::endl;
}
