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
// Basic testing of Zoltan2::TpetraRowMatrixAdapter 

/*! \file TpetraRowMatrixInput.cpp
 *  \brief Test of Zoltan2::TpetraRowMatrixAdapter class.
 *  \todo test with geometric row coordinates.
 */

#include <string>

#include <Zoltan2_TpetraRowMatrixAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::Comm;

typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> ztcrsmatrix_t;
typedef Tpetra::RowMatrix<zscalar_t, zlno_t, zgno_t, znode_t> ztrowmatrix_t;

//////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////

template <typename User>
int verifyInputAdapter(
  Zoltan2::TpetraRowMatrixAdapter<User> &ia, ztrowmatrix_t &M)
{
  typedef typename Zoltan2::InputTraits<User>::offset_t offset_t;

  RCP<const Comm<int> > comm = M.getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getLocalNumRows() != M.getNodeNumRows())
    fail = 4;

  if (M.getNodeNumRows()){
    if (!fail && ia.getLocalNumColumns() != M.getNodeNumCols())
      fail = 6;
  }

  gfail = globalFail(*comm, fail);

  const zgno_t *rowIds=NULL, *colIds=NULL;
  const offset_t *offsets=NULL;
  size_t nrows=0;

  if (!gfail){

    nrows = ia.getLocalNumRows();
    ia.getRowIDsView(rowIds);
    ia.getCRSView(offsets, colIds);

    if (nrows != M.getNodeNumRows())
      fail = 8;

    gfail = globalFail(*comm, fail);

    if (gfail == 0){
      printMatrix<offset_t>(comm, nrows, rowIds, offsets, colIds);
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

  // Create object that can give us Tpetra matrices for testing.

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

  // Input matrix and row matrix cast from it.
  RCP<ztcrsmatrix_t> tM = uinput->getUITpetraCrsMatrix();
  RCP<ztrowmatrix_t> trM = rcp_dynamic_cast<ztrowmatrix_t>(tM);

  RCP<ztrowmatrix_t> newM;   // migrated matrix

  size_t nrows = trM->getNodeNumRows();

  // To test migration in the input adapter we need a Solution object. 

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));

  int nWeights = 1;

  typedef Zoltan2::TpetraRowMatrixAdapter<ztrowmatrix_t> adapter_t;
  typedef Zoltan2::PartitioningSolution<adapter_t> soln_t;
  typedef adapter_t::part_t part_t;

  part_t *p = new part_t [nrows];
  memset(p, 0, sizeof(part_t) * nrows);
  ArrayRCP<part_t> solnParts(p, 0, nrows, true);

  soln_t solution(env, comm, nWeights);
  solution.setParts(solnParts);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::RowMatrix
  if (!gfail){ 
    if (rank==0)
      std::cout << "Input adapter for Tpetra::RowMatrix" << std::endl;
    
    RCP<const ztrowmatrix_t> ctrM = rcp_const_cast<const ztrowmatrix_t>(
                                   rcp_dynamic_cast<ztrowmatrix_t>(tM));
    RCP<adapter_t> trMInput;
  
    try {
      trMInput = rcp(new adapter_t(ctrM));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "TpetraRowMatrixAdapter ", 1);
  
    fail = verifyInputAdapter<ztrowmatrix_t>(*trMInput, *trM);
  
    gfail = globalFail(*comm, fail);
  
    if (!gfail){
      ztrowmatrix_t *mMigrate = NULL;
      try{
        trMInput->applyPartitioningSolution(*trM, mMigrate, solution);
        newM = rcp(mMigrate);
      }
      catch (std::exception &e){
        fail = 11;
        std::cout << "Error caught:  " << e.what() << std::endl;
      }

      gfail = globalFail(*comm, fail);
  
      if (!gfail){
        RCP<const ztrowmatrix_t> cnewM = 
                                 rcp_const_cast<const ztrowmatrix_t>(newM);
        RCP<adapter_t> newInput;
        try{
          newInput = rcp(new adapter_t(cnewM));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "TpetraRowMatrixAdapter 2 ", 1);
  
        if (rank==0){
          std::cout << 
           "Input adapter for Tpetra::RowMatrix migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<ztrowmatrix_t>(*newInput, *newM);
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
