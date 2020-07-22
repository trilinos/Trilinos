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
// Testing of IdentifierModel
//
//   TODO  test with BasicIdentifierAdapter with weights

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <set>
#include <bitset>

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include <Tpetra_CrsMatrix.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;

void testIdentifierModel(std::string fname, zgno_t xdim, zgno_t ydim, zgno_t zdim,
  const RCP<const Comm<int> > &comm)
{
  int rank = comm->getRank();
  int fail = 0, gfail = 0;

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags = 0;

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));

  //////////////////////////////////////////////////////////////
  // Use an Tpetra::CrsMatrix for the user data.
  //////////////////////////////////////////////////////////////
  typedef Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t> tcrsMatrix_t;
  
  UserInputForTests *uinput;
  if (fname.size() > 0)
    uinput = new UserInputForTests(testDataFilePath, fname, comm, true);
  else
    uinput = new UserInputForTests(xdim,ydim,zdim,string(""),comm, true, true);

  RCP<tcrsMatrix_t > M = uinput->getUITpetraCrsMatrix();
  zlno_t nLocalIds = M->getNodeNumRows();
  zgno_t nGlobalIds =  M->getGlobalNumRows();

  ArrayView<const zgno_t> idList = M->getRowMap()->getNodeElementList();
  std::set<zgno_t> idSet(idList.begin(), idList.end());

  //////////////////////////////////////////////////////////////
  // Create an IdentifierModel with this input
  //////////////////////////////////////////////////////////////

  typedef Zoltan2::XpetraCrsMatrixAdapter<tcrsMatrix_t> adapter_t;
  typedef Zoltan2::MatrixAdapter<tcrsMatrix_t> base_adapter_t;
  typedef Zoltan2::StridedData<zlno_t, zscalar_t> input_t;

  RCP<const adapter_t> ia = Teuchos::rcp(new adapter_t(M));
  
  Zoltan2::IdentifierModel<base_adapter_t> *model = NULL;
  RCP<const base_adapter_t> base_ia = 
                            Teuchos::rcp_dynamic_cast<const base_adapter_t>(ia);

  try{
    model = new Zoltan2::IdentifierModel<base_adapter_t>(
      base_ia, env, comm, modelFlags);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1;
  }

  gfail = globalFail(*comm, fail);

  if (gfail)
    printFailureCode(*comm, fail);
  
  // Test the IdentifierModel interface

  if (model->getLocalNumIdentifiers() != size_t(nLocalIds)) {
    std::cerr << rank << ") getLocalNumIdentifiers "
              << model->getLocalNumIdentifiers() << " "
              << nLocalIds << std::endl;
    fail = 2;
  }

  if (!fail && model->getGlobalNumIdentifiers() != size_t(nGlobalIds)) {
    std::cerr << rank << ") getGlobalNumIdentifiers "
              << model->getGlobalNumIdentifiers() << " "
              << nGlobalIds << std::endl;
    fail = 3;
  }

  gfail = globalFail(*comm, fail);

  if (gfail)
    printFailureCode(*comm, fail);
  
  ArrayView<const zgno_t> gids;
  ArrayView<input_t> wgts;
  
  model->getIdentifierList(gids, wgts);

  if (!fail && gids.size() != nLocalIds) {
    std::cerr << rank << ") getIdentifierList IDs "
              << gids.size() << " "
              << nLocalIds << std::endl;
    fail = 5;
  }

  if (!fail && wgts.size() != 0) {
    std::cerr << rank << ") getIdentifierList Weights "
              << wgts.size() << " "
              << 0 << std::endl;
    fail = 6;
  }

  for (zlno_t i=0; !fail && i < nLocalIds; i++){
    std::set<zgno_t>::iterator next = idSet.find(gids[i]);
    if (next == idSet.end()) {
      std::cerr << rank << ") getIdentifierList gid not found "
              << gids[i] << std::endl;
      fail = 7;
    }
  }

  gfail = globalFail(*comm, fail);

  if (gfail)
    printFailureCode(*comm, fail);

  delete model;
  delete uinput;
}

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();

  string fname("simple");


  if (rank == 0){
    std::cout << fname;
  }
  testIdentifierModel(fname, 0,0,0,comm);

  if (rank==0) std::cout << "PASS" << std::endl;

  return 0;
}
