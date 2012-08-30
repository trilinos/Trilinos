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
//   TODO  test with BasicIdentifierInput with weights

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_BasicIdentifierInput.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <set>
#include <bitset>

#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include <Tpetra_CrsMatrix.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::DefaultComm;

void testIdentifierModel(std::string fname, gno_t xdim, gno_t ydim, gno_t zdim,
  const RCP<const Comm<int> > &comm, bool consecutiveIds)
{
  int rank = comm->getRank();
  int fail = 0, gfail = 0;

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags = 0;
  if (consecutiveIds)
    modelFlags.set(Zoltan2::IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  //////////////////////////////////////////////////////////////
  // Use an Tpetra::CrsMatrix for the user data.
  //////////////////////////////////////////////////////////////
  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t> tcrsMatrix_t;
  
  UserInputForTests *input;
  if (fname.size() > 0)
    input = new UserInputForTests(testDataFilePath, fname, comm, true);
  else
    input = new UserInputForTests(xdim,ydim,zdim,string(""),comm, true);

  RCP<tcrsMatrix_t > M = input->getTpetraCrsMatrix();
  lno_t nLocalIds = M->getNodeNumRows();
  gno_t nGlobalIds =  M->getGlobalNumRows();

  ArrayView<const gno_t> idList = M->getRowMap()->getNodeElementList();
  std::set<gno_t> idSet(idList.begin(), idList.end());

  //////////////////////////////////////////////////////////////
  // Create an IdentifierModel with this input
  //////////////////////////////////////////////////////////////

  typedef Zoltan2::XpetraCrsMatrixInput<tcrsMatrix_t> adapter_t;
  typedef Zoltan2::MatrixInput<tcrsMatrix_t> base_adapter_t;
  typedef Zoltan2::StridedData<lno_t, scalar_t> input_t;

  RCP<const adapter_t> ia = Teuchos::rcp(new adapter_t(M));
  
  Zoltan2::IdentifierModel<base_adapter_t> *model = NULL;
  const base_adapter_t *base_ia = ia.get();

  try{
    model = new Zoltan2::IdentifierModel<base_adapter_t>(
      base_ia, env, comm, modelFlags);
  }
  catch (std::exception &e){
    std::cerr << rank << ") " << e.what() << std::endl;
    fail = 1;
  }

  gfail = globalFail(comm, fail);

  if (gfail)
    printFailureCode(comm, fail);
  
  // Test the IdentifierModel interface

  if (model->getLocalNumIdentifiers() != size_t(nLocalIds))
    fail = 2;

  if (!fail && model->getGlobalNumIdentifiers() != size_t(nGlobalIds))
    fail = 3;

  // For now, MatrixInput does not implement weights
  if (!fail && model->getIdentifierWeightDim() !=  0)
    fail = 4;

  gfail = globalFail(comm, fail);

  if (gfail)
    printFailureCode(comm, fail);
  
  ArrayView<const gno_t> gids;
  ArrayView<input_t> wgts;
  
  model->getIdentifierList(gids, wgts);

  if (!fail && gids.size() != nLocalIds)
    fail = 5;

  if (!fail && wgts.size() != 0)
    fail = 6;

  for (lno_t i=0; !fail && i < nLocalIds; i++){
    std::set<gno_t>::iterator next = idSet.find(gids[i]);
    if (next == idSet.end())
      fail = 7;
  }

  if (!fail && consecutiveIds){
    bool inARow = Zoltan2::IdentifierTraits<gno_t>::areConsecutive(
      gids.getRawPtr(), nLocalIds);

    if (!inARow)
      fail = 8;
  }

  gfail = globalFail(comm, fail);

  if (gfail)
    printFailureCode(comm, fail);

  delete model;
  delete input;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  string fname("simple");

  bool wishConsecutiveIds = true;

  if (rank == 0){
    std::cout << fname;
    std::cout << ", consecutive IDs not requested" << std::endl;
  }
  testIdentifierModel(fname, 0,0,0,comm, !wishConsecutiveIds);

  if (rank == 0){
    std::cout << fname;
    std::cout << ", consecutive IDs are requested" << std::endl;
  }
  testIdentifierModel(fname, 0,0,0,comm,  wishConsecutiveIds);

  if (rank==0) std::cout << "PASS" << std::endl;

  return 0;
}
