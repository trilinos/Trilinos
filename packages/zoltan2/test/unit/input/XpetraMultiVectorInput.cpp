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
// @HEADER
// ***********************************************************************
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
// ***********************************************************************

/*! \file XpetraMultiVectorInput.cpp
 *  \brief Test of Zoltan2::XpetraMultiVectorAdapter
 *  \todo test with weights
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

typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tvector_t;
typedef Xpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> xvector_t;

/////////////////////////////////////////////////////////////////////////
template <typename User>
int verifyInputAdapter(
  Zoltan2::XpetraMultiVectorAdapter<User> &ia, tvector_t &vector, int nvec,
    int wdim, zscalar_t **weights, int *strides)
{
  RCP<const Comm<int> > comm = vector.getMap()->getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getNumEntriesPerID() !=nvec) 
    fail = 42;

  if (!fail && ia.getNumWeightsPerID() !=wdim) 
    fail = 41;

  size_t length = vector.getLocalLength();

  if (!fail && ia.getLocalNumIDs() != length)
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

    for (int v=0; v < nvec; v++){
      auto vecdata = vector.getData(v);

      ia.getEntriesView(vals, stride, v);

      if (!fail && stride != 1)
        fail = 10;

      // check the values returned
      for (size_t i = 0; i < length; i++)
        if (vals[i*stride] != vecdata[i]) {
          fail = 104;
        }
    }

    gfail = globalFail(*comm, fail);
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

/////////////////////////////////////////////////////////////////////////

template <typename User>
int verifyGenerateFiles(
  Zoltan2::VectorAdapter<User> &ia, 
  const char *fileprefixInput,
  const Teuchos::Comm<int> &comm
)
{
  int fail = 0, gfail=0;

  const char *fileprefixGen = "unitTestOutput";
  ia.generateFiles(fileprefixGen, comm);

  // Only rank zero needs to check the resulting files
  if (comm.getRank() == 0) {

    size_t nIDsGen, nIDsInput;
    size_t nEdgeGen, nEdgeInput;
    char codeGen[4], codeInput[4];

    std::ifstream fpGen, fpInput;
    std::string graphFilenameGen = fileprefixGen;
    graphFilenameGen = graphFilenameGen + ".graph";
    std::string graphFilenameInput = fileprefixInput;
    graphFilenameInput = graphFilenameInput + ".graph";

    // Read header info from generated file
    fpGen.open(graphFilenameGen.c_str(), std::ios::in);
    std::string lineGen;
    std::getline(fpGen, lineGen);
    std::istringstream issGen(lineGen);
    issGen >> nIDsGen >> nEdgeGen >> codeGen;

    // Read header info from input file
    fpInput.open(graphFilenameInput.c_str(), std::ios::in);
    std::string lineInput;
    std::getline(fpInput, lineInput);
    while (lineInput[0]=='#') std::getline(fpInput, lineInput); // skip comments
    std::istringstream issInput(lineGen);
    issInput >> nIDsInput >> nEdgeInput >> codeInput;

    // input file and generated file should have same number of IDs
    if (nIDsGen != nIDsInput) {
      std::cout << "GenerateFiles:  nIDsGen " << nIDsGen
                << " != nIDsInput " << nIDsInput << std::endl;
      fail = 222;
    }

    // Vector adapters don't have edges
    if (!fail && nEdgeGen != 0) {
      std::cout << "GenerateFiles:  nEdgeGen " << nEdgeGen << " != 0" << std::endl;
      fail = 222;
    }

    // Check the weights, if any
    if (!fail && !strcmp(codeGen, "010")) {
      // TODO
      // If input file has weights, compare weights
      // Otherwise, just make sure there are enough weights in file
    }

    fpGen.close();
    fpInput.close();
    
    // check coordinate files
    if (!fail) {
      std::string coordsFilenameGen = fileprefixGen;
      coordsFilenameGen = coordsFilenameGen + ".coords";
      std::string coordsFilenameInput = fileprefixInput;
      coordsFilenameInput = coordsFilenameInput + ".coords";

      fpGen.open(coordsFilenameGen.c_str(), std::ios::in);
      fpInput.open(coordsFilenameInput.c_str(), std::ios::in);

      size_t cnt = 0;
      for (; std::getline(fpGen,lineGen) && std::getline(fpInput,lineInput);) { 

        cnt++;

        // Check each token
        issGen.str(lineGen);
        issInput.str(lineInput);

        while (issGen && issInput) {
          double xGen, xInput;
          issGen >> xGen;
          issInput >> xInput;

          if (xGen != xInput) {
            std::cout << "Coordinates " << xGen << " != " << xInput 
                      << std::endl;
            fail = 333;
          }
        }

        // Check same number of tokens:  are there any left in either line?
        if (issGen || issInput) {
          std::cout << "Dimension of generated file != dimension of input file"
                    << std::endl;
          fail = 334;
        }
      }

      // Did we have the correct number of coordinates?
      if (!fail && cnt != nIDsGen) {
        std::cout << "Number of coordinates read " << cnt 
                  << " != number of IDs " << nIDsGen << std::endl;
        fail = 444;
      }

      fpGen.close();
      fpInput.close();
    }
  }

  gfail = globalFail(comm, fail);
  return gfail;
}

/////////////////////////////////////////////////////////////////////////

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
  const char *inputFilePrefix = "simple";
  params.set("input file", inputFilePrefix);
  params.set("file type", "Chaco");

  try{
    uinput = rcp(new UserInputForTests(params, comm));
  }
  catch(std::exception &e){
    aok = false;
    std::cout << e.what() << std::endl;
  }
  TEST_FAIL_AND_EXIT(*comm, aok, "input ", 1);

  RCP<tvector_t> inputMVector;     // original vector (for checking)
  RCP<tvector_t> migratedMVector;   // migrated vector

  int nVec = 2;

  inputMVector = uinput->getUICoordinates();
  size_t vlen = inputMVector->getLocalLength();

  // To test migration in the input adapter we need a Solution object.

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));

  int nWeights = 1;

  typedef Zoltan2::XpetraMultiVectorAdapter<tvector_t> ia_t;
  typedef Zoltan2::PartitioningSolution<ia_t> soln_t;
  typedef ia_t::part_t part_t;

  part_t *p = new part_t [vlen];
  memset(p, 0, sizeof(part_t) * vlen);
  ArrayRCP<part_t> solnParts(p, 0, vlen, true);

  soln_t solution(env, comm, nWeights);
  solution.setParts(solnParts);

  std::vector<const zscalar_t *> emptyWeights;
  std::vector<int> emptyStrides;

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::MultiVector, no weights
  if (!gfail){ 
    if (rank==0)
      std::cout << "Constructed with Tpetra::MultiVector" << std::endl;
    
    RCP<const tvector_t> ctV = rcp_const_cast<const tvector_t>(inputMVector);
    RCP<Zoltan2::XpetraMultiVectorAdapter<tvector_t> > tVInput;
  
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
  
    fail = verifyInputAdapter<tvector_t>(*tVInput, *inputMVector, nVec, 0, NULL, NULL);
    fail = verifyGenerateFiles(*tVInput, inputFilePrefix, *comm);
  
    gfail = globalFail(*comm, fail);
  
    if (!gfail){
      tvector_t *vMigrate = NULL;
      try{
        tVInput->applyPartitioningSolution(*inputMVector, vMigrate, solution);
        migratedMVector = rcp(vMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(*comm, fail);
  
      if (!gfail){
        RCP<const tvector_t> cnewV = rcp_const_cast<const tvector_t>(migratedMVector);
        RCP<Zoltan2::XpetraMultiVectorAdapter<tvector_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraMultiVectorAdapter<tvector_t>(
            cnewV, emptyWeights, emptyStrides));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter 2 ", 1);
  
        if (rank==0){
          std::cout << "Constructed with ";
          std::cout << "Tpetra::MultiVector migrated to proc 0" << std::endl;
        }
        fail = verifyInputAdapter<tvector_t>(*newInput, *migratedMVector, nVec, 0, NULL, NULL);
        if (fail) fail += 100;
        gfail = globalFail(*comm, fail);
      }
    }
    if (gfail){
      printFailureCode(*comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::MultiVector
  if (!gfail){ 
    if (rank==0)
      std::cout << "Constructed with Xpetra::MultiVector" << std::endl;

    RCP<xvector_t> xV = Zoltan2::XpetraTraits<tvector_t>::convertToXpetra(inputMVector);
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
  
    fail = verifyInputAdapter<xvector_t>(*xVInput, *inputMVector, nVec, 0, NULL, NULL);
  
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
            rcp(new Zoltan2::XpetraMultiVectorAdapter<xvector_t>(
              cnewV, emptyWeights, emptyStrides));
        }
        catch (std::exception &e){
          aok = false;
          std::cout << e.what() << std::endl;
        }
        TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter 4 ", 1);
  
        if (rank==0){
          std::cout << "Constructed with ";
          std::cout << "Xpetra::MultiVector migrated to proc 0" << std::endl;
        }
        fail = verifyInputAdapter<xvector_t>(*newInput, *migratedMVector, nVec, 0, NULL, NULL);
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
  // User object is Epetra_MultiVector
  typedef Epetra_MultiVector evector_t;
  if (!gfail){ 
    if (rank==0)
      std::cout << "Constructed with Epetra_MultiVector" << std::endl;

    RCP<evector_t> eV = 
        rcp(new Epetra_MultiVector(uinput->getUIEpetraCrsGraph()->RowMap(),
                                   nVec));
    for (int v = 0; v < nVec; v++) {
      auto inV = inputMVector->getData(v);
      for (int i = 0; i < eV->MyLength(); i++)
        eV->ReplaceMyValue(i, v, inV[i]);
    }
        
    RCP<xvector_t> xV = Zoltan2::XpetraTraits<evector_t>::convertToXpetra(eV);
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
      fail = verifyInputAdapter<evector_t>(*eVInput, *inputMVector, nVec, 0, NULL, NULL);
  
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
            std::cout << "Epetra_MultiVector migrated to proc 0" << std::endl;
          }
          fail = verifyInputAdapter<evector_t>(*newInput, *migratedMVector, nVec, 0, NULL, NULL);
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
