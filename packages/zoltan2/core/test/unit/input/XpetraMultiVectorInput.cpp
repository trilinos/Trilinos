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
  Zoltan2::XpetraMultiVectorAdapter<User> &ia, tvector_t &vector, int nvec)
{
  RCP<const Comm<int> > comm = vector.getMap()->getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getNumEntriesPerID() !=nvec) 
    fail = 42;

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


  return gfail;
}


template <typename User>
int verifyInputAdapter(
  Zoltan2::XpetraMultiVectorAdapter<User> &ia, tvector_t &vector, int nvec,
    std::vector<const zscalar_t *> &weights, std::vector<int> &strides)
{
  // Check the input adapter 
  int gfail = verifyInputAdapter(ia, vector, nvec);
  int fail = gfail;

  RCP<const Comm<int> > comm = vector.getMap()->getComm();
  int wdim = weights.size();

  // Check the input adapter weights
  if (!gfail && (ia.getNumWeightsPerID() != wdim)) fail = 41;

  gfail = globalFail(*comm, fail);

  if (!gfail && wdim) {
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

void readChacoGraphHeaderInfo(
  std::ifstream &fp, 
  size_t &nIDs, 
  size_t &nEdges, 
  char *code, 
  int &nWgts)
{
  // Read the header info from a Chaco .graph file
  std::string line;
  std::getline(fp, line);
  while (line[0]=='#') std::getline(fp, line); // skip comments
  std::istringstream issHeader(line);
  issHeader >> nIDs >> nEdges >> code;
  if (!strcmp(code, "010") || !strcmp(code, "011")) {
    if (!(issHeader >> nWgts)) nWgts = 1;
  }
}

template <typename User>
int verifyGenerateFiles(
  Zoltan2::VectorAdapter<User> &ia, 
  const char *fileprefixInp,
  const Teuchos::Comm<int> &comm
)
{
  // Compares files generated by generateFiles (in fileprefixGen.*)
  // to input files (in fileprefixInp.*).

  int fail = 0, gfail=0;

  std::string tmp(fileprefixInp);
  tmp = tmp + "_generated";
  const char *fileprefixGen = tmp.c_str();

  ia.generateFiles(fileprefixGen, comm);

  // Only rank zero needs to check the resulting files
  if (comm.getRank() == 0) {

    size_t nIDsGen, nIDsInp;
    size_t nEdgesGen, nEdgesInp;
    char codeGen[4], codeInp[4];
    int nWgtsGen = 0, nWgtsInp = 0;
    std::string lineGen, lineInp;

    std::ifstream fpGen, fpInp;
    std::string graphFilenameGen = fileprefixGen;
    graphFilenameGen = graphFilenameGen + ".graph";
    std::string graphFilenameInp = fileprefixInp;
    graphFilenameInp = graphFilenameInp + ".graph";

    // Read header info from generated file
    fpGen.open(graphFilenameGen.c_str(), std::ios::in);
    readChacoGraphHeaderInfo(fpGen, nIDsGen, nEdgesGen, codeGen, nWgtsGen);

    // Read header info from input file
    fpInp.open(graphFilenameInp.c_str(), std::ios::in);
    readChacoGraphHeaderInfo(fpInp, nIDsInp, nEdgesInp, codeInp, nWgtsInp);

    // input file and generated file should have same number of IDs
    if (nIDsGen != nIDsInp) {
      std::cout << "GenerateFiles:  nIDsGen " << nIDsGen
                << " != nIDsInp " << nIDsInp << std::endl;
      fail = 2222;
    }

    // Vector adapters don't have edges
    if (!fail && nEdgesGen != 0) {
      std::cout << "GenerateFiles:  nEdgesGen " << nEdgesGen << " != 0" 
                << std::endl;
      fail = 2223;
    }

    // Check the weights, if any
    if (!fail && nWgtsGen) {
      if (nWgtsInp) {
        // input file has weights; compare weights
        size_t cntWgtLines;
        for (cntWgtLines = 0; cntWgtLines < nIDsGen &&
                              std::getline(fpGen, lineGen) &&
                              std::getline(fpInp, lineInp); cntWgtLines++) {

          std::istringstream issGen(lineGen);
          std::istringstream issInp(lineInp);

          int nw = 0;
          double wgtGen, wgtInp;

          while (issGen >> wgtGen) {

            if (nw < nWgtsInp) {
              issInp >> wgtInp;
              if (wgtGen != wgtInp) fail = 2224;  // weights don't match
            }
            nw++;
          }

          if (nw != nWgtsGen) fail = 2225;  // Too many or few weights on line
        }
        if (cntWgtLines != nIDsGen) fail = 2226;  // not enough input lines
      }
      else {
        // input file does not have weights; 
        // just make sure generated file has enough weights
        size_t cntWgtLines = 0;
        for (cntWgtLines = 0; cntWgtLines < nIDsGen && 
                              std::getline(fpGen, lineGen); cntWgtLines++) {
          std::istringstream issGen(lineGen);
          int nw = 0;
          double wgtGen;
          while (issGen) { issGen >> wgtGen; nw++; }
          if (nw != nWgtsGen) fail = 2227;  // Too many or few weights on line
        }
        if (cntWgtLines != nIDsGen) fail = 2228;  // not enough input lines
      }       
    }

    fpGen.close();
    fpInp.close();
    
    // check coordinate files
    if (!fail) {
      std::string coordsFilenameGen = fileprefixGen;
      coordsFilenameGen = coordsFilenameGen + ".coords";
      std::string coordsFilenameInp = fileprefixInp;
      coordsFilenameInp = coordsFilenameInp + ".coords";

      fpGen.open(coordsFilenameGen.c_str(), std::ios::in);
      fpInp.open(coordsFilenameInp.c_str(), std::ios::in);

      size_t cnt;
      for (cnt = 0; std::getline(fpGen,lineGen) &&
                    std::getline(fpInp,lineInp); cnt++) { 

        // Check each token
        std::istringstream issGen(lineGen);
        std::istringstream issInp(lineInp);

        double xGen, xInp;
        issGen >> xGen;
        issInp >> xInp;
        while (issGen && issInp) {

          if (xGen != xInp) {
            std::cout << "Coordinates " << xGen << " != " << xInp 
                      << std::endl;
            fail = 333;
          }
          issGen >> xGen;
          issInp >> xInp;
        }

        // Check same number of tokens:  are there any left in either line?
        if (issGen || issInp) {
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
      fpInp.close();
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

  // Read run-time options.
  std::string inputFilePrefix("simple");
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("file", &inputFilePrefix,
                 "chaco file (prefix) to be used in test");
  cmdp.parse(narg, arg);

  // Create object that can give us test Tpetra, Xpetra
  // and Epetra vectors for testing.

  RCP<UserInputForTests> uinput;
  Teuchos::ParameterList params;
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

  // Vertex weights
  std::vector<const zscalar_t *> IDWeights;
  std::vector<int> IDWeightsStrides;

  int nWeights = 0;
  if (uinput->hasUIWeights()) {

    auto uiweights = uinput->getUIWeights();

    nWeights = uiweights->getNumVectors();
    IDWeights.resize(nWeights);
    IDWeightsStrides.resize(nWeights);

    for (int w = 0; w < nWeights; w++) {
      IDWeights[w] = uiweights->getData(w).getRawPtr();
      IDWeightsStrides[w] = 1;
    }
  }

  // To test migration in the input adapter we need a Solution object.

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));

  typedef Zoltan2::XpetraMultiVectorAdapter<tvector_t> ia_t;
  typedef Zoltan2::PartitioningSolution<ia_t> soln_t;
  typedef ia_t::part_t part_t;

  part_t *p = new part_t [vlen];
  memset(p, 0, sizeof(part_t) * vlen);
  ArrayRCP<part_t> solnParts(p, 0, vlen, true);

  soln_t solution(env, comm, nWeights);
  solution.setParts(solnParts);


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
          IDWeights, IDWeightsStrides));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter ", 1);
  
    fail = verifyInputAdapter<tvector_t>(*tVInput, *inputMVector, nVec,
                                         IDWeights, IDWeightsStrides);
    fail = verifyGenerateFiles(*tVInput, inputFilePrefix.c_str(), *comm);
  
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
        RCP<const tvector_t> cnewV = 
            rcp_const_cast<const tvector_t>(migratedMVector);
        RCP<Zoltan2::XpetraMultiVectorAdapter<tvector_t> > newInput;
        try{
          newInput = 
             rcp(new Zoltan2::XpetraMultiVectorAdapter<tvector_t>(cnewV));
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
        // We didn't apply partitioning soln to weights,
        //  so don't check them when verifying adapter
        fail = verifyInputAdapter<tvector_t>(*newInput, *migratedMVector, nVec);
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

    RCP<xvector_t> xV =
        Zoltan2::XpetraTraits<tvector_t>::convertToXpetra(inputMVector);
    RCP<const xvector_t> cxV = rcp_const_cast<const xvector_t>(xV);
    RCP<Zoltan2::XpetraMultiVectorAdapter<xvector_t> > xVInput;
  
    try {
      xVInput = 
        rcp(new Zoltan2::XpetraMultiVectorAdapter<xvector_t>(cxV, 
                                                             IDWeights,
                                                             IDWeightsStrides));
    }
    catch (std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "XpetraMultiVectorAdapter 3 ", 1);
  
    fail = verifyInputAdapter<xvector_t>(*xVInput, *inputMVector, nVec, 
                                         IDWeights, IDWeightsStrides);
  
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
            rcp(new Zoltan2::XpetraMultiVectorAdapter<xvector_t>(cnewV));
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
        fail = verifyInputAdapter<xvector_t>(*newInput, *migratedMVector, nVec);
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
        
    RCP<const evector_t> ceV = rcp_const_cast<const evector_t>(eV);
    RCP<Zoltan2::XpetraMultiVectorAdapter<evector_t> > eVInput;
  
    bool goodAdapter = true;
    try {
      eVInput = 
        rcp(new Zoltan2::XpetraMultiVectorAdapter<evector_t>(ceV,
              IDWeights, IDWeightsStrides));
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
      fail = verifyInputAdapter<evector_t>(*eVInput, *inputMVector, nVec,
                                           IDWeights, IDWeightsStrides);
  
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
              rcp(new Zoltan2::XpetraMultiVectorAdapter<evector_t>(cnewV));
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
          fail = verifyInputAdapter<evector_t>(*newInput, *migratedMVector, 
                                               nVec);
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
