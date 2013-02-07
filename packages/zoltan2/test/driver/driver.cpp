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

/*! \file  driver.cpp
    \brief Test driver for Zoltan2 library

    Command line argument is the data file of test definitions.
    If omitted, we use "./driver.xml".

    The goal of this test is to test algorithms.  Components like
    InputAdapter and Model classes are tested in unit tests. 

    \todo There are no defaults.  Everything about the test has to
              be given in the XML input file.
    \todo There are many other useful options that could be added
              to the test options in the input file.
    \todo Add timing and memory profiling.
    \todo write createInput()
    \todo write evaluateSuccess()
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolutionQuality.hpp>
#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_BasicIdentifierInput.hpp>
#include <Zoltan2_BasicVectorInput.hpp>
#include <Zoltan2_XpetraCrsGraphInput.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_XpetraMultiVectorInput.hpp>
#include <Zoltan2_XpetraVectorInput.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <sstream>
#include <string>
#include <iostream>
#include <vector>

using Teuchos::ParameterList;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::XMLObject;

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::exception;
using std::ostringstream;
using std::vector;

using Zoltan2::PartitioningProblem;

typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;
typedef Tpetra::CrsGraph<lno_t, gno_t, node_t> tcrsGraph_t;
typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> tVector_t;
typedef Tpetra::Map<lno_t, gno_t, node_t> tmap_t;

typedef Zoltan2::BasicCoordinateAdapter<tcrsMatrix_t> bci_t;
typedef Zoltan2::BasicIdentifierAdapter<tcrsMatrix_t> bii_t;
typedef Zoltan2::BasicVectorAdapter<tcrsMatrix_t>     bvi_t; 
typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t>  xgi_t;
typedef Zoltan2::XpetraCrsMatrixAdapter<tcrsMatrix_t> xmi_t;
typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t>  xmvi_t;
typedef Zoltan2::XpetraVectorAdapter<tVector_t>    xvi_t;

#define ERRMSG(msg) if (rank == 0){ cerr << "FAIL: " << msg << endl; }
#define EXC_ERRMSG(msg, e) \
  if (rank==0){ cerr << "FAIL: " << msg << endl << e.what() << endl;}

template <typename Adapter>
  int runAlgorithm(int rank, Adapter *ia, ParameterList &pList,
    ArrayRCP<const Zoltan2::MetricValues<scalar_t> > &quality)
{
  int fail = 0;

  PartitioningProblem<Adapter> *p = NULL;

  try{
    p = new PartitioningProblem<Adapter>(ia, &pList, MPI_COMM_WORLD);
  }
  catch (std::exception &e){
    EXC_ERRMSG("Test error: Problem build", e);
    fail = 1;
  }

  if (!fail){
    try{
      p->solve();
    }
    catch (std::exception &e){
      EXC_ERRMSG("Test error: Solve", e);
      fail = 1;
    }
  }

  if (!fail)
    quality = p->getMetrics();

  return fail;
}
        
void createInput(const XMLObject &testInfo, 
  tcrsMatrix_t *&M, tMVector_t *&objWgts, tMVector_t *&edgeWgts, 
  tMVector_t * &coords)
{
  M = NULL;
  objWgts = NULL;
  edgeWgts = NULL;
  coords = NULL;

  // UserInputForTests or GeometryGenerator can create
  // the input.
}

int evaluateSuccess(const XMLObject &criteria,
    ArrayRCP<const Zoltan2::MetricValues<scalar_t> > &quality)
{
  int fail = 0;

  return fail;
}

int main(int argc, char *argv[])
{
  ///////////////////////////////////
  // Set up parallel environment.

  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  ///////////////////////////////////
  // Read the input file of tests to run.

  string inFileName("driver.xml");
  if (argc > 1)
    inFileName = string(argv[1]);

  Teuchos::FileInputSource inFile(inFileName);
  XMLObject xmlIn;

  try{
    xmlIn = inFile.getObject();
  }
  catch (exception &e){
    EXC_ERRMSG("Driver error: reading ", e);
    return 1;
  }

  ///////////////////////////////////
  // Obtain the list of valid Zoltan2 parameters.

  ParameterList zoltan2Parameters;
  Zoltan2::createAllParameters(zoltan2Parameters);

  ///////////////////////////////////
  // Loop over tests and execute them.

  const string testsTag = xmlIn.getTag();

  if (testsTag != string("Tests")){
    ERRMSG("Driver error: xml file does not define \"Tests\"");
    return 1;
  }

  const string testSuiteName = xmlIn.getAttribute(string("name"));
  int numTests = xmlIn.numChildren();

  if (rank == 0){
    cout << "Test suite: " << testSuiteName << ", ";
    cout << numTests << " tests" << endl;
  }

  int numPass = 0;

  for (int testNum=0; testNum < numTests; testNum++){

    ////
    // Get the next test.
    ////

    const XMLObject nextTest = xmlIn.getChild(testNum);
    const string tName = xmlIn.getAttribute(string("name"));

    if (rank==0)
      cout << "Test " << tName << endl;

    int tpLoc = nextTest.findFirstChild("TestParameters");
    int zpLoc = nextTest.findFirstChild("Zoltan2Parameters");

    if (tpLoc < 0){
      ostringstream msg;
      msg << "Driver error: no TestParameters for " << tName; 
      ERRMSG(msg.str());
    }

    if (zpLoc < 0){
      ostringstream msg;
      msg << "Driver error: no Zoltan2Parameters for " << tName; 
      ERRMSG(msg.str());
      continue;
    }

    ////
    // Create the Zoltan2 parameter list.
    ////

    const XMLObject z2Param = nextTest.getChild(zpLoc);
    ParameterList pList;
    int numParams = z2Param.numChildren();
    for (int i=0; i < numParams; i++){
      const XMLObject param = z2Param.getChild(i);
      const string &pname = param.getTag();
      const string &pvalue = param.getAttribute(string("value"));
      // Problem:  For some of the parameters it is probably
      // illegal to set the value as a string.  We should
      // fix this.
      pList.set(pname, pvalue);
    }

    pList.set("compute_metrics", "yes");

    ////
    // Create the input.
    ////

    const XMLObject testParam = nextTest.getChild(tpLoc);

    int iaLoc = testParam.findFirstChild("inputAdapter");
    if (iaLoc < 0){
      ERRMSG("Driver error: no inputAdapter tag");
      continue;
    }

    const XMLObject iaNode = testParam.getChild(iaLoc);
    if (!iaNode.hasAttribute("name")){
      ERRMSG("Driver error: no \"name=\" attribute for inputAdapter");
      continue;
    }

    const string &iaType = iaNode.getAttribute(string("name"));

    tcrsMatrix_t *M = NULL;
    tMVector_t *objectWeights = NULL;
    tMVector_t *edgeWeights = NULL;
    tMVector_t *objectCoordinates = NULL;

    try{
      createInput(testParam, M, objectWeights, edgeWeights,objectCoordinates);
    }
    catch(exception &e){
      ostringstream msg;
      msg << "Driver error: can't create input for test " << tName;
      EXC_ERRMSG(msg.str(), e);
      continue; 
    }

#ifdef OMIT_UNTIL_CREATEINPUT_IS_DONE

    size_t numRows = M->getNodeNumRows();
    int vWeightDim = objectWeights->getNumVectors();
    int eWeightDim = edgeWeights->getNumVectors();
    int coordDim = objectCoordinates->getNumVectors();

    const RCP<const tmap_t> rowMap = M->getRowMap();
    const RCP<const tmap_t> colMap = M->getColMap();

    const gno_t *gids = rowMap->getNodeElementList().getRawPtr();

    vector<const scalar_t *> coords;
    for (int i=0; i < coordDim; i++)
      coords.push_back(objectCoordinates->getData(i).getRawPtr());

    vector<const scalar_t *> vweights;
    for (int i=0; i < vWeightDim; i++)
      coords.push_back(objectWeights->getData(i).getRawPtr());

    vector<const scalar_t *> eweights;
    for (int i=0; i < eWeightDim; i++)
      coords.push_back(objectWeights->getData(i).getRawPtr());

    std::vector<int> strides;  // empty vector implies all are "1"

    ////
    // Create the requested input adapter, run the
    // algorithm, and obtain the quality metrics,
    // and evaluate success.
    ////

    bool fail = false;
    ArrayRCP<const Zoltan2::MetricValues<scalar_t> > quality;

    if (iaType == string("BasicCoordinateInput")){
      bci_t *ia = NULL;
      try{
        ia = new bci_t(numRows, gids, coords, strides, vweights, strides);
      }
      catch (std::exception &e){
        EXC_ERRMSG("Test error: InputAdapter build", e);
        fail = true;
      }

      int err = runAlgorithm<bci_t>(rank, ia, pList, quality);

      if (err != 0)
        fail = true;
    }
    else if (iaType == string("BasicIdentifierInput")){
      bii_t *ia = NULL;
      try{
        ia = new bii_t(numRows, gids, vweights, strides);
      }
      catch (std::exception &e){
        EXC_ERRMSG("Test error: InputAdapter build", e);
        fail = true;
      }

      int err = runAlgorithm<bii_t>(rank, ia, pList, quality);

      if (err != 0)
        fail = true;
    }
    else if (iaType == string("BasicVectorInput")){
      bvi_t *ia = NULL;
      // UserInputForTests can give us a vector based on M.
      scalar_t *vec=NULL; 
      try{
        ia = new bvi_t(numRows, gids, vec, 1, vweights, strides);
      }
      catch (std::exception &e){
        EXC_ERRMSG("Test error: InputAdapter build", e);
        fail = true;
      }

      int err = runAlgorithm<bvi_t>(rank, ia, pList, quality);
      if (err != 0)
        fail = true;
    }
    else if (iaType == string("XpetraCrsGraphInput")){
      const RCP<const tcrsGraph_t> graph = M->getCrsGraph();
      xgi_t *ia = NULL;
      try{
        ia = new xgi_t(graph, vweights, strides, eweights, strides, 
          coords, strides);
      }
      catch (std::exception &e){
        EXC_ERRMSG("Test error: InputAdapter build", e);
        fail = true;
      }

      int err = runAlgorithm<xgi_t>(rank, ia, pList, quality);
      if (err != 0)
        fail = true;
    }
    else if (iaType == string("XpetraCrsMatrixInput")){
      xmi_t *ia = NULL;
      RCP<const tcrsMatrix_t> Mptr(M);
      try{
        ia = new xmi_t(Mptr, coordDim, vWeightDim);
      }
      catch (std::exception &e){
        EXC_ERRMSG("Test error: InputAdapter build", e);
        fail = true;
      }

      if (coordDim > 0){
        for (int i=0; i < coordDim; i++){
          ia->setRowCoordinates(i, coords[i], 1);
        }
      }

      if (vWeightDim > 0){
        for (int i=0; i < vWeightDim; i++){
          ia->setRowWeights(i, vweights[i], 1);
        }
      }

      int err = runAlgorithm<xmi_t>(rank, ia, pList, quality);
      if (err != 0)
        fail = true;
    }
    else if (iaType == string("XpetraMultiVectorInput")){
      // UserInputForTests can give us a multivector based on M.
      const RCP<const tMVector_t> MV;
      xmvi_t *ia = NULL;
      try{
        ia = new xmvi_t(MV, vweights, strides);
      }
      catch (std::exception &e){
        EXC_ERRMSG("Test error: InputAdapter build", e);
        fail = true;
      }

      int err = runAlgorithm<xmvi_t>(rank, ia, pList, quality);
      if (err != 0)
        fail = true;
    }
    else if (iaType == string("XpetraVectorInput")){
      // UserInputForTests can give us a vector based on M.
      const RCP<const tVector_t> V;
      xvi_t *ia = NULL;
      try{
        ia = new xvi_t(V, vweights, strides);
      }
      catch (std::exception &e){
        EXC_ERRMSG("Test error: InputAdapter build", e);
        fail = true;
      }

      int err = runAlgorithm<xvi_t>(rank, ia, pList, quality);
      if (err != 0)
        fail = true;
    }
    else{
      if (rank == 0)
        cout << "Driver error: invalid input adapter name"<< endl;
      fail = true;
    }

    if (!fail){

      if (evaluateSuccess(testParam, quality) == 0)
        numPass++;
    }
#endif

  }  // next test

  if (rank == 0)
    cout << "PASS" << endl;

  return 0;
}




