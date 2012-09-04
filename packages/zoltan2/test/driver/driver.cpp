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

    \todo There are no defaults.  Everything about the test has to
              be given in the XML input file.
*/

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_Parameters.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

using Teuchos::ParameterList;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::XMLObject;

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::exception;

typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;
typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;

void createInput(const XMLObject &testInfo, 
  tcrsMatrix_t *&M, tMVector_t *&objWgts, tMVector_t *&edgeWgts, 
  tMVector_t * &coords)
{
  M = NULL;
  objWgts = NULL;
  edgeWgts = NULL;
  coords = NULL;
}

#if 0
void createInputAdapter();

void createParameterList();

void evaluateSuccess();
#endif

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
    if (rank == 0){
      cerr << "Driver error: reading " << inFileName << endl;
    }
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
    if (rank == 0){
      cerr << "Driver error: xml file does not define \"Tests\"" << endl;
    }
  }

  const string testSuiteName = xmlIn.getAttribute(string("name"));
  int numTests = xmlIn.numChildren();

  if (rank == 0){
    cout << "Test suite: " << testSuiteName << ", ";
    cout << numTests << " tests" << endl;
  }

  for (int testNum=0; testNum < numTests; testNum++){

    // Get the next test.

    const XMLObject nextTest = xmlIn.getChild(testNum);
    const string tName = xmlIn.getAttribute(string("name"));

    int tpLoc = nextTest.findFirstChild("TestParameters");
    int zpLoc = nextTest.findFirstChild("Zoltan2Parameters");

    if (tpLoc < 0){
      if (rank == 0){
        cout << "Driver error: no TestParameters for " << tName << endl;
      }
      continue;
    }

    if (zpLoc < 0){
      if (rank == 0){
        cout << "Driver error: no Zoltan2Parameters for " << tName << endl;
      }
      continue;
    }

    const XMLObject z2Param = nextTest.getChild(zpLoc);

    // Create the Zoltan2 parameter list.

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

    // Create the specified input.

    const XMLObject testParam = nextTest.getChild(tpLoc);

    tcrsMatrix_t *M = NULL;
    tMVector_t *objectWeights = NULL;
    tMVector_t *edgeWeights = NULL;
    tMVector_t *objectCoordinates = NULL;

    try{
      createInput(testParam, M, objectWeights, edgeWeights,objectCoordinates);
    }
    catch(exception &e){
      if (rank == 0){
        cout << "Driver error: can't create input for test " << tName << endl;
        cout << e.what() << endl;
      }
      continue; 
    }

    // Create the requested input adapter.

    int iaLoc = testParam.findFirstChild("inputAdapter");
    if (iaLoc < 0){
      if (rank == 0){
        cout << "Driver error: no inputAdapter for " << tName << endl;
      }
      continue;
    }

#if 0
    createInputAdapter();

    createParameterList();

    // Add "compute_metrics" to the parameter list.

    // Create PartitioningProblem.

    // Solve.

    // Get quality metrics.

    evaluateSuccess();

#endif
  }

  return 0;
}




