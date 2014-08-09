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

/*! \file pamgenMeshAdapterTest.cpp
    \brief An example of partitioning pamgen coordinates with RCB.

    \author Created by V. Leung.

*/

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

//#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_PamgenMeshAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

//Tpetra includes
#include "Tpetra_DefaultPlatform.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Pamgen includes
#include "create_inline_mesh.h"

using namespace std;
using Teuchos::ParameterList;
using Teuchos::RCP;

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
//Tpetra typedefs
typedef Tpetra::DefaultPlatform::DefaultPlatformType            Platform;
typedef Tpetra::MultiVector<double, int, int>     tMVector_t;



/*****************************************************************************/
/******************************** MAIN ***************************************/
/*****************************************************************************/

int main(int narg, char *arg[]) {

  int numProcs=1;
  int rank=0;

  Teuchos::GlobalMPISession mpiSession(&narg, &arg,0);
  rank=mpiSession.getRank();
  numProcs=mpiSession.getNProc();


  //Get the default communicator and node for Tpetra
  //rewrite using with IFDEF for MPI/no MPI??
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > CommT = platform.getComm();
  int MyPID = CommT->getRank();

  if (MyPID == 0){
  cout 
    << "====================================================================\n" 
    << "|                                                                  |\n" 
    << "|          Example: Partition Pamgen Hexahedral Mesh               |\n" 
    << "|                                                                  |\n"
    << "|  Questions? Contact  Karen Devine      (kddevin@sandia.gov),     |\n"
    << "|                      Erik Boman        (egboman@sandia.gov),     |\n"
    << "|                      Siva Rajamanickam (srajama@sandia.gov).     |\n"
    << "|                                                                  |\n"
    << "|  Pamgen's website:   http://trilinos.sandia.gov/packages/pamgen  |\n"
    << "|  Zoltan2's website:  http://trilinos.sandia.gov/packages/zoltan2 |\n"
    << "|  Trilinos website:   http://trilinos.sandia.gov                  |\n"
    << "|                                                                  |\n"
    << "====================================================================\n";
  }


#ifdef HAVE_MPI
  if (MyPID == 0) {
    cout << "PARALLEL executable \n";
  }
#else
  if (MyPID == 0) {
    cout << "SERIAL executable \n";
  }
#endif

  /***************************************************************************/
  /*************************** GET XML INPUTS ********************************/
  /***************************************************************************/

  // default values for command-line arguments
  std::string xmlMeshInFileName("Poisson.xml");
  std::string method("rcb");
  int nParts = CommT->getSize();

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("xmlfile", &xmlMeshInFileName,
                 "XML file with PamGen specifications");
  cmdp.setOption("method", &method,
                 "Partitioning method to use:  rcb or graph");
  cmdp.setOption("nparts", &nParts,
                 "Number of parts to create");
  cmdp.parse(narg, arg);

  // Read xml file into parameter list
  ParameterList inputMeshList;

  if(xmlMeshInFileName.length()) {
    if (MyPID == 0) {
      cout << "\nReading parameter list from the XML file \""
		<<xmlMeshInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile(xmlMeshInFileName, 
					 Teuchos::inoutArg(inputMeshList));
    if (MyPID == 0) {
      inputMeshList.print(cout,2,true,true);
      cout << "\n";
    }
  }
  else {
    cout << "Cannot read input file: " << xmlMeshInFileName << "\n";
    return 0;
  }

  // Get pamgen mesh definition
  std::string meshInput = Teuchos::getParameter<std::string>(inputMeshList,
							     "meshInput");


  /***************************************************************************/
  /********************** GET CELL TOPOLOGY **********************************/
  /***************************************************************************/

  // Get dimensions
  int dim = 3;

  /***************************************************************************/
  /***************************** GENERATE MESH *******************************/
  /***************************************************************************/

  if (MyPID == 0) {
    cout << "Generating mesh ... \n\n";
  }

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);

  typedef Zoltan2::PamgenMeshAdapter<tMVector_t> inputAdapter_t;

  inputAdapter_t ia;

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("imbalance_tolerance", 1.1);
  params.set("num_global_parts", nParts);

  if (method == "rcb") {
    params.set("algorithm", "rcb");
    params.set("rectilinear", "yes");
  }
  else if (method == "graph") {
    params.set("partitioning_approach", "partition");
    params.set("algorithm", "scotch");
  }

  Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &params, CommT);

  problem.solve();

  if (CommT->getRank() == 0)
    problem.printMetrics(cout);

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  // delete mesh
  Delete_Pamgen_Mesh();
  return 0;
}
/*****************************************************************************/
/********************************* END MAIN **********************************/
/*****************************************************************************/
