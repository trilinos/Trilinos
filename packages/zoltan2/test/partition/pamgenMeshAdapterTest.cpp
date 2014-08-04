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

//Tpetra includes
#include "Tpetra_DefaultPlatform.hpp"

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
//Tpetra typedefs
typedef Tpetra::DefaultPlatform::DefaultPlatformType            Platform;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;



/*****************************************************************************/
/******************************** MAIN ***************************************/
/*****************************************************************************/

int main(int argc, char *argv[]) {


  //Get the default communicator and node for Tpetra
  //rewrite using with IFDEF for MPI/no MPI??
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > CommT = platform.getComm();
  int MyPID = CommT->getRank();

  /***************************************************************************/
  /*************************** GET XML INPUTS ********************************/
  /***************************************************************************/

  // Command line for xml file, otherwise use default
  std::string   xmlMeshInFileName;
  if(argc>=2) xmlMeshInFileName=string(argv[1]);
  else xmlMeshInFileName="Poisson.xml";

  // Read xml file into parameter list
  ParameterList inputMeshList;

  if(xmlMeshInFileName.length()) {
    if (MyPID == 0) {
      std::cout << "\nReading parameter list from the XML file \""
		<<xmlMeshInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile(xmlMeshInFileName, 
					 Teuchos::inoutArg(&inputMeshList));
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
  /***************************** GENERATE MESH *******************************/
  /***************************************************************************/

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);

  // delete mesh
  Delete_Pamgen_Mesh();

  return 0;
}
/*****************************************************************************/
/********************************* END MAIN **********************************/
/*****************************************************************************/
