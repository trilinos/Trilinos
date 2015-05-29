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
// Basic testing of Zoltan2::PamgenMeshAdapter

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_PamgenMeshAdapter.hpp>

// Teuchos includes
#include "Teuchos_XMLParameterListHelpers.hpp"

// Pamgen includes
#include "create_inline_mesh.h"

using namespace std;
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

  Teuchos::GlobalMPISession mpiSession(&narg, &arg,0);
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > CommT = platform.getComm();

  int me = CommT->getRank();
  int numProcs = CommT->getSize();
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags;

  /***************************************************************************/
  /*************************** GET XML INPUTS ********************************/
  /***************************************************************************/

  // default values for command-line arguments
  std::string xmlMeshInFileName("Poisson.xml");

  // Read run-time options.
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("xmlfile", &xmlMeshInFileName,
                 "XML file with PamGen specifications");
  cmdp.parse(narg, arg);

  // Read xml file into parameter list
  ParameterList inputMeshList;

  if(xmlMeshInFileName.length()) {
    if (me == 0) {
      cout << "\nReading parameter list from the XML file \""
		<<xmlMeshInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile(xmlMeshInFileName, 
					 Teuchos::inoutArg(inputMeshList));
    if (me == 0) {
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

  if (me == 0) cout << "Generating mesh ... \n\n";

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  Create_Pamgen_Mesh(meshInput.c_str(), dim, me, numProcs, maxInt);

  // Creating mesh adapter
  if (me == 0) cout << "Creating mesh adapter ... \n\n";

  typedef Zoltan2::PamgenMeshAdapter<tMVector_t> inputAdapter_t;
  typedef inputAdapter_t::base_adapter_t base_adapter_t;

  inputAdapter_t ia(*CommT, "region");
  inputAdapter_t::zgid_t const *adjacencyIds=NULL;
  inputAdapter_t::lno_t const *offsets=NULL;
  ia.print(me);
  Zoltan2::MeshEntityType primaryEType = ia.getPrimaryEntityType();
  Zoltan2::MeshEntityType adjEType = ia.getAdjacencyEntityType();
  Zoltan2::MeshEntityType secondAdjEType = ia.getSecondAdjacencyEntityType();

  if (ia.availAdjs(primaryEType, adjEType)) {
    if (ia.avail2ndAdjs(primaryEType, secondAdjEType)) {
      ia.get2ndAdjsView(primaryEType, secondAdjEType, offsets, adjacencyIds);
    }
    else{
      std::cout << "2nd adjacencies not available" << std::endl;
      return 2;
    }

    // Create a GraphModel based on this input data.

    if (me == 0) std::cout << "        Creating GraphModel" << std::endl;

    /*Zoltan2::GraphModel<base_adapter_t>
      model(dynamic_cast<base_adapter_t *>(&ia), env, CommT, modelFlags);*/

    
  }
  else{
    std::cout << "Adjacencies not available" << std::endl;
    return 1;
  }

  // delete mesh
  if (me == 0) cout << "Deleting the mesh ... \n\n";

  Delete_Pamgen_Mesh();

  if (me == 0)
    std::cout << "PASS" << std::endl;

  return 0;
}
/*****************************************************************************/
/********************************* END MAIN **********************************/
/*****************************************************************************/
