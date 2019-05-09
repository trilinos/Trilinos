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
#include <Zoltan2_ModelHelpers.hpp>
#include <Zoltan2_PamgenMeshAdapter.hpp>

// Teuchos includes
#include "Teuchos_XMLParameterListHelpers.hpp"

// Pamgen includes
#include "create_inline_mesh.h"

using Teuchos::RCP;

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef Zoltan2::BasicUserTypes<double> basic_user_t;



/*****************************************************************************/
/******************************** MAIN ***************************************/
/*****************************************************************************/

int main(int narg, char *arg[]) {

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > CommT = Tpetra::getDefaultComm();

  int me = CommT->getRank();
  int numProcs = CommT->getSize();

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
  Teuchos::ParameterList inputMeshList;

  if(xmlMeshInFileName.length()) {
    if (me == 0) {
      std::cout << "\nReading parameter list from the XML file \""
                <<xmlMeshInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile(xmlMeshInFileName,
                                         Teuchos::inoutArg(inputMeshList));
    if (me == 0) {
      inputMeshList.print(std::cout,2,true,true);
      std::cout << "\n";
    }
  }
  else {
    std::cout << "Cannot read input file: " << xmlMeshInFileName << "\n";
    return 5;
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

  if (me == 0) std::cout << "Generating mesh ... \n\n";

  // Generate mesh with Pamgen
  long long maxInt = std::numeric_limits<long long>::max();
  Create_Pamgen_Mesh(meshInput.c_str(), dim, me, numProcs, maxInt);

  // Creating mesh adapter
  if (me == 0) std::cout << "Creating mesh adapter ... \n\n";

  typedef Zoltan2::PamgenMeshAdapter<basic_user_t> inputAdapter_t;
  typedef inputAdapter_t::base_adapter_t base_adapter_t;

  inputAdapter_t ia(*CommT, "region");
  inputAdapter_t ia2(*CommT, "vertex");
  inputAdapter_t::gno_t const *adjacencyIds=NULL;
  inputAdapter_t::offset_t const *offsets=NULL;
  Teuchos::ArrayRCP<inputAdapter_t::offset_t> moffsets;
  Teuchos::ArrayRCP<inputAdapter_t::gno_t> madjacencyIds;
  ia.print(me);

  if (me == 0) std::cout << "REGION-BASED TEST" << std::endl;
  Zoltan2::MeshEntityType primaryEType = ia.getPrimaryEntityType();
  Zoltan2::MeshEntityType adjEType = ia.getAdjacencyEntityType();
  Zoltan2::MeshEntityType secondAdjEType = ia.getSecondAdjacencyEntityType();
  RCP<const base_adapter_t> baseInputAdapter;
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(CommT));
  std::bitset<Zoltan2::NUM_MODEL_FLAGS> modelFlags;

  if (ia.availAdjs(primaryEType, adjEType)) {
    // Create a GraphModel based on this input data.

    if (me == 0) std::cout << "        Creating GraphModel" << std::endl;

    baseInputAdapter = (rcp(dynamic_cast<const base_adapter_t *>(&ia), false));

    Zoltan2::GraphModel<base_adapter_t> graphModel(baseInputAdapter, env,
                                                   CommT, modelFlags);

    if (me == 0)
      std::cout << "        Calling get2ndAdjsViewFromAdjs" << std::endl;
    Zoltan2::get2ndAdjsViewFromAdjs(baseInputAdapter, graphModel.getComm(),
                                    primaryEType,
                                      secondAdjEType, moffsets, madjacencyIds);

    if (me == 0) std::cout << "        Checking results" << std::endl;
    if (ia.avail2ndAdjs(primaryEType, secondAdjEType)) {
      ia.get2ndAdjsView(primaryEType, secondAdjEType, offsets, adjacencyIds);
    }
    else{
      std::cout << "2nd adjacencies not available; cannot check results"
                << std::endl;
      return 2;
    }

    for (size_t telct = 0; telct < ia.getLocalNumOf(primaryEType); telct++) {
      if (offsets[telct+1]-offsets[telct]!=moffsets[telct+1]-moffsets[telct]) {
        std::cout << "Number of adjacencies do not match" << std::endl;
        return 3;
      }

      for (inputAdapter_t::offset_t j=moffsets[telct]; j<moffsets[telct+1]; j++) {
        ssize_t in_list = -1;

        for (inputAdapter_t::offset_t k=offsets[telct]; k<offsets[telct+1]; k++) {
          if (adjacencyIds[k] == madjacencyIds[j]) {
            in_list = k;
            break;
          }
        }

        if (in_list < 0) {
          std::cout << "Adjacency missing" << std::endl;
          return 4;
        }
      }
    }
  }
  else{
    std::cout << "Adjacencies not available" << std::endl;
    return 1;
  }

  if (me == 0) std::cout << "VERTEX-BASED TEST" << std::endl;
  primaryEType = ia2.getPrimaryEntityType();
  adjEType = ia2.getAdjacencyEntityType();
  secondAdjEType = ia2.getSecondAdjacencyEntityType();

  if (ia2.availAdjs(primaryEType, adjEType)) {
    if (ia2.avail2ndAdjs(primaryEType, secondAdjEType)) {
      ia2.get2ndAdjsView(primaryEType, secondAdjEType, offsets, adjacencyIds);
    }
    else{
      std::cout << "2nd adjacencies not available" << std::endl;
      return 2;
    }

    // Create a GraphModel based on this input data.

    if (me == 0) std::cout << "        Creating GraphModel" << std::endl;

    baseInputAdapter = (rcp(dynamic_cast<const base_adapter_t *>(&ia2),false));

    Zoltan2::GraphModel<base_adapter_t> graphModel2(baseInputAdapter, env,
                                                   CommT, modelFlags);

    if (me == 0)
      std::cout << "        Calling get2ndAdjsViewFromAdjs" << std::endl;
    Zoltan2::get2ndAdjsViewFromAdjs(baseInputAdapter, graphModel2.getComm(),
                                    primaryEType,
                                    secondAdjEType, moffsets, madjacencyIds);

    if (me == 0) std::cout << "        Checking results" << std::endl;

    for (size_t tnoct = 0; tnoct < ia2.getLocalNumOf(primaryEType); tnoct++) {
      if (offsets[tnoct+1]-offsets[tnoct]!=moffsets[tnoct+1]-moffsets[tnoct]) {
        std::cout << "Number of adjacencies do not match" << std::endl;
        return 3;
      }

      for (inputAdapter_t::offset_t j=moffsets[tnoct]; j<moffsets[tnoct+1]; j++) {
        ssize_t in_list = -1;

        for (inputAdapter_t::offset_t k=offsets[tnoct]; k<offsets[tnoct+1]; k++) {
          if (adjacencyIds[k] == madjacencyIds[j]) {
            in_list = k;
            break;
          }
        }

        if (in_list < 0) {
          std::cout << "Adjacency missing" << std::endl;
          return 4;
        }
      }
    }
  }
  else{
    std::cout << "Adjacencies not available" << std::endl;
    return 1;
  }

  // delete mesh
  if (me == 0) std::cout << "Deleting the mesh ... \n\n";

  Delete_Pamgen_Mesh();

  if (me == 0)
    std::cout << "PASS" << std::endl;

  return 0;
}
/*****************************************************************************/
/********************************* END MAIN **********************************/
/*****************************************************************************/
