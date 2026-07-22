// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Basic testing of Zoltan2::PamgenMeshAdapter

#include "Zoltan2_PamgenMeshAdapter.hpp"
#include "Zoltan2_componentMetrics.hpp"

// Teuchos includes
#include "Teuchos_XMLParameterListHelpers.hpp"

// Pamgen includes
#include "create_inline_mesh.h"

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef Zoltan2::BasicUserTypes<> basic_user_t;

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
  long long maxInt = 9223372036854775807LL;
  Create_Pamgen_Mesh(meshInput.c_str(), dim, me, numProcs, maxInt);

  // Creating mesh adapter
  if (me == 0) std::cout << "Creating mesh adapter ... \n\n";

  typedef Zoltan2::PamgenMeshAdapter<basic_user_t> inputAdapter_t;

  inputAdapter_t ia(*CommT, "region");
  inputAdapter_t ia2(*CommT, "vertex");
  inputAdapter_t::gno_t const *adjacencyIds=NULL;
  inputAdapter_t::offset_t const *offsets=NULL;
  ia.print(me);

  // Exercise the componentMetrics on the input; make sure the adapter works
  {
    Zoltan2::perProcessorComponentMetrics<inputAdapter_t> cc(ia, *CommT);
    std::cout << me << " Region-based: Number of components on processor = "
              << cc.getNumComponents() << std::endl;
    std::cout << me << " Region-based: Max component size on processor = "
              << cc.getMaxComponentSize() << std::endl;
    std::cout << me << " Region-based: Min component size on processor = "
              << cc.getMinComponentSize() << std::endl;
    std::cout << me << " Region-based: Avg component size on processor = "
              << cc.getAvgComponentSize() << std::endl;
  }

  Zoltan2::MeshEntityType primaryEType = ia.getPrimaryEntityType();
  Zoltan2::MeshEntityType adjEType = ia.getAdjacencyEntityType();

  int dimension, num_nodes, num_elem;
  int error = 0;
  char title[100];
  int exoid = 0;
  int num_elem_blk, num_node_sets, num_side_sets;
  error += im_ex_get_init(exoid, title, &dimension, &num_nodes, &num_elem,
                          &num_elem_blk, &num_node_sets, &num_side_sets);

  int *element_num_map = new int [num_elem];
  error += im_ex_get_elem_num_map(exoid, element_num_map);

  int *node_num_map = new int [num_nodes];
  error += im_ex_get_node_num_map(exoid, node_num_map);

  int *elem_blk_ids = new int [num_elem_blk];
  error += im_ex_get_elem_blk_ids(exoid, elem_blk_ids);

  int *num_nodes_per_elem = new int [num_elem_blk];
  int *num_attr           = new int [num_elem_blk];
  int *num_elem_this_blk  = new int [num_elem_blk];
  char **elem_type        = new char * [num_elem_blk];
  int **connect           = new int * [num_elem_blk];

  for(int i = 0; i < num_elem_blk; i++){
    elem_type[i] = new char [MAX_STR_LENGTH + 1];
    error += im_ex_get_elem_block(exoid, elem_blk_ids[i], elem_type[i],
                                  (int*)&(num_elem_this_blk[i]),
                                  (int*)&(num_nodes_per_elem[i]),
                                  (int*)&(num_attr[i]));
    delete[] elem_type[i];
  }

  delete[] elem_type;
  elem_type = NULL;
  delete[] num_attr;
  num_attr = NULL;

  for(int b = 0; b < num_elem_blk; b++) {
    connect[b] = new int [num_nodes_per_elem[b]*num_elem_this_blk[b]];
    error += im_ex_get_elem_conn(exoid, elem_blk_ids[b], connect[b]);
  }

  delete[] elem_blk_ids;
  elem_blk_ids = NULL;

  int telct = 0;

  if (ia.availAdjs(primaryEType, adjEType)) {
    ia.getAdjsView(primaryEType, adjEType, offsets, adjacencyIds);

    if ((int)ia.getLocalNumOf(primaryEType) != num_elem) {
      std::cout << "Number of elements do not match\n";
      return 2;
    }

    for (int b = 0; b < num_elem_blk; b++) {
      for (int i = 0; i < num_elem_this_blk[b]; i++) {
        if (offsets[telct + 1] - offsets[telct] != num_nodes_per_elem[b]) {
          std::cout << "Number of adjacencies do not match" << std::endl;
          return 3;
        }

        for (int j = 0; j < num_nodes_per_elem[b]; j++) {
          ssize_t in_list = -1;

          for(inputAdapter_t::offset_t k=offsets[telct];k<offsets[telct+1];k++){
            if(adjacencyIds[k] ==
               node_num_map[connect[b][i*num_nodes_per_elem[b]+j]-1]) {
              in_list = k;
              break;
            }
          }

          if (in_list < 0) {
            std::cout << "Adjacency missing" << std::endl;
            return 4;
          }
        }

        ++telct;
      }
    }

    if (telct != num_elem) {
      std::cout << "Number of elements do not match\n";
      return 2;
    }
  }
  else{
    std::cout << "Adjacencies not available" << std::endl;
    return 1;
  }

  ia2.print(me);
  {
    Zoltan2::perProcessorComponentMetrics<inputAdapter_t> cc(ia2, *CommT);
    std::cout << me << " Vertex-based: Number of components on processor = "
              << cc.getNumComponents() << std::endl;
    std::cout << me << " Vertex-based: Max component size on processor = "
              << cc.getMaxComponentSize() << std::endl;
    std::cout << me << " Vertex-based: Min component size on processor = "
              << cc.getMinComponentSize() << std::endl;
    std::cout << me << " Vertex-based: Avg component size on processor = "
              << cc.getAvgComponentSize() << std::endl;
  }

  primaryEType = ia2.getPrimaryEntityType();
  adjEType = ia2.getAdjacencyEntityType();

  if (ia2.availAdjs(primaryEType, adjEType)) {
    ia2.getAdjsView(primaryEType, adjEType, offsets, adjacencyIds);

    if ((int)ia2.getLocalNumOf(primaryEType) != num_nodes) {
      std::cout << "Number of nodes do not match\n";
      return 2;
    }

    telct = 0;
    int *num_adj = new int[num_nodes];

    for (int i = 0; i < num_nodes; i++) {
      num_adj[i] = 0;
    }

    for (int b = 0; b < num_elem_blk; b++) {
      for (int i = 0; i < num_elem_this_blk[b]; i++) {
        for (int j = 0; j < num_nodes_per_elem[b]; j++) {
          ssize_t in_list = -1;
          ++num_adj[connect[b][i * num_nodes_per_elem[b] + j] - 1];

          for(inputAdapter_t::lno_t k =
                offsets[connect[b][i * num_nodes_per_elem[b] + j] - 1];
              k < offsets[connect[b][i * num_nodes_per_elem[b] + j]]; k++) {
            if(adjacencyIds[k] == element_num_map[telct]) {
              in_list = k;
              break;
            }
          }

          if (in_list < 0) {
            std::cout << "Adjacency missing" << std::endl;
            return 4;
          }
        }

        ++telct;
      }
    }

    for (int i = 0; i < num_nodes; i++) {
      if (offsets[i + 1] - offsets[i] != num_adj[i]) {
        std::cout << "Number of adjacencies do not match" << std::endl;
        return 3;
      }
    }

    delete[] num_adj;
    num_adj = NULL;
  }
  else{
    std::cout << "Adjacencies not available" << std::endl;
    return 1;
  }

  // delete mesh
  if (me == 0) std::cout << "Deleting the mesh ... \n\n";

  Delete_Pamgen_Mesh();
  // clean up - reduce the result codes

  // make sure another process doesn't mangle the PASS output or the test will read as a fail when it should pass
  std::cout << std::flush;
  CommT->barrier();
  if (me == 0)
    std::cout << "PASS" << std::endl;

  return 0;
}
/*****************************************************************************/
/********************************* END MAIN **********************************/
/*****************************************************************************/
