// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

// TrilinosCouplings includes
#include <TrilinosCouplings_config.h>

// Intrepid includes
#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_ArrayTools.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid_RealSpaceTools.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_Utils.hpp>

// Teuchos includes
#include <Teuchos_TimeMonitor.hpp>

// Shards includes
#include <Shards_CellTopology.hpp>

// Pamgen includes
#include <create_inline_mesh.h>
#include <im_exodusII_l.h>
#include <im_ne_nemesisI_l.h>
#include <pamgen_extras.h>

// Sacado includes
#include <Sacado.hpp>

// My includes
#include "TrilinosCouplings_EpetraIntrepidPoissonExample.hpp"
#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"


namespace TrilinosCouplings {
namespace EpetraIntrepidPoissonExample {

/**********************************************************************************/
/******** FUNCTION DECLARATIONS FOR EXACT SOLUTION AND SOURCE TERMS ***************/
/**********************************************************************************/

/** \brief  User-defined exact solution.

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \return Value of the exact solution at (x,y,z)
 */
template<typename Scalar>
const Scalar
exactSolution (const Scalar& x, const Scalar& y, const Scalar& z);

/** \brief  User-defined material tensor.

    \param  material    [out]   3 x 3 material tensor evaluated at (x,y,z)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \warning Symmetric and positive definite tensor is required for every (x,y,z).
*/
template<typename Scalar>
void
materialTensor (Scalar material[][3],
                const Scalar& x,
                const Scalar& y,
                const Scalar& z);

/** \brief  Computes gradient of the exact solution. Requires user-defined exact solution.

    \param  gradExact  [out]   gradient of the exact solution evaluated at (x,y,z)
    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point
 */
template<typename Scalar>
void
exactSolutionGrad (Scalar gradExact[3],
                   const Scalar& x,
                   const Scalar& y,
                   const Scalar& z);

/** \brief Computes source term: f = -div(A.grad u).  Requires user-defined exact solution
           and material tensor.

    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point

    \return Source term corresponding to the user-defined exact solution evaluated at (x,y,z)
 */
template<typename Scalar>
const Scalar
sourceTerm (Scalar& x, Scalar& y, Scalar& z);

/** \brief Compute the material tensor at array of points in physical space.

    \param worksetMaterialValues      [out]     Rank-2, 3 or 4 array with dimensions (C,P), (C,P,D) or (C,P,D,D)
                                                with the values of the material tensor
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void
evaluateMaterialTensor (ArrayOut& worksetMaterialValues,
                        const ArrayIn& evaluationPoints);

/** \brief Compute the source term at array of points in physical space.

    \param sourceTermValues           [out]     Rank-2 (C,P) array with the values of the source term
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void
evaluateSourceTerm (ArrayOut& sourceTermValues,
                    const ArrayIn& evaluationPoints);

/** \brief Compute the exact solution at array of points in physical space.

    \param exactSolutionValues        [out]     Rank-2 (C,P) array with the values of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void
evaluateExactSolution (ArrayOut& exactSolutionValues,
                       const ArrayIn& evaluationPoints);

/** \brief Computation of the gradient of the exact solution at array of points in physical space.

    \param exactSolutionGradValues    [out]     Rank-3 (C,P,D) array with the values of the gradient of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void
evaluateExactSolutionGrad (ArrayOut& exactSolutionGradValues,
                           const ArrayIn& evaluationPoints);

void
makeMatrixAndRightHandSide (Teuchos::RCP<sparse_matrix_type>& A,
                            Teuchos::RCP<multivector_type>& B,
                            Teuchos::RCP<multivector_type>& X_exact,
                            Teuchos::RCP<multivector_type>& X,
                            const Teuchos::RCP<const Epetra_Comm>& comm,
                            const std::string& meshInput,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose,
                            const bool debug)
{
  using Teuchos::RCP;
  using Teuchos::rcp_implicit_cast;

  RCP<vector_type> b, x_exact, x;
  makeMatrixAndRightHandSide (A, b, x_exact, x, comm, meshInput,
                              out, err, verbose, debug);

  B = rcp_implicit_cast<multivector_type> (b);
  X_exact = rcp_implicit_cast<multivector_type> (x_exact);
  X = rcp_implicit_cast<multivector_type> (x);
}

void
makeMatrixAndRightHandSide (Teuchos::RCP<sparse_matrix_type>& A,
                            Teuchos::RCP<vector_type>& B,
                            Teuchos::RCP<vector_type>& X_exact,
                            Teuchos::RCP<vector_type>& X,
                            const Teuchos::RCP<const Epetra_Comm>& comm,
                            const std::string& meshInput,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose,
                            const bool debug)
{
  using namespace Intrepid;
  using Teuchos::arcp;
  using Teuchos::arcp_const_cast;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::arrayView;
  using Teuchos::as;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::TimeMonitor;
  using std::endl;
  typedef Teuchos::ArrayView<LO>::size_type size_type;
  typedef Teuchos::ScalarTraits<ST> STS;

  (void) verbose;
  (void) debug;

  typedef Epetra_Map      map_type;
  typedef Epetra_Export   export_type;
  typedef Epetra_Import   import_type;
  typedef Epetra_CrsGraph sparse_graph_type;

  // Number of independent variables fixed at 3
  typedef Sacado::Fad::SFad<ST, 3>     Fad3;
  typedef Intrepid::FunctionSpaceTools IntrepidFSTools;
  typedef Intrepid::RealSpaceTools<ST> IntrepidRSTools;
  typedef Intrepid::CellTools<ST>      IntrepidCTools;

  const int myRank = comm->MyPID ();
  const int numProcs = comm->NumProc ();
  // We'll use this to check return values of Epetra methods.
  int errCode = 0;

  *out << "makeMatrixAndRightHandSide:" << endl;
  Teuchos::OSTab tab (out);

  /**********************************************************************************/
  /***************************** GET CELL TOPOLOGY **********************************/
  /**********************************************************************************/

  *out << "Getting cell topology" << endl;

  // Get cell topology for base hexahedron
  shards::CellTopology cellType (shards::getCellTopologyData<shards::Hexahedron<8> > ());

  // Get dimensions
  int numNodesPerElem = cellType.getNodeCount();
  int spaceDim = cellType.getDimension();
  int dim = 3;

  /**********************************************************************************/
  /******************************* GENERATE MESH ************************************/
  /**********************************************************************************/

  *out << "Generating mesh" << endl;

  int error = 0; // Number of errors in generating the mesh

  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  Create_Pamgen_Mesh (meshInput.c_str (), dim, myRank, numProcs, maxInt);

  std::string msg ("Poisson: ");

  // Get mesh size info
  char title[100];
  long long numDim;
  long long numNodes;
  long long numElems;
  long long numElemBlk;
  long long numNodeSets;
  long long numSideSets;
  int id = 0;

  im_ex_get_init_l (id, title, &numDim, &numNodes, &numElems, &numElemBlk,
                    &numNodeSets, &numSideSets);

  TEUCHOS_TEST_FOR_EXCEPTION(numElems == 0, std::runtime_error,
    "The number of elements in the mesh is zero.");

  long long numNodesGlobal;
  long long numElemsGlobal;
  long long numElemBlkGlobal;
  long long numNodeSetsGlobal;
  long long numSideSetsGlobal;

  im_ne_get_init_global_l (id, &numNodesGlobal, &numElemsGlobal,
                           &numElemBlkGlobal, &numNodeSetsGlobal,
                           &numSideSetsGlobal);

  // Print mesh information
  {
    Teuchos::OSTab tab2 (out);
    *out << "Global number of elements:                     "
         << numElemsGlobal << endl
         << "Global number of nodes (incl. boundary nodes): "
         << numNodesGlobal << endl;
  }

  long long * block_ids = new long long [numElemBlk];
  error += im_ex_get_elem_blk_ids_l(id, block_ids);

  long long  *nodes_per_element   = new long long [numElemBlk];
  long long  *element_attributes  = new long long [numElemBlk];
  long long  *elements            = new long long [numElemBlk];
  char      **element_types       = new char * [numElemBlk];
  long long **elmt_node_linkage   = new long long * [numElemBlk];

  for (long long i = 0; i < numElemBlk; ++i) {
    element_types[i] = new char [MAX_STR_LENGTH + 1];
    error += im_ex_get_elem_block_l (id,
                                     block_ids[i],
                                     element_types[i],
                                     (long long*)&(elements[i]),
                                     (long long*)&(nodes_per_element[i]),
                                     (long long*)&(element_attributes[i]));
  }

  // connectivity
  for (long long b = 0; b < numElemBlk; ++b) {
    elmt_node_linkage[b] =  new long long [nodes_per_element[b]*elements[b]];
    error += im_ex_get_elem_conn_l (id,block_ids[b], elmt_node_linkage[b]);
  }

  // Get node-element connectivity
  int telct = 0;
  FieldContainer<int> elemToNode(numElems,numNodesPerElem);
  for (long long b = 0; b < numElemBlk; b++) {
    for (long long el = 0; el < elements[b]; el++) {
      for (int j = 0; j < numNodesPerElem; ++j) {
        elemToNode(telct,j) = elmt_node_linkage[b][el*numNodesPerElem + j]-1;
      }
      ++telct;
    }
  }

  // Read node coordinates and place in field container
  FieldContainer<ST> nodeCoord (numNodes, dim);
  ST * nodeCoordx = new ST [numNodes];
  ST * nodeCoordy = new ST [numNodes];
  ST * nodeCoordz = new ST [numNodes];
  im_ex_get_coord_l (id, nodeCoordx, nodeCoordy, nodeCoordz);
  for (int i=0; i<numNodes; i++) {
    nodeCoord(i,0)=nodeCoordx[i];
    nodeCoord(i,1)=nodeCoordy[i];
    nodeCoord(i,2)=nodeCoordz[i];
  }
  delete [] nodeCoordx;
  delete [] nodeCoordy;
  delete [] nodeCoordz;

  // parallel info
  long long num_internal_nodes;
  long long num_border_nodes;
  long long num_external_nodes;
  long long num_internal_elems;
  long long num_border_elems;
  long long num_node_comm_maps;
  long long num_elem_comm_maps;
  im_ne_get_loadbal_param_l( id,
                             &num_internal_nodes,
                             &num_border_nodes,
                             &num_external_nodes,
                             &num_internal_elems,
                             &num_border_elems,
                             &num_node_comm_maps,
                             &num_elem_comm_maps,
                             0/*unused*/ );

  if (num_node_comm_maps > 0) {
    node_comm_proc_ids   = new long long  [num_node_comm_maps];
    node_cmap_node_cnts  = new long long  [num_node_comm_maps];
    node_cmap_ids        = new long long  [num_node_comm_maps];
    comm_node_ids        = new long long* [num_node_comm_maps];
    comm_node_proc_ids   = new long long* [num_node_comm_maps];

    long long *  elem_cmap_ids        = new long long [num_elem_comm_maps];
    long long *  elem_cmap_elem_cnts  = new long long [num_elem_comm_maps];

    if (im_ne_get_cmap_params_l (id,
                                 node_cmap_ids,
                                 (long long*)node_cmap_node_cnts,
                                 elem_cmap_ids,
                                 (long long*)elem_cmap_elem_cnts,
                                 0/*not used proc_id*/ ) < 0) {
      ++error;
    }

    for (long long j = 0; j < num_node_comm_maps; ++j) {
      comm_node_ids[j]       = new long long [node_cmap_node_cnts[j]];
      comm_node_proc_ids[j]  = new long long [node_cmap_node_cnts[j]];
      if (im_ne_get_node_cmap_l (id,
                                 node_cmap_ids[j],
                                 comm_node_ids[j],
                                 comm_node_proc_ids[j],
                                 0/*not used proc_id*/ ) < 0) {
        ++error;
      }
      node_comm_proc_ids[j] = comm_node_proc_ids[j][0];
    }

    delete [] elem_cmap_ids;
    delete [] elem_cmap_elem_cnts;
  }

  //
  // Calculate global node ids
  //
  Array<long long> globalNodeIds (numNodes);
  // nodeIsOwned must be a raw array, because std::vector<T> (and
  // therefore Teuchos::Array<T>) has a specialization for T = bool
  // that messes up the pointer type.
  bool* nodeIsOwned = new bool [numNodes];
  calc_global_node_ids (globalNodeIds.getRawPtr (),
                        nodeIsOwned,
                        numNodes,
                        num_node_comm_maps,
                        node_cmap_node_cnts,
                        node_comm_proc_ids,
                        comm_node_ids,
                        myRank);
  //
  // Mesh cleanup
  //
  delete [] block_ids;
  block_ids = NULL;
  delete [] nodes_per_element;
  nodes_per_element = NULL;
  delete [] element_attributes;
  element_attributes = NULL;
  for (long long b = 0; b < numElemBlk; ++b) {
    delete [] elmt_node_linkage[b];
    delete [] element_types[b];
  }
  delete [] element_types;
  element_types = NULL;
  delete [] elmt_node_linkage;
  elmt_node_linkage = NULL;
  if (num_node_comm_maps > 0) {
    delete [] node_comm_proc_ids;
    node_comm_proc_ids = NULL;
    delete [] node_cmap_node_cnts;
    node_cmap_node_cnts = NULL;
    delete [] node_cmap_ids;
    node_cmap_ids = NULL;
    for (long long i = 0; i < num_node_comm_maps; ++i) {
      delete [] comm_node_ids[i];
      delete [] comm_node_proc_ids[i];
    }
    delete [] comm_node_ids;
    comm_node_ids = NULL;
    delete [] comm_node_proc_ids;
    comm_node_proc_ids = NULL;
  }
  delete [] elements;
  elements = NULL;

  // Container indicating whether a node is on the boundary (1-yes 0-no)
  FieldContainer<int> nodeOnBoundary (numNodes);

  // Get boundary (side set) information
  long long * sideSetIds = new long long [numSideSets];
  long long numSidesInSet;
  long long numDFinSet;
  im_ex_get_side_set_ids_l(id,sideSetIds);
  for (int i=0; i < numSideSets; ++i) {
    im_ex_get_side_set_param_l (id,sideSetIds[i], &numSidesInSet, &numDFinSet);
    if (numSidesInSet > 0){
      long long * sideSetElemList = new long long [numSidesInSet];
      long long * sideSetSideList = new long long [numSidesInSet];
      im_ex_get_side_set_l (id, sideSetIds[i], sideSetElemList, sideSetSideList);
      for (int j = 0; j < numSidesInSet; ++j) {
        int sideNode0 = cellType.getNodeMap(2,sideSetSideList[j]-1,0);
        int sideNode1 = cellType.getNodeMap(2,sideSetSideList[j]-1,1);
        int sideNode2 = cellType.getNodeMap(2,sideSetSideList[j]-1,2);
        int sideNode3 = cellType.getNodeMap(2,sideSetSideList[j]-1,3);

        nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode0))=1;
        nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode1))=1;
        nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode2))=1;
        nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode3))=1;
      }
      delete [] sideSetElemList;
      delete [] sideSetSideList;
    }
  }
  delete [] sideSetIds;

/**********************************************************************************/
/********************************* GET CUBATURE ***********************************/
/**********************************************************************************/

  *out << "Getting cubature" << endl;

  // Get numerical integration points and weights
  DefaultCubatureFactory<ST>  cubFactory;
  int cubDegree = 2;
  RCP<Cubature<ST> > hexCub = cubFactory.create (cellType, cubDegree);

  int cubDim       = hexCub->getDimension ();
  int numCubPoints = hexCub->getNumPoints ();

  FieldContainer<ST> cubPoints (numCubPoints, cubDim);
  FieldContainer<ST> cubWeights (numCubPoints);

  hexCub->getCubature (cubPoints, cubWeights);

  /**********************************************************************************/
  /*********************************** GET BASIS ************************************/
  /**********************************************************************************/

  *out << "Getting basis" << endl;

  // Define basis
  Basis_HGRAD_HEX_C1_FEM<ST, FieldContainer<ST> > hexHGradBasis;
  int numFieldsG = hexHGradBasis.getCardinality ();
  FieldContainer<ST> HGBValues (numFieldsG, numCubPoints);
  FieldContainer<ST> HGBGrads (numFieldsG, numCubPoints, spaceDim);

  // Evaluate basis values and gradients at cubature points
  hexHGradBasis.getValues (HGBValues, cubPoints, OPERATOR_VALUE);
  hexHGradBasis.getValues (HGBGrads, cubPoints, OPERATOR_GRAD);

  /**********************************************************************************/
  /********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
  /**********************************************************************************/

  *out << "Building Maps" << endl;

  RCP<Teuchos::Time> timerBuildGlobalMaps =
    TimeMonitor::getNewTimer ("Build global Maps and Export");
  Array<int> ownedGIDs;
  RCP<const map_type> globalMapG;
  {
    TimeMonitor timerBuildGlobalMapsL (*timerBuildGlobalMaps);
    // Count owned nodes
    int ownedNodes = 0;
    for (int i = 0; i < numNodes; ++i) {
      if (nodeIsOwned[i]) {
        ++ownedNodes;
      }
    }

    // Build a list of the OWNED global ids...
    // NTS: will need to switch back to long long
    ownedGIDs.resize(ownedNodes);
    int oidx = 0;
    for (int i = 0; i < numNodes; ++i) {
      if (nodeIsOwned[i]) {
        ownedGIDs[oidx] = as<int> (globalNodeIds[i]);
        ++oidx;
      }
    }
    globalMapG = rcp (new map_type (-1, static_cast<int> (ownedGIDs ().size ()),
                                    ownedGIDs ().getRawPtr (), 0, *comm));
  }

  /**********************************************************************************/
  /********************* BUILD MAPS FOR OVERLAPPED SOLUTION *************************/
  /**********************************************************************************/

  Array<GO> overlappedGIDs;
  RCP<const map_type> overlappedMapG;
  RCP<const export_type> exporter;
  {
    // Count owned nodes
    int overlappedNodes = numNodes;

    // Build a list of the OVERLAPPED global ids...
    overlappedGIDs.resize (overlappedNodes);
    for (int i = 0; i < numNodes; ++i) {
      overlappedGIDs[i] = as<int> (globalNodeIds[i]);
    }

    //Generate overlapped Map for nodes.
    overlappedMapG = rcp (new map_type (-1, static_cast<int> (overlappedGIDs ().size ()),
                                        overlappedGIDs ().getRawPtr (), 0, *comm));

    // Build Epetra_Export from overlapped to owned Map.
    exporter = rcp (new export_type (*overlappedMapG, *globalMapG));
  }

  /**********************************************************************************/
  /********************* BUILD GRAPH FOR OVERLAPPED SOLUTION ************************/
  /**********************************************************************************/

  *out << "Building Graph" << endl;

  RCP<sparse_graph_type> overlappedGraph;
  RCP<sparse_graph_type> ownedGraph;
  RCP<Teuchos::Time> timerBuildOverlapGraph =
    TimeMonitor::getNewTimer ("Build graphs for overlapped and owned solutions");
  {
    TimeMonitor timerBuildOverlapGraphL (*timerBuildOverlapGraph);

    // Construct Epetra_CrsGraph objects.
    overlappedGraph = rcp (new sparse_graph_type (Copy, *overlappedMapG, 0));
    ownedGraph = rcp (new sparse_graph_type (Copy, *globalMapG, 0));

    // Define desired workset size and count how many worksets
    // there are on this process's mesh block.
    int desiredWorksetSize = numElems; // change to desired workset size!
    //int desiredWorksetSize = 100;    // change to desired workset size!
    int numWorksets        = numElems/desiredWorksetSize;

    for (int workset = 0; workset < numWorksets; ++workset) {
      // Compute cell numbers where the workset starts and ends
      int worksetSize  = 0;
      int worksetBegin = (workset + 0)*desiredWorksetSize;
      int worksetEnd   = (workset + 1)*desiredWorksetSize;

      // When numElems is not divisible by desiredWorksetSize, the
      // last workset ends at numElems.
      worksetEnd = (worksetEnd <= numElems) ? worksetEnd : numElems;

      // Now we know the actual workset size and can allocate the
      // array for the cell nodes.
      worksetSize = worksetEnd - worksetBegin;

      //"WORKSET CELL" loop: local cell ordinal is relative to numElems
      for (int cell = worksetBegin; cell < worksetEnd; ++cell) {
        // Compute cell ordinal relative to the current workset
        //int worksetCellOrdinal = cell - worksetBegin;

        // "CELL EQUATION" loop for the workset cell: cellRow is
        // relative to the cell DoF numbering
        for (int cellRow = 0; cellRow < numFieldsG; cellRow++){

          int localRow  = elemToNode(cell, cellRow);
          //globalRow for Epetra_CrsGraph
          int globalRowT = as<int> (globalNodeIds[localRow]);

          // "CELL VARIABLE" loop for the workset cell: cellCol is
          // relative to the cell DoF numbering
          for (int cellCol = 0; cellCol < numFieldsG; ++cellCol) {
            int localCol  = elemToNode (cell, cellCol);
            int globalCol = as<int> (globalNodeIds[localCol]);
            //create ArrayView globalCol object for Epetra
            ArrayView<int> globalColAV = arrayView (&globalCol, 1);

            //Update Epetra overlap Graph
            const int errCode =
              overlappedGraph->InsertGlobalIndices (globalRowT,
                                                    as<int> (globalColAV.size ()),
                                                    globalColAV.getRawPtr ());
            // Positive err isn't necessarily an error, but negative
            // err always is.
            TEUCHOS_TEST_FOR_EXCEPTION(errCode < 0, std::runtime_error,
              "Epetra_CrsGraph::InsertGlobalIndices on global row "
              << globalRowT << " on process " << myRank
              << " failed with error code " << errCode << ".");
          }// *** cell col loop ***
        }// *** cell row loop ***
      }// *** workset cell loop **
    }// *** workset loop ***

    // Fill-complete overlapping distribution Graph.
    errCode = overlappedGraph->FillComplete ();
    TEUCHOS_TEST_FOR_EXCEPTION(errCode < 0, std::runtime_error,
      "Epetra_CrsGraph::FillComplete on overlapped graph on process "
      << myRank << " failed with error code " << errCode << ".");

    // Export to owned distribution Graph, and fill-complete the latter.
    errCode = ownedGraph->Export (*overlappedGraph, *exporter, Insert);
    TEUCHOS_TEST_FOR_EXCEPTION(errCode < 0, std::runtime_error,
      "Epetra_CrsGraph::Export from overlapped to owned graph on process "
      << myRank << " failed with error code " << errCode << ".");
    errCode = ownedGraph->FillComplete ();
    TEUCHOS_TEST_FOR_EXCEPTION(errCode < 0, std::runtime_error,
      "Epetra_CrsGraph::FillComplete on owned graph on process "
      << myRank << " failed with error code " << errCode << ".");
  }

  *out << "Constructing stiffness matrix and vectors" << endl;

  //
  // Construct stiffness matrix, right-hand side vector, and exact
  // solution vector.  The linear algebra objects starting with gl_
  // are for the nonoverlapped ("owned") distribution; the ones not
  // starting with gl_ are for the overlapped distribution.  Once
  // we've constructed the overlapped distribution objects, we'll
  // Export to the owned distribution objects.
  //
  // Owned distribution objects: their names start with gl_.
  //
  RCP<sparse_matrix_type> gl_StiffMatrix =
    rcp (new sparse_matrix_type (Copy, *ownedGraph));
  RCP<vector_type> gl_rhsVector = rcp (new vector_type (*globalMapG));
  // We compute the exact solution vector when we make v, and v has
  // the owned distribution, so we only need an owned distribution
  // version of the exact solution vector.
  RCP<vector_type> gl_exactSolVector = rcp (new vector_type (*globalMapG));

  //
  // Overlapped distribution objects: their names don't start with gl_.
  //
  RCP<sparse_matrix_type> StiffMatrix =
    rcp (new sparse_matrix_type (Copy, *overlappedGraph));
  errCode = StiffMatrix->PutScalar (STS::zero ());
  TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
    "Epetra_CrsMatrix::PutScalar on overlapped matrix on process "
    << myRank << " failed with error code " << errCode << ".");
  RCP<vector_type> rhsVector = rcp (new vector_type (*overlappedMapG));
  errCode = rhsVector->PutScalar (STS::zero ());
  TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
    "Epetra_MultiVector::PutScalar on process "
    << myRank << " failed with error code " << errCode << ".");

  /**********************************************************************************/
  /************************** DIRICHLET BC SETUP ************************************/
  /**********************************************************************************/

  *out << "Setting up Dirichlet boundary conditions" << endl;

  // Timer for Dirichlet BC setup
  RCP<Teuchos::Time> timerDirichletBC =
    TimeMonitor::getNewTimer ("Get Dirichlet boundary values: Total Time");
  int numBCNodes = 0;
  RCP<vector_type> v;
  Array<int> BCNodes;
  {
    TimeMonitor timerDirichletBCL(*timerDirichletBC);
    for (int inode = 0; inode < numNodes; ++inode) {
      if (nodeOnBoundary(inode) && nodeIsOwned[inode]) {
        ++numBCNodes;
      }
    }
    // v: Vector for use in applying Dirichlet BCs.
    // v uses the owned (nonoverlapped) distribution.
    v = rcp (new vector_type (*globalMapG, true));
    errCode = v->PutScalar (STS::zero ());
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_MultiVector::PutScalar on process "
      << myRank << " failed with error code " << errCode << ".");

    // Set v to boundary values on Dirichlet nodes.  While we're
    // iterating through all the nodes that this process can see
    // (which includes the nodes which my process owns, as well as
    // nodes on the interface between my process and another process),
    // also set the exact solution values.
    BCNodes.resize (numBCNodes);
    int indbc = 0;
    int iOwned = 0;
    for (int inode = 0; inode < numNodes; ++inode) {
      if (nodeIsOwned[inode]) {
        const ST x = nodeCoord(inode, 0);
        const ST y = nodeCoord(inode, 1);
        const ST z = nodeCoord(inode, 2);

        if (nodeOnBoundary (inode)) {
          BCNodes[indbc]=iOwned;
          ++indbc;
          // We assume that the exact solution is defined on the
          // boundary, and evaluate the exact solution function to get
          // boundary values.
          const ST val = exactSolution (x, y, z);
          errCode = v->ReplaceMyValues (1, &val, &iOwned);
          TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
            "Epetra_Vector::ReplaceMyValues on local row " << iOwned
            << " on process " << myRank << " failed with error code "
            << errCode << ".");
        } // if node inode is on the boundary

        const ST val = exactSolution (x, y, z);
        errCode = gl_exactSolVector->ReplaceMyValues (1, &val, &iOwned);
        TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
          "Epetra_Vector::ReplaceMyValues on local row " << iOwned
          << " on process " << myRank << " failed with error code "
          << errCode << ".");
        ++iOwned;
      } // if node inode is owned by my process
    } // for each node inode that my process can see
  } // Done setting up v (for Dirichlet BCs) and gl_exactSolVector

  /**********************************************************************************/
  /******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
  /**********************************************************************************/

  *out << "Building discretization matrix and right hand side" << endl;

  // Define desired workset size and count how many worksets there are
  // on this processor's mesh block
  int desiredWorksetSize = numElems; // change to desired workset size!
  //int desiredWorksetSize = 100;    // change to desired workset size!
  int numWorksets        = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase
  // workset count by 1
  if (numWorksets*desiredWorksetSize < numElems) {
    numWorksets += 1;
  }

  {
    Teuchos::OSTab tab2 (out);
    *out << "Desired workset size:             " << desiredWorksetSize << endl
         << "Number of worksets (per process): " << numWorksets << endl;
  }

  for (int workset = 0; workset < numWorksets; ++workset) {
    // Compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // When numElems is not divisible by desiredWorksetSize, the last
    // workset ends at numElems.
    worksetEnd = (worksetEnd <= numElems) ? worksetEnd : numElems;

    // Now we know the actual workset size and can allocate the array
    // for the cell nodes.
    worksetSize = worksetEnd - worksetBegin;
    FieldContainer<ST> cellWorkset (worksetSize, numNodesPerElem, spaceDim);

    // Copy coordinates into cell workset
    int cellCounter = 0;
    for (int cell = worksetBegin; cell < worksetEnd; ++cell) {
      for (int node = 0; node < numNodesPerElem; ++node) {
        cellWorkset(cellCounter, node, 0) = nodeCoord( elemToNode(cell, node), 0);
        cellWorkset(cellCounter, node, 1) = nodeCoord( elemToNode(cell, node), 1);
        cellWorkset(cellCounter, node, 2) = nodeCoord( elemToNode(cell, node), 2);
      }
      ++cellCounter;
    }

    /**********************************************************************************/
    /*                                Allocate arrays                                 */
    /**********************************************************************************/

    // Containers for Jacobians, integration measure & cubature points in workset cells
    FieldContainer<ST> worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<ST> worksetJacobInv  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<ST> worksetJacobDet  (worksetSize, numCubPoints);
    FieldContainer<ST> worksetCubWeights(worksetSize, numCubPoints);
    FieldContainer<ST> worksetCubPoints (worksetSize, numCubPoints, cubDim);

    // Containers for basis values transformed to workset cells and
    // them multiplied by cubature weights
    FieldContainer<ST> worksetHGBValues        (worksetSize, numFieldsG, numCubPoints);
    FieldContainer<ST> worksetHGBValuesWeighted(worksetSize, numFieldsG, numCubPoints);
    FieldContainer<ST> worksetHGBGrads         (worksetSize, numFieldsG, numCubPoints, spaceDim);
    FieldContainer<ST> worksetHGBGradsWeighted (worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for diffusive & advective fluxes & non-conservative
    // adv. term and reactive terms
    FieldContainer<ST> worksetDiffusiveFlux(worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for material values and source term. Require
    // user-defined functions
    FieldContainer<ST> worksetMaterialVals (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<ST> worksetSourceTerm   (worksetSize, numCubPoints);

    // Containers for workset contributions to the discretization
    // matrix and the right hand side
    FieldContainer<ST> worksetStiffMatrix (worksetSize, numFieldsG, numFieldsG);
    FieldContainer<ST> worksetRHS         (worksetSize, numFieldsG);

    /**********************************************************************************/
    /*                                Calculate Jacobians                             */
    /**********************************************************************************/

    IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
    IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
    IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

    /**********************************************************************************/
    /*          Cubature Points to Physical Frame and Compute Data                    */
    /**********************************************************************************/

    // Map cubature points to physical frame.
    IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);

    // Evaluate the material tensor A at cubature points.
    evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

    // Evaluate the source term at cubature points.
    evaluateSourceTerm (worksetSourceTerm, worksetCubPoints);

    /**********************************************************************************/
    /*                         Compute Stiffness Matrix                               */
    /**********************************************************************************/

    // Transform basis gradients to physical frame:
    IntrepidFSTools::HGRADtransformGRAD<ST> (worksetHGBGrads,   // DF^{-T}(grad u)
                                             worksetJacobInv,
                                             HGBGrads);
    // Compute integration measure for workset cells:
    IntrepidFSTools::computeCellMeasure<ST> (worksetCubWeights, // Det(DF)*w = J*w
                                             worksetJacobDet,
                                             cubWeights);
    // Multiply transformed (workset) gradients with weighted measure
    IntrepidFSTools::multiplyMeasure<ST> (worksetHGBGradsWeighted, // DF^{-T}(grad u)*J*w
                                          worksetCubWeights,
                                          worksetHGBGrads);
    // Compute the diffusive flux:
    IntrepidFSTools::tensorMultiplyDataField<ST> (worksetDiffusiveFlux, // A*(DF^{-T}(grad u)
                                                  worksetMaterialVals,
                                                  worksetHGBGrads);
    // Integrate to compute workset diffusion contribution to global matrix:
    IntrepidFSTools::integrate<ST> (worksetStiffMatrix, // (DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
                                    worksetHGBGradsWeighted,
                                    worksetDiffusiveFlux,
                                    COMP_BLAS);

    /**********************************************************************************/
    /*                                   Compute RHS                                  */
    /**********************************************************************************/

    // Transform basis values to physical frame:
    IntrepidFSTools::HGRADtransformVALUE<ST> (worksetHGBValues, // clones basis values (u)
                                              HGBValues);
    // Multiply transformed (workset) values with weighted measure
    IntrepidFSTools::multiplyMeasure<ST> (worksetHGBValuesWeighted, // (u)*J*w
                                          worksetCubWeights,
                                          worksetHGBValues);
    // Integrate worksetSourceTerm against weighted basis function set
    IntrepidFSTools::integrate<ST> (worksetRHS, // f.(u)*J*w
                                    worksetSourceTerm,
                                    worksetHGBValuesWeighted,
                                    COMP_BLAS);

    /**********************************************************************************/
    /*                         Assemble into Global Matrix                            */
    /**********************************************************************************/

    RCP<Teuchos::Time> timerAssembleGlobalMatrix =
      TimeMonitor::getNewTimer ("Assemble overlapped global matrix and RHS");
    {
      TimeMonitor timerAssembleGlobalMatrixL (*timerAssembleGlobalMatrix);

      // "WORKSET CELL" loop: local cell ordinal is relative to numElems
      for (int cell = worksetBegin; cell < worksetEnd; ++cell) {

        // Compute cell ordinal relative to the current workset
        const int worksetCellOrdinal = cell - worksetBegin;

        // "CELL EQUATION" loop for the workset cell: cellRow is
        // relative to the cell DoF numbering.
        for (int cellRow = 0; cellRow < numFieldsG; ++cellRow) {
          int localRow  = elemToNode (cell, cellRow);
          int globalRow = as<int> (globalNodeIds[localRow]);
          ST sourceTermContribution = worksetRHS (worksetCellOrdinal, cellRow);
          ArrayView<ST> sourceTermContributionAV =
            arrayView (&sourceTermContribution, 1);

	  {
	    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Assembly: Element, RHS", elem_rhs);
	    errCode = rhsVector->SumIntoGlobalValues (1, sourceTermContributionAV.getRawPtr (), &globalRow);
	  }
          TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
            "Epetra_Vector::SumIntoGlobalValues on global row " << globalRow
            << " on process " << myRank << " failed with error code "
            << errCode << ".");

          // "CELL VARIABLE" loop for the workset cell: cellCol is
          // relative to the cell DoF numbering.
          for (int cellCol = 0; cellCol < numFieldsG; cellCol++){
            const int localCol  = elemToNode(cell, cellCol);
            int globalCol = as<int> (globalNodeIds[localCol]);
            ArrayView<int> globalColAV = arrayView<int> (&globalCol, 1);
            ST operatorMatrixContribution =
              worksetStiffMatrix (worksetCellOrdinal, cellRow, cellCol);
            ArrayView<ST> operatorMatrixContributionAV =
              arrayView<ST> (&operatorMatrixContribution, 1);
	    {
	      TEUCHOS_FUNC_TIME_MONITOR_DIFF("Assembly: Element, Matrix", 
					     elem_matrix);
	      errCode = StiffMatrix->SumIntoGlobalValues (
		globalRow,
		as<int> (operatorMatrixContributionAV.size ()),
		operatorMatrixContributionAV.getRawPtr (),
		globalColAV.getRawPtr ());
	    }
            TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
              "Epetra_CrsMatrix::SumIntoGlobalValues on global row "
              << globalRow << " on process " << myRank << " failed with "
              "error code " << errCode << ".");

          }// *** cell col loop ***
        }// *** cell row loop ***
      }// *** workset cell loop **
    } // *** stop timer ***
  }// *** workset loop ***

  //////////////////////////////////////////////////////////////////////////////
  // Export sparse matrix and right-hand side from overlapping row Map
  // to owned (nonoverlapping) row Map.
  //////////////////////////////////////////////////////////////////////////////

  *out << "Exporting matrix and right-hand side from overlapped to owned Map" << endl;

  RCP<Teuchos::Time> timerAssembMultProc =
    TimeMonitor::getNewTimer ("Export from overlapped to owned");
  {
    TimeMonitor timerAssembMultProcL (*timerAssembMultProc);
    errCode = gl_StiffMatrix->PutScalar (STS::zero ());
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_CrsMatrix::PutScalar on process " << myRank << " failed with "
      "error code " << errCode << ".");

    errCode = gl_StiffMatrix->Export (*StiffMatrix, *exporter, InsertAdd);
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_CrsMatrix::Export on process " << myRank << " failed with "
      "error code " << errCode << ".");

    // If target of export has static graph, no need to do
    // PutScalar(0.0); export will clobber values.
    errCode = gl_StiffMatrix->FillComplete ();
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_CrsMatrix::FillComplete on process " << myRank
       << "failed with error code " << errCode << ".");

    errCode = gl_rhsVector->PutScalar (STS::zero ());
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_Vector::PutScalar on process " << myRank << " failed with "
      "error code " << errCode << ".");

    errCode = gl_rhsVector->Export (*rhsVector, *exporter, Add);
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_Vector::Export on process " << myRank << " failed with "
      "error code " << errCode << ".");
  }

  //////////////////////////////////////////////////////////////////////////////
  // Adjust matrix for boundary conditions
  //////////////////////////////////////////////////////////////////////////////

  *out << "Adjusting matrix and right-hand side for BCs" << endl;

  RCP<Teuchos::Time> timerAdjustMatrixBC =
    TimeMonitor::getNewTimer ("Adjust owned matrix and RHS for BCs");
  {
    TimeMonitor timerAdjustMatrixBCL (*timerAdjustMatrixBC);

    // Apply owned stiffness matrix to v: rhs := A*v
    RCP<multivector_type> rhsDir =
      rcp (new multivector_type (*globalMapG, 1, true));
    errCode = gl_StiffMatrix->Apply (*v.getConst (), *rhsDir);
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_Operator::Apply with global stiffness matrix on process "
      << myRank << " failed with error code " << errCode << ".");

    // Update right-hand side
    errCode = gl_rhsVector->Update (as<ST> (-STS::one ()), *rhsDir, STS::one ());
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_MultiVector::Apply with global right-hand side vector on "
      "process " << myRank << " failed with error code " << errCode << ".");

    // Epetra_(Multi)Vector's operator[] returns a raw pointer to the
    // (multi)vector's local data.  Wrap it in a nonowning ArrayRCP.
    ArrayRCP<const ST> vArrayRCP;
    {
      const ST* v_data = &((*v)[0]);
      const ArrayRCP<const ST>::size_type lowerOffset = 0;
      const ArrayRCP<const ST>::size_type size = v->MyLength ();
      const bool ownsMem = false;
      vArrayRCP = arcp<const ST> (v_data, lowerOffset, size, ownsMem);
    //arcp_const_cast<const ST> (arcp<ST> (v_data, lowerOffset, size, ownsMem));
    }

    // Adjust rhs due to Dirichlet boundary conditions.
    for (int inode = 0; inode < numNodes; ++inode) {
      if (nodeIsOwned[inode]) {
        if (nodeOnBoundary (inode)) {
          // Get global node ID
          const GO gni = as<GO> (globalNodeIds[inode]);
          const LO lidT = globalMapG->LID (gni);
          ST v_valT = vArrayRCP[lidT];
          errCode = gl_rhsVector->ReplaceGlobalValues (1, &v_valT, &gni);
          TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
            "Epetra_Vector::ReplaceGlobalValues with global right-hand side "
            "vector on process " << myRank << " failed with error code "
            << errCode << ".");
        }
      }
    }

    // Zero out rows and columns of stiffness matrix corresponding to
    // Dirichlet edges and add one to diagonal.  The following is the
    // "pure Epetra" analog of Apply_OAZToMatrix().
    //
    // Epetra_CrsMatrix, unlike Tpetra::CrsMatrix, doesn't require a
    // "resumeFill" call before modifying values.
    //gl_StiffMatrix->resumeFill ();

    // Find the local column numbers to nuke
    RCP<const map_type> ColMap = rcpFromRef (gl_StiffMatrix->ColMap ());
    RCP<const map_type> globalMap =
      rcp (new map_type (gl_StiffMatrix->NumGlobalCols (), 0, *comm));

    // Create the exporter from this process' column Map to the global
    // 1-1 column map. (???)
    RCP<const export_type> bdyExporter =
      rcp (new export_type (*ColMap, *globalMap));
    // Create a vector of global column indices to which we will export
    RCP<Epetra_IntVector> globColsToZeroT =
      rcp (new Epetra_IntVector (*globalMap));
    // Create a vector of local column indices from which we will export
    RCP<Epetra_IntVector> myColsToZeroT = rcp (new Epetra_IntVector (*ColMap));
    errCode = myColsToZeroT->PutValue (0);
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_IntVector::PutValue failed with error code " << errCode << ".");

    // Flag (set to 1) all local columns corresponding to the local
    // rows specified.
    for (int i = 0; i < numBCNodes; ++i) {
      const GO globalRow = gl_StiffMatrix->RowMap ().GID (BCNodes[i]);
      const LO localCol = gl_StiffMatrix->ColMap ().LID (globalRow);
      // Epetra_IntVector has an operator[] which takes an LID.
      (*myColsToZeroT)[localCol] = 1;
    }

    // Export to the global column map.
    errCode = globColsToZeroT->Export (*myColsToZeroT, *bdyExporter, Add);
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_IntVector::Export failed with error code " << errCode << ".");

    // Import from the global column map to the local column map.
    errCode = myColsToZeroT->Import (*globColsToZeroT, *bdyExporter, Insert);
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
      "Epetra_IntVector::Import (with exporter) failed with error code "
      << errCode << ".");

    Array<ST> values;
    Array<int> indices;
    // Epetra_IntVector's Values() method returns a raw pointer to the
    // (multi)vector's local data.  Wrap it in a nonowning ArrayRCP.
    ArrayRCP<const int> myColsToZeroArrayRCP =
      arcp_const_cast<const int> (arcp<int> (myColsToZeroT->Values (), 0,
                                             v->MyLength (), false));
    //size_t NumEntries = 0;
    int NumEntries = 0;

    // Zero the columns corresponding to Dirichlet BCs.
    for (LO i = 0; i < gl_StiffMatrix->NumMyRows (); ++i) {
      errCode = gl_StiffMatrix->NumMyRowEntries (i, NumEntries);
      TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
        "Epetra_RowMatrix::NumMyRowEntries on process " << myRank
        << " and local row " << i << "failed with error code "
        << errCode << ".");

      values.resize (NumEntries);
      indices.resize (NumEntries);
      errCode = gl_StiffMatrix->ExtractMyRowCopy (i, NumEntries, NumEntries,
                                              values ().getRawPtr (),
                                              indices ().getRawPtr ());
      TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
        "Epetra_RowMatrix::ExtractMyRowCopy on process " << myRank
        << " and local row " << i << "failed with error code "
        << errCode << ".");

      for (int j = 0; j < NumEntries; ++j) {
        if (myColsToZeroArrayRCP[indices[j]] == 1)
          values[j] = STS::zero ();
      }
      errCode = gl_StiffMatrix->ReplaceMyValues (i, as<int> (values ().size ()),
                                             values ().getRawPtr (),
                                             indices ().getRawPtr ());
      TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
        "Epetra_CrsMatrix::ReplaceMyValues on process " << myRank
        << " and local row " << i << "failed with error code "
        << errCode << ".");
    } // for each (local) row of the global stiffness matrix

    // Zero the rows and add ones to diagonal.
    for (int i = 0; i < numBCNodes; ++i) {
      errCode = gl_StiffMatrix->NumMyRowEntries (BCNodes[i], NumEntries);
      TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
        "Epetra_RowMatrix::NumMyRowEntries on process " << myRank
        << " and local row " << BCNodes[i] << "failed with error code "
        << errCode << ".");

      indices.resize (NumEntries);
      values.resize (NumEntries);
      errCode = gl_StiffMatrix->ExtractMyRowCopy (BCNodes[i], NumEntries,
                                              NumEntries,
                                              values ().getRawPtr (),
                                              indices ().getRawPtr ());
      TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
        "Epetra_RowMatrix::ExtractMyRowCopy on process " << myRank
        << " and local row " << BCNodes[i] << "failed with error code "
        << errCode << ".");

      const GO globalRow = gl_StiffMatrix->RowMap ().GID (BCNodes[i]);
      const LO localCol = gl_StiffMatrix->ColMap ().LID (globalRow);
      for (int j = 0; j < NumEntries; ++j) {
        values[j] = STS::zero ();
        if (indices[j] == localCol) {
          values[j] = STS::one ();
        }
      } // for each entry in the current row
      errCode = gl_StiffMatrix->ReplaceMyValues (BCNodes[i],
                                             as<int> (values ().size ()),
                                             values ().getRawPtr (),
                                             indices ().getRawPtr ());
      TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
        "Epetra_CrsMatrix::ReplaceMyValues on process " << myRank
        << " and local row " << BCNodes[i] << "failed with error code "
        << errCode << ".");

    } // for each BC node
  }

  //*out << "Calling fillComplete() on owned-Map stiffness matrix" << endl;

  // We're done modifying the owned stiffness matrix.  Epetra allows
  // modifications to values (not structure) after calling
  // FillComplete(), so we didn't have to call "resumeFill" above, and
  // we don't need to call FillComplete() again here (in fact, it's
  // not allowed to call FillComplete() more than once).
  // gl_StiffMatrix->FillComplete ();

  //
  // We're done with assembly, so we can delete the mesh.
  //
  delete [] nodeIsOwned;
  nodeIsOwned = NULL;
  Delete_Pamgen_Mesh ();

  // Create vector to store approximate solution, and set initial guess.
  RCP<vector_type> gl_approxSolVector = rcp (new vector_type (*globalMapG));
  errCode = gl_approxSolVector->PutScalar (STS::zero ());
  TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
    "Epetra_Vector::PutScalar on process " << myRank
    << "failed with error code " << errCode << ".");

  // Store the output pointers.
  A = gl_StiffMatrix;
  B = gl_rhsVector;
  X_exact = gl_exactSolVector;
  X = gl_approxSolVector;
}

std::vector<Teuchos::ScalarTraits<ST>::magnitudeType>
exactResidualNorm (const Teuchos::RCP<const sparse_matrix_type>& A,
                   const Teuchos::RCP<const vector_type>& B,
                   const Teuchos::RCP<const vector_type>& X_exact)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef STS::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;

  RCP<vector_type> R = rcp (new vector_type (B->Map ()));
  // R := A*X_exact
  int errCode = A->Apply (*X_exact, *R);
  TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error, "exactResidualNorm: "
    "Epetra_Operator::Apply failed with error code " << errCode << ".");

  // R := -1*R + 1*B == B - A*X_exact
  errCode = R->Update (-STS::one(), *B, STS::one());
  TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error, "exactResidualNorm: "
    "Epetra_MultiVector::Update failed with error code " << errCode << ".");

  MT R_norm = STM::zero();
  MT B_norm = STM::zero();
  errCode = R->Norm2 (&R_norm);
  TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error, "exactResidualNorm: "
    "Epetra_MultiVector::Norm2 failed with error code " << errCode << ".");
  errCode = B->Norm2 (&B_norm);
  TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error, "exactResidualNorm: "
    "Epetra_MultiVector::Norm2 failed with error code " << errCode << ".");

  std::vector<MT> results (3);
  results[0] = R_norm;
  results[1] = B_norm;
  results[2] = A->NormFrobenius ();
  return results;
}

/**********************************************************************************/
/************ USER DEFINED FUNCTIONS FOR EXACT SOLUTION ***************************/
/**********************************************************************************/

template<typename Scalar>
const Scalar
exactSolution (const Scalar& x, const Scalar& y, const Scalar& z)
{
  // Patch test: tri-linear function is in the FE space and should be recovered
  return 1. + x + y + z + x*y + x*z + y*z + x*y*z;

  // Analytic solution with homogeneous Dirichlet boundary data
  // return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);

  // Analytic solution with inhomogeneous Dirichlet boundary data
  // return exp(x + y + z)/(1. + x*y + y*z + x*y*z);
}


template<typename Scalar>
void
materialTensor (Scalar material[][3],
                const Scalar& x,
                const Scalar& y,
                const Scalar& z)
{
  typedef Teuchos::ScalarTraits<Scalar> STS;

  const bool illConditioned = false;
  if (illConditioned) {
    const Scalar zero = STS::zero ();
    const Scalar one = STS::one ();
    const Scalar two = one + one;
    const Scalar four = two + two;
    const Scalar eight = four + four;
    (void) four;
    (void) eight;

    if (false) {
      material[0][0] = one;
      material[0][1] = one;
      material[0][2] = zero;

      material[1][0] = one;
      material[1][1] = one;
      material[1][2] = one;

      material[2][0] = one;
      material[2][1] = one;
      material[2][2] = zero;
    }
    else {
      material[0][0] = one;
      material[0][1] = one - one / two;
      material[0][2] = zero;

      material[1][0] = one;
      material[1][1] = one - one / two;
      material[1][2] = one;

      material[2][0] = one;
      material[2][1] = one - one / two;
      material[2][2] = zero;
    }
  }
  else {
    material[0][0] = STS::one();
    material[0][1] = STS::zero();
    material[0][2] = STS::zero();

    material[1][0] = STS::zero();
    material[1][1] = STS::one();
    material[1][2] = STS::zero();

    material[2][0] = STS::zero();
    material[2][1] = STS::zero();
    material[2][2] = STS::one();
  }
}

/**********************************************************************************/
/************** AUXILIARY FUNCTIONS FROM EXACT SOLUTION ***************************/
/**********************************************************************************/

/// \brief Compute gradient of exact solution at (x,y,z).
///
/// \param gradExact [out] The gradient at (x,y,z).
/// \param x [in]
/// \param y [in]
/// \param z [in]
template<typename Scalar>
void
exactSolutionGrad (Scalar gradExact[3],
                   const Scalar& x,
                   const Scalar& y,
                   const Scalar& z)
{
  // To enable taking derivatives of the gradient (that is, second
  // derivatives of the exact solution), we need 2 levels of fad
  // types.
  Sacado::Fad::SFad<Scalar,3> fad_x = x;
  Sacado::Fad::SFad<Scalar,3> fad_y = y;
  Sacado::Fad::SFad<Scalar,3> fad_z = z;
  Sacado::Fad::SFad<Scalar,3> u;

  // Indicate the independent variables
  fad_x.diff(0,3);
  fad_y.diff(1,3);
  fad_z.diff(2,3);

  u = exactSolution(fad_x, fad_y, fad_z);

  gradExact[0] = u.dx(0);
  gradExact[1] = u.dx(1);
  gradExact[2] = u.dx(2);
}

/// \brief Compute the source term (right-hand side) at (x,y,z).
///
/// We compute the source term from the exact solution, its gradient,
/// and the material tensor at (x,y,z).
template<typename Scalar>
const Scalar
sourceTerm (Scalar& x, Scalar& y, Scalar& z)
{
  Scalar u;
  Scalar grad_u[3];
  Scalar flux[3];
  Scalar material[3][3];
  Scalar f = 0.;

  // Indicate the independent variables
  x.diff(0,3);
  y.diff(1,3);
  z.diff(2,3);

  // Get exact solution and its gradient
  u = exactSolution (x, y, z);
  exactSolutionGrad (grad_u, x, y, z);

  // Get material tensor
  materialTensor<Scalar> (material, x, y, z);

  // Compute total flux = (A.grad u)
  for (int i = 0; i < 3; ++i) {
    // Add diffusive flux
    for (int j = 0; j < 3; ++j) {
      flux[i] += material[i][j]*grad_u[j];
    }
  }

  // Compute source term (right hand side): f = -div(A.grad u)
  f = -(flux[0].dx(0) + flux[1].dx(1) + flux[2].dx(2));

  return f;
}

/**********************************************************************************/
/*************************** EVALUATION METHODS ***********************************/
/**********************************************************************************/

//! Compute the material tensor over a workset.
template<class ArrayOut, class ArrayIn>
void
evaluateMaterialTensor (ArrayOut&      matTensorValues,
                        const ArrayIn& evaluationPoints)
{
  typedef typename ArrayOut::scalar_type scalar_type;

  const int numWorksetCells  = evaluationPoints.dimension(0);
  const int numPoints        = evaluationPoints.dimension(1);
  const int spaceDim         = evaluationPoints.dimension(2);

  scalar_type material[3][3];

  for (int cell = 0; cell < numWorksetCells; ++cell) {
    for (int pt = 0; pt < numPoints; ++pt) {
      scalar_type x = evaluationPoints(cell, pt, 0);
      scalar_type y = evaluationPoints(cell, pt, 1);
      scalar_type z = evaluationPoints(cell, pt, 2);

      materialTensor<scalar_type> (material, x, y, z);

      for (int row = 0; row < spaceDim; ++row) {
        for(int col = 0; col < spaceDim; ++col) {
          matTensorValues(cell, pt, row, col) = material[row][col];
        }
      }
    }
  }
}

/************ Source Term (RHS) ****************/
template<class ArrayOut, class ArrayIn>
void
evaluateSourceTerm (ArrayOut&      sourceTermValues,
                    const ArrayIn& evaluationPoints)
{
  typedef typename ArrayOut::scalar_type scalar_type;

  const int numWorksetCells  = evaluationPoints.dimension(0);
  const int numPoints = evaluationPoints.dimension(1);

  for (int cell = 0; cell < numWorksetCells; ++cell) {
    for (int pt = 0; pt < numPoints; ++pt) {
      Sacado::Fad::SFad<scalar_type,3> x = evaluationPoints(cell, pt, 0);
      Sacado::Fad::SFad<scalar_type,3> y = evaluationPoints(cell, pt, 1);
      Sacado::Fad::SFad<scalar_type,3> z = evaluationPoints(cell, pt, 2);

      sourceTermValues(cell, pt) =
        sourceTerm<Sacado::Fad::SFad<scalar_type,3> > (x, y, z).val ();
    }
  }
}

/************ Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void
evaluateExactSolution (ArrayOut&      exactSolutionValues,
                       const ArrayIn& evaluationPoints)
{
  typedef typename ArrayOut::scalar_type scalar_type;

  const int numWorksetCells  = evaluationPoints.dimension(0);
  const int numPoints = evaluationPoints.dimension(1);

  for (int cell = 0; cell < numWorksetCells; ++cell) {
    for (int pt = 0; pt < numPoints; ++pt) {
      const scalar_type x = evaluationPoints(cell, pt, 0);
      const scalar_type y = evaluationPoints(cell, pt, 1);
      const scalar_type z = evaluationPoints(cell, pt, 2);
      exactSolutionValues(cell, pt) = exactSolution<scalar_type> (x, y, z);
    }
  }
}

/************ Grad of Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void
evaluateExactSolutionGrad (ArrayOut&      exactSolutionGradValues,
                           const ArrayIn& evaluationPoints)
{
  typedef typename ArrayOut::scalar_type scalar_type;

  const int numWorksetCells  = evaluationPoints.dimension(0);
  const int numPoints = evaluationPoints.dimension(1);
  const int spaceDim  = evaluationPoints.dimension(2);

  scalar_type gradient[3];

  for (int cell = 0; cell < numWorksetCells; ++cell) {
    for (int pt = 0; pt < numPoints; ++pt) {
      const scalar_type x = evaluationPoints(cell, pt, 0);
      const scalar_type y = evaluationPoints(cell, pt, 1);
      const scalar_type z = evaluationPoints(cell, pt, 2);

      exactSolutionGrad<scalar_type> (gradient, x, y, z);

      for (int row = 0; row < spaceDim; ++row) {
        exactSolutionGradValues(cell, pt, row) = gradient[row];
      }
    }
  }
}

} // namespace EpetraIntrepidPoissonExample
} // namespace TrilinosCouplings
