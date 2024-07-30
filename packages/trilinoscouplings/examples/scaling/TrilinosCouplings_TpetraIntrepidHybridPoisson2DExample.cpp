// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
//#include <Teuchos_TimeMonitor.hpp>

// Shards includes
#include <Shards_CellTopology.hpp>

// Pamgen includes
#include <create_inline_mesh.h>
#include <pamgen_im_exodusII_l.h>
#include <pamgen_im_ne_nemesisI_l.h>
#include <pamgen_extras.h>

#ifdef HAVE_INTREPID_KOKKOS
#include "Sacado.hpp"
#else
// Sacado includes
#include <Sacado_No_Kokkos.hpp>
#endif

// My includes
#include "TrilinosCouplings_TpetraIntrepidHybridPoisson2DExample.hpp"
#include "TrilinosCouplings_Pamgen_Utils.hpp"
#include "TrilinosCouplings_IntrepidPoissonExampleHelpers.hpp"

namespace TrilinosCouplings {
namespace TpetraIntrepidPoissonExample {

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
                            Teuchos::Array<Teuchos::Array<ST> >& coordArray,
                            Teuchos::Array<LO>& lNodesPerDim,
                            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const std::string& meshInput,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose,
                            const bool debug)
{
  using Teuchos::RCP;
  using Teuchos::rcp_implicit_cast;

  RCP<vector_type> b, x_exact, x;
  makeMatrixAndRightHandSide (A, b, x_exact, x, coordArray, lNodesPerDim, comm,
                              meshInput, out, err, verbose, debug);

  B = rcp_implicit_cast<multivector_type> (b);
  X_exact = rcp_implicit_cast<multivector_type> (x_exact);
  X = rcp_implicit_cast<multivector_type> (x);
}

void
makeMatrixAndRightHandSide (Teuchos::RCP<sparse_matrix_type>& A,
                            Teuchos::RCP<vector_type>& B,
                            Teuchos::RCP<vector_type>& X_exact,
                            Teuchos::RCP<vector_type>& X,
                            Teuchos::Array<Teuchos::Array<ST> >& coordArray,
                            Teuchos::Array<LO>& lNodesPerDim,
                            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const std::string& meshInput,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose,
                            const bool debug)
{
  using namespace Intrepid;
  using Tpetra::global_size_t;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::arrayView;
  using Teuchos::as;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using std::endl;
  //typedef Teuchos::ArrayView<LO>::size_type size_type; // unused
  typedef Teuchos::ScalarTraits<ST> STS;

  (void) verbose;
  (void) debug;

  Tpetra::global_size_t INVALID_GO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  //
  // mfh 19 Apr 2012: If you want to change the template parameters of
  // these typedefs, modify the typedefs (ST, LO, GO, Node) near the
  // top of this file.
  //
  typedef Tpetra::Map<LO, GO, Node>         map_type;
  typedef Tpetra::Export<LO, GO, Node>      export_type;
  //typedef Tpetra::Import<LO, GO, Node>      import_type; // unused
  typedef Tpetra::CrsGraph<LO, GO, Node>    sparse_graph_type;

  // Number of independent variables fixed at 3
  //typedef Sacado::Fad::SFad<ST, 3>     Fad3; // unused
  typedef Intrepid::FunctionSpaceTools IntrepidFSTools;
  //typedef Intrepid::RealSpaceTools<ST> IntrepidRSTools; // unused
  typedef Intrepid::CellTools<ST>      IntrepidCTools;

  const int numProcs = comm->getSize ();
  const int myRank = comm->getRank ();

  *out << "makeMatrixAndRightHandSide:" << endl;
  Teuchos::OSTab tab (out);

  /**********************************************************************************/
  /***************************** GET CELL TOPOLOGY **********************************/
  /**********************************************************************************/

  *out << "Getting cell topology" << endl;

  // Get cell topology for base hexahedron
  shards::CellTopology cellType(shards::getCellTopologyData<shards::Quadrilateral<4> >() );

  // Get dimensions
  int numNodesPerElem = cellType.getNodeCount();
  int spaceDim = cellType.getDimension();
  int dim = 2;

  /**********************************************************************************/
  /******************************* GENERATE MESH ************************************/
  /**********************************************************************************/

  *out << "Generating mesh (tpetra)" << endl;

  int error = 0; // Number of errors in generating the mesh

  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  long long cr_result = Create_Pamgen_Mesh (meshInput.c_str (), dim, myRank, numProcs, maxInt);
  TrilinosCouplings::pamgen_error_check(*out,cr_result);

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

  long long local_num_ijk[3];
  im_ne_get_local_num_ijk_l (id, local_num_ijk);
  for(int dim = 0; dim < 3; ++dim) {
    lNodesPerDim[dim] = local_num_ijk[dim];
  }
  // Ensure we have the correct nodes per dim for a 2D problem
  lNodesPerDim[2] = 1;

  Teuchos::Array<LO> nodesPerDim(3);
  long long num_ijk[3];
  im_ne_get_num_ijk_l (id, num_ijk);
  for(int dim = 0; dim < 3; ++dim) {
    nodesPerDim[dim] = num_ijk[dim];
  }

  // Print mesh information
  {
    Teuchos::OSTab tab2 (out);
    *out << "Global number of elements:                     "
         << numElemsGlobal << endl
         << "Global number of nodes (incl. boundary nodes): "
         << numNodesGlobal << endl
         << "lNodesPerDim: " << lNodesPerDim << endl
         << "nodesPerDim: " << nodesPerDim << endl;
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
        int sideNode0 = cellType.getNodeMap(1,sideSetSideList[j]-1,0);
        int sideNode1 = cellType.getNodeMap(1,sideSetSideList[j]-1,1);

        nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode0))=1;
        nodeOnBoundary(elemToNode(sideSetElemList[j]-1,sideNode1))=1;
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
  Basis_HGRAD_QUAD_C1_FEM<ST,FieldContainer<ST> > hexHGradBasis;

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

  Array<GO> ownedGIDs;
  RCP<const map_type> globalMapG;
  Teuchos::Array<ST> coordXArray;
  Teuchos::Array<ST> coordYArray;
  Teuchos::Array<ST> coordZArray;
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Build global maps", build_maps);

    // Count owned nodes
    int ownedNodes = 0;
    for (int i = 0; i < numNodes; ++i) {
      if (nodeIsOwned[i]) {
        ++ownedNodes;
      }
    }

    // Build a list of the OWNED global ids...
    // NTS: will need to switch back to long long
    // Get owned coordinates
    ownedGIDs.resize(ownedNodes);
    coordXArray.resize(ownedNodes);
    coordYArray.resize(ownedNodes);
    coordZArray.resize(ownedNodes);

    // It is expected that some nodes are duplicated locally
    // so the values extracted with lNodePerDim are possibly
    // wrong... Check that node i=0 and node i=nx are owned.
    Teuchos::Array<LO> lNodesPerDimOld(3);
    for(int k = 0; k < 3; ++k) {lNodesPerDimOld[k] = lNodesPerDim[k];}
//    if(lNodesPerDimOld[0] > 2 && lNodesPerDimOld[1] > 2 && lNodesPerDimOld[2] > 2) {
//    if(lNodesPerDimOld[0] >= 2 && lNodesPerDimOld[1] >= 2 && lNodesPerDimOld[2] >= 2) {
      LO idx = 0;
      Teuchos::Array<LO> factors(3);
      for(int k = 0; k < 3; ++k) {factors[k] = lNodesPerDimOld[k] / 2;}

      // check ilo face
      idx = factors[2]*lNodesPerDimOld[1]*lNodesPerDimOld[0]
        + factors[1]*lNodesPerDimOld[0];
      if(!nodeIsOwned[idx]) {--lNodesPerDim[0];}

      // check ihi face
      idx = factors[2]*lNodesPerDimOld[1]*lNodesPerDimOld[0]
        + factors[1]*lNodesPerDimOld[0] + lNodesPerDimOld[0] - 1;
      if(!nodeIsOwned[idx]) {--lNodesPerDim[0];}

      // check jlo face
      idx = factors[2]*lNodesPerDimOld[1]*lNodesPerDimOld[0]
        + 0*lNodesPerDimOld[0] + factors[0];
      if(!nodeIsOwned[idx]) {--lNodesPerDim[1];}

      // check jhi face
      idx = factors[2]*lNodesPerDimOld[1]*lNodesPerDimOld[0]
        + (lNodesPerDimOld[1] - 1)*lNodesPerDimOld[0] + factors[0];
      if(!nodeIsOwned[idx]) {--lNodesPerDim[1];}

      // check klo face
      idx = 0*lNodesPerDimOld[1]*lNodesPerDimOld[0]
        + factors[1]*lNodesPerDimOld[0] + factors[0];
      if(!nodeIsOwned[idx]) {--lNodesPerDim[2];}

      // check khi face
      idx = (lNodesPerDimOld[2] - 1)*lNodesPerDimOld[1]*lNodesPerDimOld[0]
        + factors[1]*lNodesPerDimOld[0] + factors[0];
      if(!nodeIsOwned[idx]) {--lNodesPerDim[2];}
//    }

    int oidx = 0;
    for (int i = 0; i < numNodes; ++i) {
      if (nodeIsOwned[i]) {
        ownedGIDs[oidx] = as<GO> (globalNodeIds[i]);

        coordXArray[oidx] = nodeCoord(i,0);
        coordYArray[oidx] = nodeCoord(i,1);
        coordZArray[oidx] = nodeCoord(i,2);

        ++oidx;
      }
    }
    coordArray[0] = coordXArray;
    coordArray[1] = coordYArray;
    coordArray[2] = coordZArray;
    globalMapG = rcp (new map_type (INVALID_GO, ownedGIDs (), 0, comm));
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
    overlappedMapG = rcp (new map_type (INVALID_GO, overlappedGIDs (), 0, comm));

    // Build Tpetra Export from overlapped to owned Map.
    exporter = rcp (new export_type (overlappedMapG, globalMapG));
  }

  /**********************************************************************************/
  /********************* BUILD GRAPH FOR OVERLAPPED SOLUTION ************************/
  /**********************************************************************************/

  *out << "Building Graph" << endl;

  RCP<sparse_graph_type> overlappedGraph;
  RCP<sparse_graph_type> ownedGraph;
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Build matrix graph: 0-Total", graph_total);

    // Construct Tpetra::CrsGraph objects.
    overlappedGraph = rcp (new sparse_graph_type (overlappedMapG, 0));
    ownedGraph = rcp (new sparse_graph_type (globalMapG, 0));

    // Define desired workset size and count how many worksets
    // there are on this process's mesh block.
    //int desiredWorksetSize = numElems; // change to desired workset size!
    int desiredWorksetSize = 100;    // change to desired workset size!
    int numWorksets        = (numElems+desiredWorksetSize-1)/desiredWorksetSize;

    for (int workset = 0; workset < numWorksets; ++workset) {
      // Compute cell numbers where the workset starts and ends
      int worksetBegin = (workset + 0)*desiredWorksetSize;
      int worksetEnd   = (workset + 1)*desiredWorksetSize;

      // When numElems is not divisible by desiredWorksetSize, the
      // last workset ends at numElems.
      worksetEnd = (worksetEnd <= numElems) ? worksetEnd : numElems;

      // Now we know the actual workset size and can allocate the
      // array for the cell nodes.
      //int worksetSize = worksetEnd - worksetBegin;

      //"WORKSET CELL" loop: local cell ordinal is relative to numElems
      for (int cell = worksetBegin; cell < worksetEnd; ++cell) {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF(
          "Build matrix graph: 1-Insert global indices", graph_insert);

        // Compute cell ordinal relative to the current workset
        //int worksetCellOrdinal = cell - worksetBegin;

        // "CELL EQUATION" loop for the workset cell: cellRow is
        // relative to the cell DoF numbering
        for (int cellRow = 0; cellRow < numFieldsG; cellRow++){

          int localRow  = elemToNode(cell, cellRow);
          //globalRow for Tpetra Graph
          global_size_t globalRowT = as<global_size_t> (globalNodeIds[localRow]);

          // "CELL VARIABLE" loop for the workset cell: cellCol is
          // relative to the cell DoF numbering
          for (int cellCol = 0; cellCol < numFieldsG; ++cellCol) {
            int localCol  = elemToNode (cell, cellCol);
            GO globalCol = as<GO> (globalNodeIds[localCol]);
            //create ArrayView globalCol object for Tpetra
            ArrayView<GO> globalColAV = arrayView (&globalCol, 1);

            //Update Tpetra overlap Graph
            overlappedGraph->insertGlobalIndices (globalRowT, globalColAV);
          }// *** cell col loop ***
        }// *** cell row loop ***
      }// *** workset cell loop **
    }// *** workset loop ***

    // Fill-complete overlapping distribution Graph.
    {
      TEUCHOS_FUNC_TIME_MONITOR_DIFF(
        "Build matrix graph: 2-Overlapped fill complete",
        overlapped_fill_complete);
      overlappedGraph->fillComplete ();
    }

    // Export to owned distribution Graph, and fill-complete the latter.
    {
      TEUCHOS_FUNC_TIME_MONITOR_DIFF(
        "Build matrix graph: 3-Owned graph export", graph_export);
      ownedGraph->doExport (*overlappedGraph, *exporter, Tpetra::INSERT);
    }
    {
      TEUCHOS_FUNC_TIME_MONITOR_DIFF(
        "Build matrix graph: 4-Owned fill complete", owned_fill_complete);
      ownedGraph->fillComplete ();
    }
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
    rcp (new sparse_matrix_type (ownedGraph.getConst ()));
  RCP<vector_type> gl_rhsVector = rcp (new vector_type (globalMapG));
  // We compute the exact solution vector when we make v, and v has
  // the owned distribution, so we only need an owned distribution
  // version of the exact solution vector.
  RCP<vector_type> gl_exactSolVector = rcp (new vector_type (globalMapG));

  //
  // Overlapped distribution objects: their names don't start with gl_.
  //
  RCP<sparse_matrix_type> StiffMatrix =
    rcp (new sparse_matrix_type (overlappedGraph.getConst ()));
  StiffMatrix->setAllToScalar (STS::zero ());
  RCP<vector_type> rhsVector = rcp (new vector_type (overlappedMapG));
  rhsVector->putScalar (STS::zero ());

  /**********************************************************************************/
  /************************** DIRICHLET BC SETUP ************************************/
  /**********************************************************************************/

  *out << "Setting up Dirichlet boundary conditions" << endl;

  int numBCNodes = 0;
  RCP<vector_type> v;
  Array<int> BCNodes;
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Dirichlet BC setup", bc_setup);
    for (int inode = 0; inode < numNodes; ++inode) {
      if (nodeOnBoundary(inode) && nodeIsOwned[inode]) {
        ++numBCNodes;
      }
    }
    // v: Vector for use in applying Dirichlet BCs.
    // v uses the owned (nonoverlapped) distribution.
    v = rcp (new vector_type (globalMapG, true));
    v->putScalar (STS::zero ());

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
          v->replaceLocalValue (iOwned, exactSolution (x,y,z));
        } // if node inode is on the boundary

        gl_exactSolVector->replaceLocalValue (iOwned, exactSolution (x,y,z));
        ++iOwned;
      } // if node inode is owned by my process
    } // for each node inode that my process can see
  } // Done setting up v (for Dirichlet BCs) and gl_exactSolVector

  /**********************************************************************************/
  /******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
  /**********************************************************************************/

  *out << "Building discretization matrix and right hand side" << endl;

  {
  TEUCHOS_FUNC_TIME_MONITOR_DIFF(
    "Matrix/RHS fill:  0-Total", fill_total);

  // Define desired workset size and count how many worksets there are
  // on this processor's mesh block
  //int desiredWorksetSize = numElems; // change to desired workset size!
  int desiredWorksetSize = 100;    // change to desired workset size!
  int numWorksets        = (numElems+desiredWorksetSize-1)/desiredWorksetSize;

  {
    Teuchos::OSTab tab2 (out);
    *out << "Desired workset size:             " << desiredWorksetSize << endl
         << "Number of worksets (per process): " << numWorksets << endl;
  }

  for (int workset = 0; workset < numWorksets; ++workset) {

    // Compute cell numbers where the workset starts and ends
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // When numElems is not divisible by desiredWorksetSize, the last
    // workset ends at numElems.
    worksetEnd = (worksetEnd <= numElems) ? worksetEnd : numElems;

    // Now we know the actual workset size and can allocate the array
    // for the cell nodes.
    int worksetSize = worksetEnd - worksetBegin;
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

    {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF(
      "Matrix/RHS fill:  1-Element discretization", elem_discretization);

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

    } // Element discretization timer

    /**********************************************************************************/
    /*                         Assemble into Global Matrix                            */
    /**********************************************************************************/

    {
      TEUCHOS_FUNC_TIME_MONITOR_DIFF(
        "Matrix/RHS fill:  2-Element assembly", elem_assembly_total);

      // "WORKSET CELL" loop: local cell ordinal is relative to numElems
      for (int cell = worksetBegin; cell < worksetEnd; ++cell) {

        // Compute cell ordinal relative to the current workset
        const int worksetCellOrdinal = cell - worksetBegin;

        // "CELL EQUATION" loop for the workset cell: cellRow is
        // relative to the cell DoF numbering.
        for (int cellRow = 0; cellRow < numFieldsG; ++cellRow) {
          int localRow  = elemToNode (cell, cellRow);
          ST sourceTermContribution = worksetRHS (worksetCellOrdinal, cellRow);

          {
            // TEUCHOS_FUNC_TIME_MONITOR_DIFF(
            //   "Matrix/RHS fill:  2a-RHS sum into local values",
            //   elem_rhs);
            rhsVector->sumIntoLocalValue (localRow, sourceTermContribution);
          }

          // "CELL VARIABLE" loop for the workset cell: sum entire element
          // stiff matrix contribution in one function call
          ArrayView<int> localColAV =
            arrayView<int> (&elemToNode(cell,0), numFieldsG);
          ArrayView<ST> operatorMatrixContributionAV =
            arrayView<ST> (&worksetStiffMatrix(worksetCellOrdinal,cellRow,0),
                           numFieldsG);
          {
            // TEUCHOS_FUNC_TIME_MONITOR_DIFF(
            //   "Matrix/RHS fill:  2b-Matrix sum into local values",
            //   elem_matrix);
            StiffMatrix->sumIntoLocalValues (localRow, localColAV,
                                             operatorMatrixContributionAV);
          }
        }// *** cell row loop ***
      }// *** workset cell loop **
    } // *** stop timer ***
  }// *** workset loop ***

  //////////////////////////////////////////////////////////////////////////////
  // Export sparse matrix and right-hand side from overlapping row Map
  // to owned (nonoverlapping) row Map.
  //////////////////////////////////////////////////////////////////////////////

  *out << "Exporting matrix and right-hand side from overlapped to owned Map" << endl;

  gl_StiffMatrix->setAllToScalar (STS::zero ());
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF(
      "Matrix/RHS fill:  3-Matrix export", matrix_export);
    gl_StiffMatrix->doExport (*StiffMatrix, *exporter, Tpetra::ADD);
  }
  // If target of export has static graph, no need to do
  // setAllToScalar(0.0); export will clobber values.
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF(
      "Matrix/RHS fill:  4-Matrix fill complete", matrix_fill_complete);
    gl_StiffMatrix->fillComplete ();
  }

  gl_rhsVector->putScalar (STS::zero ());
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF(
      "Matrix/RHS fill:  5-RHS export", rhs_export);
    gl_rhsVector->doExport (*rhsVector, *exporter, Tpetra::ADD);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Adjust matrix for boundary conditions
  //////////////////////////////////////////////////////////////////////////////

  *out << "Adjusting matrix and right-hand side for BCs" << endl;

  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF(
      "Matrix/RHS fill:  6-Dirichlet BC assembly", bc_assembly);

    // Apply owned stiffness matrix to v: rhs := A*v
    RCP<multivector_type> rhsDir =
      rcp (new multivector_type (globalMapG, true));
    gl_StiffMatrix->apply (*v.getConst (), *rhsDir);
    // Update right-hand side
    gl_rhsVector->update (as<ST> (-STS::one ()), *rhsDir, STS::one ());

    // Get a const view of the vector's local data.
    ArrayRCP<const ST> vArrayRCP = v->getData (0);

    // Adjust rhs due to Dirichlet boundary conditions.
    for (int inode = 0; inode < numNodes; ++inode) {
      if (nodeIsOwned[inode]) {
        if (nodeOnBoundary (inode)) {
          // Get global node ID
          const GO gni = as<GO> (globalNodeIds[inode]);
          const LO lidT = globalMapG->getLocalElement (gni);
          ST v_valT = vArrayRCP[lidT];
          gl_rhsVector->replaceGlobalValue (gni, v_valT);
        }
      }
    }

    // Zero out rows and columns of stiffness matrix corresponding to
    // Dirichlet edges and add one to diagonal.  The following is the
    // Tpetra analog of Apply_OAZToMatrix().
    //
    // Reenable changes to the values and structure of the global
    // stiffness matrix.
    gl_StiffMatrix->resumeFill ();

    // Find the local column numbers to nuke
    RCP<const map_type> ColMap = gl_StiffMatrix->getColMap ();
    RCP<const map_type> globalMap =
      rcp (new map_type (gl_StiffMatrix->getGlobalNumCols (), 0, comm,
                         Tpetra::GloballyDistributed));

    // Create the exporter from this process' column Map to the global
    // 1-1 column map. (???)
    RCP<const export_type> bdyExporter =
      rcp (new export_type (ColMap, globalMap));
    // Create a vector of global column indices to which we will export
    RCP<Tpetra::Vector<int, LO, GO, Node> > globColsToZeroT =
      rcp (new Tpetra::Vector<int, LO, GO, Node> (globalMap));
    // Create a vector of local column indices from which we will export
    RCP<Tpetra::Vector<int, LO, GO, Node> > myColsToZeroT =
      rcp (new Tpetra::Vector<int, LO, GO, Node> (ColMap));
    myColsToZeroT->putScalar (0);

    // Flag (set to 1) all local columns corresponding to the local
    // rows specified.
    for (int i = 0; i < numBCNodes; ++i) {
      const GO globalRow = gl_StiffMatrix->getRowMap ()->getGlobalElement (BCNodes[i]);
      const LO localCol = gl_StiffMatrix->getColMap ()->getLocalElement (globalRow);
      // Tpetra::Vector<int, ...> works just like
      // Tpetra::Vector<double, ...>.  Epetra has a separate
      // Epetra_IntVector class for ints.
      myColsToZeroT->replaceLocalValue (localCol, 1);
    }

    // Export to the global column map.
    globColsToZeroT->doExport (*myColsToZeroT, *bdyExporter, Tpetra::ADD);
    // Import from the global column map to the local column map.
    myColsToZeroT->doImport (*globColsToZeroT, *bdyExporter, Tpetra::INSERT);

    ArrayRCP<const int> myColsToZeroArrayRCP = myColsToZeroT->getData(0);
    size_t NumEntries = 0;

    // Zero the columns corresponding to Dirichlet BCs.
    typename sparse_matrix_type::nonconst_local_inds_host_view_type indices("indices", 1);
    typename sparse_matrix_type::nonconst_values_host_view_type values("values", 1);
    for (LO i = 0; i < as<int> (gl_StiffMatrix->getLocalNumRows ()); ++i) {
      NumEntries = gl_StiffMatrix->getNumEntriesInLocalRow (i);
      Kokkos::resize(indices, NumEntries);
      Kokkos::resize(values, NumEntries);
      gl_StiffMatrix->getLocalRowCopy (i, indices, values, NumEntries);
      for (int j = 0; j < as<int> (NumEntries); ++j) {
        if (myColsToZeroArrayRCP[indices(j)] == 1)
          values(j) = STS::zero ();
      }
      gl_StiffMatrix->replaceLocalValues (i, indices, values);
    } // for each (local) row of the global stiffness matrix

    // Zero the rows and add ones to diagonal.
    for (int i = 0; i < numBCNodes; ++i) {
      NumEntries = gl_StiffMatrix->getNumEntriesInLocalRow (BCNodes[i]);
      Kokkos::resize(indices, NumEntries);
      Kokkos::resize(values, NumEntries);
      gl_StiffMatrix->getLocalRowCopy (BCNodes[i], indices, values, NumEntries);
      const GO globalRow = gl_StiffMatrix->getRowMap ()->getGlobalElement (BCNodes[i]);
      const LO localCol = gl_StiffMatrix->getColMap ()->getLocalElement (globalRow);
      for (int j = 0; j < as<int> (NumEntries); ++j) {
        values(j) = STS::zero ();
        if (indices(j) == localCol) {
          values(j) = STS::one ();
        }
      } // for each entry in the current row
      gl_StiffMatrix->replaceLocalValues (BCNodes[i], indices, values);
    } // for each BC node
  }

  *out << "Calling fillComplete() on owned-Map stiffness matrix" << endl;

  // We're done modifying the owned stiffness matrix.
  gl_StiffMatrix->fillComplete ();

  } // Matrix/RHS fill timing

  //
  // We're done with assembly, so we can delete the mesh.
  //
  delete [] nodeIsOwned;
  nodeIsOwned = NULL;
  Delete_Pamgen_Mesh ();

  // Create vector to store approximate solution, and set initial guess.
  RCP<vector_type> gl_approxSolVector = rcp (new vector_type (globalMapG));
  gl_approxSolVector->putScalar (STS::zero ());

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

  RCP<vector_type> R = rcp (new vector_type (*B)); // R := B
  // R := 1.0*R - 1.0*A*X_exact.
  A->apply (*X_exact, *R, Teuchos::NO_TRANS, -STS::one(), STS::one());

  std::vector<MT> results (3);
  results[0] = R->norm2 ();
  results[1] = B->norm2 ();
  results[2] = A->getFrobeniusNorm ();
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


using ::TrilinosCouplings::IntrepidPoissonExample::getMaterialTensorOffDiagonalValue;

/** \brief  User-defined material tensor.

    Evaluate the tensor using operator().  Its arguments are:

    \param  material    [out]   3 x 3 material tensor evaluated at (x,y,z)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \warning Symmetric and positive definite tensor is required for every (x,y,z).
*/
template<typename Scalar>
class MaterialTensor {
public:
  MaterialTensor (const double offDiagVal) :
    offDiagVal_ (Scalar (offDiagVal))
  {}

  void
  operator () (Scalar material[][3],
               const Scalar& x,
               const Scalar& y,
               const Scalar& z) const
  {
    typedef Teuchos::ScalarTraits<Scalar> STS;

    // We go through this trouble to make numbers, because Scalar
    // isn't necessarily double.  It could be some automatic
    // differentiation type.
    const Scalar zero = STS::zero ();
    const Scalar one = STS::one ();

    // You can use the value of offDiagVal_ to control the iteration
    // count.  The iteration counts below are for Belos' GMRES with no
    // preconditioning, using the default problem size.
    // setMaterialTensorOffDiagonalValue() sets the value of this
    // parameter.
    //
    // Classical elasticity assumes a symmetric material tensor.  I
    // suppose one could solve Poisson's equation with an unsymmetric
    // material tensor, but I'm not sure what that would mean.

    // offDiagVal_ = -5/4: 209 iterations (CG breaks!)
    // offDiagVal_ = -1/2: 47 iterations
    // offDiagVal_ = 0: 40 iterations (CG works)
    // offDiagVal_ = 1/2: 46 iterations
    // offDiagVal_ = 3/4: 47 iterations
    // offDiagVal_ = 1: 59 iterations
    // offDiagVal_ = 5/4: 183 iterations
    // offDiagVal_ = 3/2: 491 iterations
    // offDiagVal_ = 2: 939 iterations (CG breaks!)
    material[0][0] = one;
    material[0][1] = zero;
    material[0][2] = offDiagVal_;

    material[1][0] = zero;
    material[1][1] = one;
    material[1][2] = zero;

    material[2][0] = offDiagVal_;
    material[2][1] = zero;
    material[2][2] = one;
  }

private:
  const Scalar offDiagVal_;
};

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
  Scalar flux[3] = {0.0, 0.0, 0.0};
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
  MaterialTensor<Scalar> matTens (getMaterialTensorOffDiagonalValue ());
  matTens (material, x, y, z);

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

  MaterialTensor<scalar_type> matTens (getMaterialTensorOffDiagonalValue ());
  for (int cell = 0; cell < numWorksetCells; ++cell) {
    for (int pt = 0; pt < numPoints; ++pt) {
      scalar_type x = evaluationPoints(cell, pt, 0);
      scalar_type y = evaluationPoints(cell, pt, 1);
      scalar_type z = evaluationPoints(cell, pt, 2);

      matTens (material, x, y, z);

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

} // namespace TpetraIntrepidPoissonExample
} // namespace TrilinosCouplings
