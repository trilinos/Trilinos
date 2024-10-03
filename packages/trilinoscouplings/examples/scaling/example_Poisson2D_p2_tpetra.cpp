// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   example_Poisson.cpp
    \brief  Example solution of a Poisson equation on a quad mesh using
    nodal (Hgrad) elements.

    This example uses the following Trilinos packages:
    \li     Pamgen to generate a Quad mesh.
    \li     Sacado to form the source term from user-specified manufactured solution.
    \li     Intrepid to build the discretization matrix and right-hand side.
    \li     Tpetra to handle the global matrix and vector.
    \li     Isorropia to partition the matrix. (Optional)
    \li     MueLu to solve the linear system.


    \verbatim

    Poisson system:

    div A grad u = f in Omega
    u = g on Gamma

    where
    A is a symmetric, positive definite material tensor
    f is a given source term

    Corresponding discrete linear system for nodal coefficients(x):

    Kx = b

    K - HGrad stiffness matrix
    b - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal, K. Peterson, D. Hensinger, C. Siefert.

    \remark Usage:
    \code   ./TrilinosCouplings_examples_scaling_example_Poisson.exe \endcode

    \remark Example driver requires input file named Poisson2D.xml with Pamgen
    formatted mesh description and settings for Isorropia (a version
    is included in the Trilinos repository with this driver).

    \remark The exact solution (u) and material tensor (A) are set in the
    functions "exactSolution" and "materialTensor" and may be
    modified by the user.

*/

/*** Uncomment if you would like output data for plotting ***/
//#define DUMP_DATA

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

// TrilinosCouplings includes
#include "TrilinosCouplings_config.h"
#include "TrilinosCouplings_Pamgen_Utils.hpp"

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <MatrixMarket_Tpetra.hpp>

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// Pamgen includes
#include "create_inline_mesh.h"
#include "pamgen_im_exodusII_l.h"
#include "pamgen_im_ne_nemesisI_l.h"
#include "pamgen_extras.h"

#ifdef HAVE_TRILINOSCOUPLINGS_MUELU
#  include "MueLu.hpp"
#  include "MueLu_ParameterListInterpreter.hpp"
#  include "MueLu_TpetraOperator.hpp"
#  include "MueLu_Utilities.hpp"
#  include "MueLu_HierarchyManager.hpp"
#  include "MueLu_FactoryManagerBase.hpp"
#endif // HAVE_TRILINOSCOUPLINGS_MUELU

#ifdef HAVE_INTREPID_KOKKOS
#include "Sacado.hpp"
#else
// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#endif

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosFixedPointSolMgr.hpp>

#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp

using namespace std;
using namespace Intrepid;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::TimeMonitor;
using Teuchos::ParameterList;



/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef double scalar_type;
typedef Teuchos::ScalarTraits<scalar_type> ScalarTraits;
using local_ordinal_type = Tpetra::Map<>::local_ordinal_type;
using global_ordinal_type = Tpetra::Map<>::global_ordinal_type;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;
typedef Sacado::Fad::SFad<double,2>      Fad2; //# ind. vars fixed at 2
typedef Intrepid::FunctionSpaceTools     IntrepidFSTools;
typedef Intrepid::RealSpaceTools<double> IntrepidRSTools;
typedef Intrepid::CellTools<double>      IntrepidCTools;

typedef Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,NO>    operator_type;
typedef Tpetra::CrsMatrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>   crs_matrix_type;
typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,NO>      vector_type;
typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,NO> multivector_type;
typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,NO>                     driver_map_type;
Tpetra::global_size_t INVALID_GO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

typedef Belos::LinearProblem<scalar_type, multivector_type, operator_type> linear_problem_type;

typedef MueLu::TpetraOperator<scalar_type,local_ordinal_type,global_ordinal_type,NO> muelu_tpetra_operator;

//typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;

#define TC_sumAll(rcpComm, in, out) \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out))
#define TC_minAll(rcpComm, in, out) \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out))
#define TC_maxAll(rcpComm, in, out) \
    Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out))


// forward declarations

void PromoteMesh(const int degree, const FieldContainer<int> & P1_elemToNode, const FieldContainer<double> & P1_nodeCoord, const FieldContainer<double> & P1_edgeCoord,  const FieldContainer<int> & P1_elemToEdge,  const FieldContainer<int> & P1_elemToEdgeOrient, const FieldContainer<int> & P1_nodeOnBoundary,
                 FieldContainer<int> & P2_elemToNode, FieldContainer<double> & P2_nodeCoord,FieldContainer<int> & P2_nodeOnBoundary, std::vector<int> & P2_edgeNodes, std::vector<int> & P2_cellNodes);

void GenerateEdgeEnumeration(const FieldContainer<int> & elemToNode, const FieldContainer<double> & nodeCoord, FieldContainer<int> & elemToEdge, FieldContainer<int> & elemToEdgeOrient, FieldContainer<double> & edgeCoord);

void CreateP1MeshFromP2Mesh(const FieldContainer<int> & P2_elemToNode, FieldContainer<int> &aux_P1_elemToNode);

void CreateLinearSystem(int numWorkSets,
                        int desiredWorksetSize,
                        FieldContainer<int> const &elemToNode,
                        FieldContainer<double> const &nodeCoord,
                        FieldContainer<double> const &cubPoints,
                        FieldContainer<double> const &cubWeights,
                        FieldContainer<double> const &HGBGrads,
                        FieldContainer<double> const &HGBValues,
                        std::vector<global_ordinal_type>       const &globalNodeIds,
                        shards::CellTopology const &cellType,
                        crs_matrix_type &StiffMatrix,
                        RCP<multivector_type> &rhsVector,
                        std::string &msg
                        );

void GenerateLinearCoarsening_p2_to_p1(const FieldContainer<int> & P2_elemToNode, RCP<driver_map_type> &P1_map, RCP<driver_map_type> & P2_map,RCP<crs_matrix_type> & P);

void GenerateIdentityCoarsening_p2_to_p1(const FieldContainer<int> & P2_elemToNode,
                      const FieldContainer<int> & P1_elemToNode,
                      RCP<const driver_map_type> const & P1_map_aux, RCP<const driver_map_type> const &P2_map,
                      RCP<crs_matrix_type> & Interpolation,
                      RCP<crs_matrix_type> & Restriction);

void Apply_Dirichlet_BCs(std::vector<int> &BCNodes, crs_matrix_type & A, multivector_type & x, multivector_type & b, const multivector_type & soln);

/**********************************************************************************/
/***************** FUNCTION DECLARATION FOR MUELU PRECONDITIONER *********************/
/**********************************************************************************/

/** \brief  MueLu Preconditioner

    \param  ProblemType        [in]    problem type
    \param  amgList            [in]    MueLu parameter list
    \param  A                  [in]    discrete operator matrix
    \param  xexact             [in]    exact solution
    \param  b                  [in]    right-hand-side vector
    \param  uh                 [out]   solution vector
    \param  TotalErrorResidual [out]   error residual
    \param  TotalErrorExactSol [out]   error in uh

*/
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
                                 ParameterList   & amgList,
                                 RCP<crs_matrix_type>   const & A,
                                 RCP<crs_matrix_type>   const & P,
                                 RCP<crs_matrix_type>   const & R,
                                 RCP<multivector_type> const & xexact,
                                 RCP<multivector_type> & b,
                                 RCP<multivector_type> & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol,
				 std::string &amgType,
				 std::string &solveType);



/**********************************************************************************/
/******** FUNCTION DECLARATIONS FOR EXACT SOLUTION AND SOURCE TERMS ***************/
/**********************************************************************************/

/** \brief  User-defined exact solution.

    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point

    \return Value of the exact solution at (x,y)
*/
template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y);

/** \brief  User-defined material tensor.

    \param  material    [out]   3 x 3 material tensor evaluated at (x,y)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point

    \warning Symmetric and positive definite tensor is required for every (x,y).
*/
template<typename Scalar>
void materialTensor(Scalar material[][2], const Scalar&  x, const Scalar&  y);

/** \brief  Computes gradient of the exact solution. Requires user-defined exact solution.

    \param  gradExact  [out]   gradient of the exact solution evaluated at (x,y)
    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
*/
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[2], const Scalar& x, const Scalar& y);


/** \brief Computes source term: f = -div(A.grad u).  Requires user-defined exact solution
    and material tensor.

    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point

    \return Source term corresponding to the user-defined exact solution evaluated at (x,y)
*/
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y);


/** \brief Computation of the material tensor at array of points in physical space.

    \param worksetMaterialValues      [out]     Rank-2, 3 or 4 array with dimensions (C,P), (C,P,D) or (C,P,D,D)
    with the values of the material tensor
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        worksetMaterialValues,
                            const ArrayIn &   evaluationPoints);


/** \brief Computation of the source term at array of points in physical space.

    \param sourceTermValues           [out]     Rank-2 (C,P) array with the values of the source term
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints);

/** \brief Computation of the exact solution at array of points in physical space.

    \param exactSolutionValues        [out]     Rank-2 (C,P) array with the values of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints);


/** \brief Computation of the gradient of the exact solution at array of points in physical space.

    \param exactSolutionGradValues    [out]     Rank-3 (C,P,D) array with the values of the gradient of the exact solution
    \param evaluationPoints           [in]      Rank-3 (C,P,D) array with the evaluation points in physical frame
*/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints);


/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/

int main(int argc, char *argv[]) {

  int error = 0;
  int numProcs=1;
  int rank=0;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  RCP<const Teuchos::Comm<int> > Comm = Teuchos::DefaultComm<int>::getComm();
  int MyPID = Comm->getRank();

  RCP<TimeMonitor> tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Pamgen Setup")));

  //Check number of arguments
  if (argc > 3) {
    std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
    std::cout <<"Usage:\n\n";
    std::cout <<"  ./TrilinosCouplings_examples_scaling_Example_Poisson.exe [meshfile.xml] [solver.xml]\n\n";
    std::cout <<"   meshfile.xml(optional) - xml file with description of Pamgen mesh\n\n";
    std::cout <<"   solver.xml(optional) - xml file with MueLu solver options\n\n";
    exit(1);
  }

  if (MyPID == 0){
    std::cout \
      << "===============================================================================\n" \
      << "|                                                                             |\n" \
      << "|          Example: Solve Poisson Equation on Hexahedral Mesh                 |\n" \
      << "|                                                                             |\n" \
      << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                    |\n" \
      << "|                      Denis Ridzal  (dridzal@sandia.gov),                    |\n" \
      << "|                      Kara Peterson (kjpeter@sandia.gov).                    |\n" \
      << "|                                                                             |\n" \
      << "|  Intrepid's website:   http://trilinos.sandia.gov/packages/intrepid         |\n" \
      << "|  Pamgen's website:     http://trilinos.sandia.gov/packages/pamgen           |\n" \
      << "|  ML's website:         http://trilinos.sandia.gov/packages/ml               |\n" \
      << "|  Isorropia's website:  http://trilinos.sandia.gov/packages/isorropia        |\n" \
      << "|  Trilinos website:     http://trilinos.sandia.gov                           |\n" \
      << "|                                                                             |\n" \
      << "===============================================================================\n";
  }


#ifdef HAVE_MPI
  if (MyPID == 0) {
    std::cout << "PARALLEL executable \n";
  }
#else
  if (MyPID == 0) {
    std::cout << "SERIAL executable \n";
  }
#endif

  /**********************************************************************************/
  /********************************** GET XML INPUTS ********************************/
  /**********************************************************************************/

  // Command line for xml file, otherwise use default
  std::string   xmlMeshInFileName, xmlSolverInFileName;
  if(argc>=2) xmlMeshInFileName=string(argv[1]);
  else xmlMeshInFileName="Poisson2D.xml";
  if(argc>=3) xmlSolverInFileName=string(argv[2]);

  // Read xml file into parameter list
  ParameterList inputMeshList;
  ParameterList inputSolverList;

  if(xmlMeshInFileName.length()) {
    if (MyPID == 0) {
      std::cout << "\nReading parameter list from the XML file \""<<xmlMeshInFileName<<"\" ...\n\n";
    }
    Teuchos::updateParametersFromXmlFile (xmlMeshInFileName, Teuchos::ptr (&inputMeshList));
    if (MyPID == 0) {
      inputMeshList.print(std::cout,2,true,true);
      std::cout << "\n";
    }
  }
  else
    {
      std::cout << "Cannot read input file: " << xmlMeshInFileName << "\n";
      return 0;
    }

  if(xmlSolverInFileName.length()) {
    if (MyPID == 0)
      std::cout << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n\n";
    Teuchos::updateParametersFromXmlFile(xmlSolverInFileName, Teuchos::inoutArg(inputSolverList));
  } else if (MyPID == 0) std::cout << "Using default solver values ..." << std::endl;

  // Get pamgen mesh definition
  std::string meshInput = Teuchos::getParameter<std::string>(inputMeshList,"meshInput");

  // Get Isorropia and Zoltan parameters.
  ParameterList iso_paramlist = inputMeshList.sublist
    ("Isorropia Input") ;
  if (MyPID == 0) {
    std::cout << "Isorropia/Zoltan parameters" << std::endl;
    iso_paramlist.print(std::cout,2,true,true);
  }


  /**********************************************************************************/
  /***************************** GET CELL TOPOLOGY **********************************/
  /**********************************************************************************/

  // Get cell topology for base hexahedron
  shards::CellTopology P1_cellType(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  shards::CellTopology P2_cellType(shards::getCellTopologyData<shards::Quadrilateral<9> >() );
  assert(P1_cellType.getDimension() == P2_cellType.getDimension());

  // Get dimensions
  int P1_numNodesPerElem = P1_cellType.getNodeCount();
  int P2_numNodesPerElem = P2_cellType.getNodeCount();
  int spaceDim = P1_cellType.getDimension();
  int dim = 2;

  /**********************************************************************************/
  /******************************* GENERATE MESH ************************************/
  /**********************************************************************************/

  if (MyPID == 0) {
    std::cout << "Generating mesh ... \n\n";
  }

  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  long long cr_result = Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);
  TrilinosCouplings::pamgen_error_check(std::cout,cr_result);

  string msg("Poisson: ");
  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Mesh queries")));

  // Get mesh size info
  char title[100];
  long long numDim;
  long long numNodes;
  long long numElems;
  long long numElemBlk;
  long long numNodeSets;
  long long numSideSets;
  int id = 0;

  im_ex_get_init_l(id, title, &numDim, &numNodes,
                   &numElems, &numElemBlk, &numNodeSets,
                   &numSideSets);

  long long numNodesGlobal;
  long long numElemsGlobal;
  long long numElemBlkGlobal;
  long long numNodeSetsGlobal;
  long long numSideSetsGlobal;

  im_ne_get_init_global_l(id, &numNodesGlobal, &numElemsGlobal,
                          &numElemBlkGlobal, &numNodeSetsGlobal,
                          &numSideSetsGlobal);

  // Print mesh information
  if (MyPID == 0){
    std::cout << " Number of Global Elements: " << numElemsGlobal << " \n";
    std::cout << " Number of Global Nodes: " << numNodesGlobal << " \n\n";
  }

  long long * block_ids = new long long [numElemBlk];
  error += im_ex_get_elem_blk_ids_l(id, block_ids);


  long long  *nodes_per_element   = new long long [numElemBlk];
  long long  *element_attributes  = new long long [numElemBlk];
  long long  *elements            = new long long [numElemBlk];
  char      **element_types       = new char * [numElemBlk];
  long long **elmt_node_linkage   = new long long * [numElemBlk];


  for(long long i = 0; i < numElemBlk; i ++){
    element_types[i] = new char [MAX_STR_LENGTH + 1];
    error += im_ex_get_elem_block_l(id,
                                    block_ids[i],
                                    element_types[i],
                                    (long long*)&(elements[i]),
                                    (long long*)&(nodes_per_element[i]),
                                    (long long*)&(element_attributes[i]));
  }

  /*connectivity*/
  for(long long b = 0; b < numElemBlk; b++){
    elmt_node_linkage[b] =  new long long [nodes_per_element[b]*elements[b]];
    error += im_ex_get_elem_conn_l(id,block_ids[b],elmt_node_linkage[b]);
  }

  // Get node-element connectivity
  int telct = 0;
  FieldContainer<int> P1_elemToNode(numElems,P1_numNodesPerElem);
  for(long long b = 0; b < numElemBlk; b++){
    for(long long el = 0; el < elements[b]; el++){
      for (int j=0; j<P1_numNodesPerElem; j++) {
        P1_elemToNode(telct,j) = elmt_node_linkage[b][el*P1_numNodesPerElem + j]-1;
      }
      telct ++;
    }
  }

  // Read node coordinates and place in field container
  FieldContainer<double> P1_nodeCoord(numNodes,dim);
  double * nodeCoordx = new double [numNodes];
  double * nodeCoordy = new double [numNodes];
  im_ex_get_coord_l(id,nodeCoordx,nodeCoordy,0);
  for (int i=0; i<numNodes; i++) {
    P1_nodeCoord(i,0)=nodeCoordx[i];
    P1_nodeCoord(i,1)=nodeCoordy[i];
  }

  /*parallel info*/
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

  if(num_node_comm_maps > 0){
    node_comm_proc_ids   = new long long  [num_node_comm_maps];
    node_cmap_node_cnts  = new long long  [num_node_comm_maps];
    node_cmap_ids        = new long long  [num_node_comm_maps];
    comm_node_ids        = new long long* [num_node_comm_maps];
    comm_node_proc_ids   = new long long* [num_node_comm_maps];

    long long *  elem_cmap_ids        = new long long [num_elem_comm_maps];
    long long *  elem_cmap_elem_cnts  = new long long [num_elem_comm_maps];


    if ( im_ne_get_cmap_params_l( id,
                                  node_cmap_ids,
                                  (long long*)node_cmap_node_cnts,
                                  elem_cmap_ids,
                                  (long long*)elem_cmap_elem_cnts,
                                  0/*not used proc_id*/ ) < 0 )++error;

    for(long long j = 0; j < num_node_comm_maps; j++) {
      comm_node_ids[j]       = new long long [node_cmap_node_cnts[j]];
      comm_node_proc_ids[j]  = new long long [node_cmap_node_cnts[j]];
      if ( im_ne_get_node_cmap_l( id,
                                  node_cmap_ids[j],
                                  comm_node_ids[j],
                                  comm_node_proc_ids[j],
                                  0/*not used proc_id*/ ) < 0 )++error;
      node_comm_proc_ids[j] = comm_node_proc_ids[j][0];
    }

    delete [] elem_cmap_ids;
    delete [] elem_cmap_elem_cnts;
  }

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Global Node Nums")));


  //Calculate global node ids
  long long * P1_globalNodeIds = new long long[numNodes];
  bool * P1_nodeIsOwned = new bool[numNodes];

  calc_global_node_ids(P1_globalNodeIds,
                       P1_nodeIsOwned,
                       numNodes,
                       num_node_comm_maps,
                       node_cmap_node_cnts,
                       node_comm_proc_ids,
                       comm_node_ids,
                       rank);


  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Boundary Conds")));

  // Container indicating whether a node is on the boundary (1-yes 0-no)
  FieldContainer<int> P1_nodeOnBoundary(numNodes);

  // Get boundary (side set) information
  long long * sideSetIds = new long long [numSideSets];
  long long numSidesInSet;
  long long numDFinSet;
  im_ex_get_side_set_ids_l(id,sideSetIds);
  for (int i=0; i<numSideSets; i++) {
    im_ex_get_side_set_param_l(id,sideSetIds[i],&numSidesInSet,&numDFinSet);
    if (numSidesInSet > 0){
      long long * sideSetElemList = new long long [numSidesInSet];
      long long * sideSetSideList = new long long [numSidesInSet];
      im_ex_get_side_set_l(id,sideSetIds[i],sideSetElemList,sideSetSideList);
      for (int j=0; j<numSidesInSet; j++) {

        int sideNode0 = P1_cellType.getNodeMap(1,sideSetSideList[j]-1,0);
        int sideNode1 = P1_cellType.getNodeMap(1,sideSetSideList[j]-1,1);

        P1_nodeOnBoundary(P1_elemToNode(sideSetElemList[j]-1,sideNode0))=1;
        P1_nodeOnBoundary(P1_elemToNode(sideSetElemList[j]-1,sideNode1))=1;
      }
      delete [] sideSetElemList;
      delete [] sideSetSideList;
    }
  }
  delete [] sideSetIds;

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Getting cubature")));



  // Enumerate edges 
  // NOTE: Only correct in serial
  FieldContainer<int> P1_elemToEdge(numElems,4);// Because quads
  FieldContainer<int> P1_elemToEdgeOrient(numElems,4);
  FieldContainer<double> P1_edgeCoord(1,dim);//will be resized  
  GenerateEdgeEnumeration(P1_elemToNode, P1_nodeCoord, P1_elemToEdge,P1_elemToEdgeOrient,P1_edgeCoord);


  // Generate higher order mesh
  // NOTE: Only correct in serial
  int P2_numNodes = numNodes + P1_edgeCoord.dimension(0) + numElems;
  FieldContainer<int>    elemToNode(numElems,9); //because quads
  FieldContainer<double> nodeCoord(P2_numNodes,dim);
  FieldContainer<int>   nodeOnBoundary(P2_numNodes);
  std::vector<int> P2_edgeNodes;
  std::vector<int> P2_cellNodes;
  PromoteMesh(2,P1_elemToNode,P1_nodeCoord,P1_edgeCoord,P1_elemToEdge,P1_elemToEdgeOrient,P1_nodeOnBoundary,
              elemToNode, nodeCoord, nodeOnBoundary, P2_edgeNodes, P2_cellNodes);

  long long numElems_aux = numElems*4;  //4 P1 elements per P2 element in auxiliary mesh
  FieldContainer<int> aux_P1_elemToNode(numElems_aux,P1_numNodesPerElem); //4 P1 elements per P2 element
  CreateP1MeshFromP2Mesh(elemToNode, aux_P1_elemToNode);

  // Only works in serial
  std::vector<bool>P2_nodeIsOwned(P2_numNodes,true);
  std::vector<global_ordinal_type>P2_globalNodeIds(P2_numNodes);
  for(int i=0; i<P2_numNodes; i++)
    P2_globalNodeIds[i]=i;

  std::vector<double> P2_nodeCoordx(P2_numNodes);
  std::vector<double> P2_nodeCoordy(P2_numNodes);
  for (int i=0; i<P2_numNodes; i++) {
    P2_nodeCoordx[i] = nodeCoord(i,0);
    P2_nodeCoordy[i] = nodeCoord(i,1);
  }

  // Reset constants
  int P1_numNodes =numNodes;
  numNodes = P2_numNodes;
  numNodesGlobal = numNodes;

  // Print mesh information
  if (MyPID == 0){
    std::cout << " Number of P2 Global Elements: " << numElemsGlobal << " \n";
    std::cout << " Number of P2 Global Nodes: " << numNodesGlobal << " \n";
    std::cout << " Number of faux P1 Global Elements: " << aux_P1_elemToNode.dimension(0) << " \n\n";
  }


  /**********************************************************************************/
  /********************************* GET CUBATURE ***********************************/
  /**********************************************************************************/

  // Get numerical integration points and weights
  DefaultCubatureFactory<double>  cubFactory;
  int cubDegree = 4;
  RCP<Cubature<double> > myCub = cubFactory.create(P2_cellType, cubDegree);

  int cubDim       = myCub->getDimension();
  int numCubPoints = myCub->getNumPoints();

  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);

  myCub->getCubature(cubPoints, cubWeights);

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Getting basis")));



  /**********************************************************************************/
  /*********************************** GET BASIS ************************************/
  /**********************************************************************************/

  // Define basis
  Basis_HGRAD_QUAD_C2_FEM<double, FieldContainer<double> > myHGradBasis;
  int numFieldsG = myHGradBasis.getCardinality();
  FieldContainer<double> HGBValues(numFieldsG, numCubPoints);
  FieldContainer<double> HGBGrads(numFieldsG, numCubPoints, spaceDim);

  // Evaluate basis values and gradients at cubature points
  myHGradBasis.getValues(HGBValues, cubPoints, OPERATOR_VALUE);
  myHGradBasis.getValues(HGBGrads, cubPoints, OPERATOR_GRAD);

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Build global maps")));


  /**********************************************************************************/
  /********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
  /**********************************************************************************/
  // Count owned nodes (P2)
  int P2_ownedNodes=0;
  for(int i=0;i<numNodes;i++)
    if(P2_nodeIsOwned[i]) P2_ownedNodes++;

  // Build a list of the OWNED global ids...
  // NTS: will need to switch back to long long
  std::vector<global_ordinal_type> P2_ownedGIDs(P2_ownedNodes);
  int oidx=0;
  for(int i=0;i<numNodes;i++)
    if(P2_nodeIsOwned[i]){
      P2_ownedGIDs[oidx]=(global_ordinal_type)P2_globalNodeIds[i];
      oidx++;
    }

  // Count owned nodes (P1)
  int P1_ownedNodes=0;
  for(int i=0;i<P1_numNodes;i++)
    if(P1_nodeIsOwned[i]) P1_ownedNodes++;
  std::vector<global_ordinal_type> P1_ownedGIDs(P1_ownedNodes);
  oidx=0;
  for(int i=0;i<P1_numNodes;i++)
    if(P1_nodeIsOwned[i]){
      P1_ownedGIDs[oidx]=(global_ordinal_type)P1_globalNodeIds[i];
      //nodeSeeds[oidx] = (int)P1_globalNodeIds[i];
      oidx++;
    }

  //seed points for block relaxation
  // numNodes is the number of nodes in the P2 mesh
  ArrayRCP<global_ordinal_type> nodeSeeds(numNodes,Teuchos::OrdinalTraits<local_ordinal_type>::invalid());
  //P1 nodes
  oidx=0;
  for(int i=0;i<P1_numNodes;i++) {
    if(P1_nodeIsOwned[i]){
      nodeSeeds[(int)P1_globalNodeIds[i]] = oidx;
      oidx++;
    }
  }
  int numNodeSeeds = oidx;

  //edge nodes only
  ArrayRCP<global_ordinal_type> edgeSeeds(numNodes,Teuchos::OrdinalTraits<local_ordinal_type>::invalid());
  for(size_t i=0;i<P2_edgeNodes.size();i++) {
    edgeSeeds[P2_edgeNodes[i]] = i;
  }
  int numEdgeSeeds = P2_edgeNodes.size();

  //cell nodes only
  ArrayRCP<global_ordinal_type> cellSeeds(numNodes,Teuchos::OrdinalTraits<local_ordinal_type>::invalid());
  for(size_t i=0;i<P2_cellNodes.size();i++) {
    cellSeeds[P2_cellNodes[i]] = i;
  }
  int numCellSeeds = P2_cellNodes.size();

       
  // Generate map for nodes
  RCP<driver_map_type> globalMapG = rcp(new driver_map_type(INVALID_GO,&P2_ownedGIDs[0],P2_ownedNodes,0,Comm));
    
  // Generate p1 map
  RCP<driver_map_type> P1_globalMap = rcp(new driver_map_type(INVALID_GO,&P1_ownedGIDs[0],P1_ownedNodes,0,Comm));

  // Genetrate P2-to-P1 coarsening.
  if (inputSolverList.isParameter("aux P1") && inputSolverList.isParameter("linear P1"))
    throw std::runtime_error("Can only specify \"aux P1\" or \"linear P1\", not both.");
  RCP<crs_matrix_type> P_linear;
  if (inputSolverList.isParameter("linear P1")) {
    GenerateLinearCoarsening_p2_to_p1(elemToNode,P1_globalMap,globalMapG,P_linear);
    inputSolverList.remove("linear P1"); //even though LevelWrap happily accepts this parameter
  }

  crs_matrix_type StiffMatrix(globalMapG, 20*numFieldsG);
  RCP<multivector_type> rhsVector = rcp(new multivector_type(globalMapG,1));
  RCP<multivector_type> femCoefficients = rcp(new multivector_type(globalMapG,1));

  tm = Teuchos::null;


#ifdef DUMP_DATA_OLD
  /**********************************************************************************/
  /**** PUT COORDINATES AND NODAL VALUES IN ARRAYS FOR OUTPUT (FOR PLOTTING ONLY) ***/
  /**********************************************************************************/

  // Put coordinates in multivector for output
  multivector_type nCoord(*globalMapG,dim);
  multivector_type nBound(*globalMapG,1);

  int indOwned = 0;
  for (int inode=0; inode<numNodes; inode++) {
    if (nodeIsOwned[inode]) {
      nCoord[0][indOwned]=nodeCoord(inode,0);
      nCoord[1][indOwned]=nodeCoord(inode,1);
      nBound[0][indOwned]=nodeOnBoundary(inode);
      indOwned++;
    }
  }
  EpetraExt::MultiVectorToMatrixMarketFile("coords.dat",nCoord,0,0,false);
  EpetraExt::MultiVectorToMatrixMarketFile("nodeOnBound.dat",nBound,0,0,false);

  // Put element to node mapping in multivector for output
  //driver_map_type   globalMapElem(numElemsGlobal, numElems, 0, Comm);
  driver_map_type   globalMapElem(numElemsGlobal, numElems, 0, Comm, NO);
  multivector_type elem2nodeMV(globalMapElem, numNodesPerElem);
  for (int ielem=0; ielem<numElems; ielem++) {
    for (int inode=0; inode<numNodesPerElem; inode++) {
      elem2nodeMV[inode][ielem]=globalNodeIds[elemToNode(ielem,inode)];
    }
  }
  EpetraExt::MultiVectorToMatrixMarketFile("elem2node.dat",elem2nodeMV,0,0,false);

  if(MyPID==0) {Time.ResetStartTime();}

#endif

  /**********************************************************************************/
  /************************** DIRICHLET BC SETUP ************************************/
  /**********************************************************************************/
  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Get Dirichlet boundary values")));

  int numBCNodes = 0;
  for (int inode = 0; inode < numNodes; inode++){
    if (nodeOnBoundary(inode) && P2_nodeIsOwned[inode]){
      numBCNodes++;
    }
  }


  // Vector for use in applying BCs
  multivector_type v(globalMapG,true);
  ArrayRCP<scalar_type> vdata = v.getDataNonConst(0);

  // Set v to boundary values on Dirichlet nodes
  //int * BCNodes = new int [numBCNodes];
  std::vector<int> BCNodes(numBCNodes);
  int indbc=0;
  int iOwned=0;
  for (int inode=0; inode<numNodes; inode++){
    if (P2_nodeIsOwned[inode]){
      if (nodeOnBoundary(inode)){
        BCNodes[indbc]=iOwned;
        indbc++;
        double x  = nodeCoord(inode, 0);
        double y  = nodeCoord(inode, 1);
        vdata[iOwned]=exactSolution(x, y);
      }
      iOwned++;
    }
  }

    
  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Global assembly")));

  /**********************************************************************************/
  /******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
  /**********************************************************************************/

  // Define desired workset size and count how many worksets there are on this processor's mesh block
  int desiredWorksetSize = numElems;                      // change to desired workset size!
  int numWorksets        = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksets*desiredWorksetSize < numElems) numWorksets += 1;

  if (MyPID == 0) {
    std::cout << "Building discretization matrix and right hand side... \n\n";
    std::cout << "\tDesired workset size:                 " << desiredWorksetSize <<"\n";
    std::cout << "\tNumber of worksets (per processor):   " << numWorksets <<"\n\n";
  }

  // Create P2 matrix and RHS
  CreateLinearSystem(numWorksets,
                     desiredWorksetSize,
                     elemToNode,
                     nodeCoord,
                     cubPoints,
                     cubWeights,
                     HGBGrads,
                     HGBValues,
                     P2_globalNodeIds,
                     P2_cellType,
                     StiffMatrix,
                     rhsVector,
                     msg
                     );

  /**********************************************************************************/
  /********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
  /**********************************************************************************/

  StiffMatrix.fillComplete();
  //writer_type::writeSparseFile("A0_p2_tpetra.m", rcpFromRef(StiffMatrix));
  //writer_type::writeDenseFile("rhs_tpetra.m", rhsVector);

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Getting cubature for auxiliary P1 mesh")));

  /////////////////////////////////////////////////////////////////////
  // Create P1 matrix and RHS to be used in preconditioner
  /////////////////////////////////////////////////////////////////////

  //Cubature

  // Get numerical integration points and weights
  int cubDegree_aux = 2; //TODO  This was 3 for P=2.  I think this should be 2 now .... Ask Chris.
  RCP<Cubature<double> >  myCub_aux = cubFactory.create(P1_cellType, cubDegree_aux);

  int cubDim_aux       = myCub_aux->getDimension();
  int numCubPoints_aux = myCub_aux->getNumPoints();

  FieldContainer<double> cubPoints_aux(numCubPoints_aux, cubDim_aux);
  FieldContainer<double> cubWeights_aux(numCubPoints_aux);
  myCub_aux->getCubature(cubPoints_aux, cubWeights_aux);
  
  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Getting basis for auxiliary P1 mesh")));

  //Basis

  // Define basis
  Basis_HGRAD_QUAD_C1_FEM<double, FieldContainer<double> > myHGradBasis_aux;
  int numFieldsG_aux = myHGradBasis_aux.getCardinality();
  FieldContainer<double> HGBValues_aux(numFieldsG_aux, numCubPoints_aux);
  FieldContainer<double> HGBGrads_aux(numFieldsG_aux, numCubPoints_aux, spaceDim);

  // Evaluate basis values and gradients at cubature points
  myHGradBasis_aux.getValues(HGBValues_aux, cubPoints_aux, OPERATOR_VALUE);
  myHGradBasis_aux.getValues(HGBGrads_aux, cubPoints_aux, OPERATOR_GRAD);

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Global assembly (auxiliary system)")));

  crs_matrix_type StiffMatrix_aux(globalMapG, 20*numFieldsG_aux);
  RCP<multivector_type> rhsVector_aux = rcp(new multivector_type(globalMapG,1));

  desiredWorksetSize = numElems_aux;
  numWorksets        = numElems_aux/desiredWorksetSize;
  if(numWorksets*desiredWorksetSize < numElems_aux) numWorksets += 1;
  CreateLinearSystem(numWorksets,
                     desiredWorksetSize,
                     aux_P1_elemToNode,
                     nodeCoord,
                     cubPoints_aux,
                     cubWeights_aux,
                     HGBGrads_aux,
                     HGBValues_aux,
                     P2_globalNodeIds,
                     P1_cellType,
                     StiffMatrix_aux,
                     rhsVector_aux,
                     msg
                     );

  /**********************************************************************************/
  /********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
  /**********************************************************************************/

  StiffMatrix_aux.fillComplete();
  //writer_type::writeSparseFile("A1_p1_aux_tpetra.m", rcpFromRef(StiffMatrix_aux));

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs")));

  // Generate P2-to-P1 identity coarsening (base mesh to auxiliary mesh).
  RCP<crs_matrix_type> P_identity, R_identity;
  if (inputSolverList.isParameter("aux P1")) {
    GenerateIdentityCoarsening_p2_to_p1(elemToNode, aux_P1_elemToNode, StiffMatrix_aux.getDomainMap(), StiffMatrix.getRangeMap(), P_identity, R_identity);
    //writer_type::writeSparseFile("IdentityP_tpetra.m", P_identity);
    //writer_type::writeSparseFile("IdentityR_tpetra.m", R_identity);
    inputSolverList.remove("aux P1"); //even though LevelWrap happily accepts this parameter
  }

/**********************************************************************************/
/******************************* ADJUST MATRIX DUE TO BC **************************/
/**********************************************************************************/

  //writer_type::writeDenseFile("rhs_tpetra_bcs.m", rhsVector);

  // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
  //  and add one to diagonal.
  Apply_Dirichlet_BCs(BCNodes,StiffMatrix,*femCoefficients,*rhsVector,v);

  // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
  //  and add one to diagonal.
  Apply_Dirichlet_BCs(BCNodes,StiffMatrix_aux,*rhsVector_aux,*rhsVector_aux,v);

  tm = Teuchos::null;

  //writer_type::writeSparseFile("A1_p1_aux_bcs_tpetra.m", rcpFromRef(StiffMatrix_aux));

  /**********************************************************************************/
  /*********************************** SOLVE ****************************************/
  /**********************************************************************************/

  // Which Solver?
  std::string solveType("cg");
  if(inputSolverList.isParameter("solver")) {
    solveType = inputSolverList.get<std::string>("solver");
  }

  // Run the solver
  std::string amgType("MueLu");

  ParameterList amgList;
  std::string seedType = inputSolverList.get("seed type","node");
  inputSolverList.remove("seed type");
  if (inputSolverList.isSublist("MueLu"))
    amgList = inputSolverList.sublist("MueLu");
  else
    amgList = inputSolverList;

  std::string lev0List = "level 0";
  if (amgList.isSublist(lev0List)) {
    std::cout << "found \"" << lev0List << "\" sublist" << std::endl;
    ParameterList &sl = amgList.sublist(lev0List);
    std::string smooType = sl.get<std::string>("smoother: type");
    if (smooType == "BLOCK RELAXATION" && sl.isParameter("smoother: params")) {
      ParameterList &ssl = sl.sublist("smoother: params");
      std::cout << "found \"smoother: params\" for block relaxation" << std::endl;
      int numLocalParts;
      std::cout << "setting \"partitioner: map\"" << std::endl;
      if (seedType == "node") {
        ssl.set("partitioner: map", nodeSeeds);
        numLocalParts = numNodeSeeds;
      } else if (seedType == "edge") {
        ssl.set("partitioner: map", edgeSeeds);
        numLocalParts = numEdgeSeeds;
      } else if (seedType == "cell") {
        ssl.set("partitioner: map", cellSeeds);
        numLocalParts = numCellSeeds;
      }
      std::cout << "setting \"partitioner: local parts\" = " << numLocalParts << std::endl;
      ssl.set("partitioner: local parts", numLocalParts);
    }
  }

  // /////////////////////////////////////////////////////////////////////// //

  // Shove coordinates in multivector for saving to file.
  multivector_type coordinates(StiffMatrix.getRowMap(),2);
  ArrayRCP<scalar_type> xdat = coordinates.getDataNonConst(0);
  ArrayRCP<scalar_type> ydat = coordinates.getDataNonConst(1);
  for (size_t i=0; i<P2_nodeCoordx.size(); ++i) {
    xdat[i] = P2_nodeCoordx[i];
    ydat[i] = P2_nodeCoordy[i];
  }
  xdat = Teuchos::null; ydat = Teuchos::null;
  //writer_type::writeDenseFile("coords.m", coordinates);

  multivector_type P1coordinates(P1_globalMap,2);
  xdat = P1coordinates.getDataNonConst(0);
  ydat = P1coordinates.getDataNonConst(1);
  for (size_t i=0; i<P1coordinates.getLocalLength(); i++) {
    xdat[i] = nodeCoordx[i];
    ydat[i] = nodeCoordy[i];
  }
  xdat = Teuchos::null; ydat = Teuchos::null;
  //writer_type::writeDenseFile("P1coords.m", P1coordinates);

  multivector_type nBound(globalMapG,true);

  RCP<multivector_type> exactNodalVals = rcp(new multivector_type(globalMapG,1));
  double TotalErrorResidual = 0.0;
  double TotalErrorExactSol = 0.0;

  // Get exact solution at nodes
  for (int i = 0; i<numNodes; i++) {
    if (P2_nodeIsOwned[i]){
      double x = nodeCoord(i,0);
      double y = nodeCoord(i,1);
      double exactu = exactSolution(x, y);

      int rowindex=P2_globalNodeIds[i];
      exactNodalVals->sumIntoGlobalValue(rowindex, 0, exactu);
    }
  }

  char probType[10] = "laplace";

  RCP<crs_matrix_type> interpolationMatrix, restrictionMatrix;
  if (P_identity != Teuchos::null) {
    amgList.set("user coarse matrix",(crs_matrix_type*)&StiffMatrix_aux);
    interpolationMatrix = P_identity;
    restrictionMatrix = R_identity;
  }
  if (P_linear != Teuchos::null) {
    interpolationMatrix = P_linear;
  }

  TestMultiLevelPreconditionerLaplace(probType, amgList,
                                      rcpFromRef(StiffMatrix),
                                      interpolationMatrix, restrictionMatrix, exactNodalVals,
                                      rhsVector,            femCoefficients,
                                      TotalErrorResidual,   TotalErrorExactSol,
                                      amgType,
				      solveType);

  /**********************************************************************************/
  /**************************** CALCULATE ERROR *************************************/
  /**********************************************************************************/

  tm.reset();
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Calculate error")));

  double L2err = 0.0;
  double L2errTot = 0.0;
  double H1err = 0.0;
  double H1errTot = 0.0;
  double Linferr = 0.0;
  double LinferrTot = 0.0;

#ifdef HAVE_MPI
  // Import solution onto current processor
  RCP<driver_map_type>  solnMap = rcp(new driver_map_type(static_cast<int>(numNodesGlobal), static_cast<int>(numNodesGlobal), 0, Comm));
  Tpetra::Import<local_ordinal_type, global_ordinal_type, NO> solnImporter(globalMapG,solnMap);
  multivector_type  uCoeff(solnMap,1);
  uCoeff.doImport(*femCoefficients, solnImporter, Tpetra::INSERT);
#endif

  // Define desired workset size
  desiredWorksetSize = numElems;
  int numWorksetsErr    = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksetsErr*desiredWorksetSize < numElems) numWorksetsErr += 1;

  // Get cubature points and weights for error calc (may be different from previous)
  Intrepid::DefaultCubatureFactory<double>  cubFactoryErr;
  int cubDegErr = 6;
  RCP<Intrepid::Cubature<double> > cellCubatureErr = cubFactoryErr.create(P2_cellType, cubDegErr);
  int cubDimErr       = cellCubatureErr->getDimension();
  int numCubPointsErr = cellCubatureErr->getNumPoints();
  Intrepid::FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
  Intrepid::FieldContainer<double> cubWeightsErr(numCubPointsErr);
  cellCubatureErr->getCubature(cubPointsErr, cubWeightsErr);

  // Evaluate basis values and gradients at cubature points
  Intrepid::FieldContainer<double> uhGVals(numFieldsG, numCubPointsErr);
  Intrepid::FieldContainer<double> uhGrads(numFieldsG, numCubPointsErr, spaceDim);
  myHGradBasis.getValues(uhGVals, cubPointsErr, Intrepid::OPERATOR_VALUE);
  myHGradBasis.getValues(uhGrads, cubPointsErr, Intrepid::OPERATOR_GRAD);

  // Loop over worksets
  for(int workset = 0; workset < numWorksetsErr; workset++){

    // compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // when numElems is not divisible by desiredWorksetSize, the last workset ends at numElems
    worksetEnd   = (worksetEnd <= numElems) ? worksetEnd : numElems;

    // now we know the actual workset size and can allocate the array for the cell nodes
    worksetSize  = worksetEnd - worksetBegin;
    Intrepid::FieldContainer<double> cellWorksetEr(worksetSize, P2_numNodesPerElem, spaceDim);
    Intrepid::FieldContainer<double> worksetApproxSolnCoef(worksetSize, P2_numNodesPerElem);

    // loop over cells to fill arrays with coordinates and discrete solution coefficient
    int cellCounter = 0;
    ArrayRCP<const scalar_type> uCoeffData = uCoeff.getData(0);
    ArrayRCP<const scalar_type> femCoeffData = femCoefficients->getData(0);
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      for (int node = 0; node < P2_numNodesPerElem; node++) {
        cellWorksetEr(cellCounter, node, 0) = nodeCoord( elemToNode(cell, node), 0);
        cellWorksetEr(cellCounter, node, 1) = nodeCoord( elemToNode(cell, node), 1);

        int rowIndex  = P2_globalNodeIds[elemToNode(cell, node)];
#ifdef HAVE_MPI
        worksetApproxSolnCoef(cellCounter, node) = uCoeffData[rowIndex];
#else
        worksetApproxSolnCoef(cellCounter, node) = femCoeffData[rowIndex];
#endif
      }

      cellCounter++;

    } // end cell loop

      // Containers for Jacobian
    Intrepid::FieldContainer<double> worksetJacobianE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
    Intrepid::FieldContainer<double> worksetJacobInvE(worksetSize, numCubPointsErr, spaceDim, spaceDim);
    Intrepid::FieldContainer<double> worksetJacobDetE(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetCubWeightsE(worksetSize, numCubPointsErr);

    // Containers for basis values and gradients in physical space
    Intrepid::FieldContainer<double> uhGValsTrans(worksetSize,numFieldsG, numCubPointsErr);
    Intrepid::FieldContainer<double> uhGradsTrans(worksetSize, numFieldsG, numCubPointsErr, spaceDim);

    // compute cell Jacobians, their inverses and their determinants
    IntrepidCTools::setJacobian(worksetJacobianE, cubPointsErr, cellWorksetEr, P2_cellType);
    IntrepidCTools::setJacobianInv(worksetJacobInvE, worksetJacobianE );
    IntrepidCTools::setJacobianDet(worksetJacobDetE, worksetJacobianE );

    // map cubature points to physical frame
    Intrepid::FieldContainer<double> worksetCubPoints(worksetSize, numCubPointsErr, cubDimErr);
    IntrepidCTools::mapToPhysicalFrame(worksetCubPoints, cubPointsErr, cellWorksetEr, P2_cellType);

    // evaluate exact solution and gradient at cubature points
    Intrepid::FieldContainer<double> worksetExactSoln(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetExactSolnGrad(worksetSize, numCubPointsErr, spaceDim);
    evaluateExactSolution(worksetExactSoln, worksetCubPoints);
    evaluateExactSolutionGrad(worksetExactSolnGrad, worksetCubPoints);

    // transform basis values to physical coordinates
    IntrepidFSTools::HGRADtransformVALUE<double>(uhGValsTrans, uhGVals);
    IntrepidFSTools::HGRADtransformGRAD<double>(uhGradsTrans, worksetJacobInvE, uhGrads);

    // compute weighted measure
    IntrepidFSTools::computeCellMeasure<double>(worksetCubWeightsE, worksetJacobDetE, cubWeightsErr);

    // evaluate the approximate solution and gradient at cubature points
    Intrepid::FieldContainer<double> worksetApproxSoln(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetApproxSolnGrad(worksetSize, numCubPointsErr, spaceDim);
    IntrepidFSTools::evaluate<double>(worksetApproxSoln, worksetApproxSolnCoef, uhGValsTrans);
    IntrepidFSTools::evaluate<double>(worksetApproxSolnGrad, worksetApproxSolnCoef, uhGradsTrans);

    // get difference between approximate and exact solutions
    Intrepid::FieldContainer<double> worksetDeltaSoln(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetDeltaSolnGrad(worksetSize, numCubPointsErr, spaceDim);
    IntrepidRSTools::subtract(worksetDeltaSoln, worksetApproxSoln, worksetExactSoln);
    IntrepidRSTools::subtract(worksetDeltaSolnGrad, worksetApproxSolnGrad, worksetExactSolnGrad);

    // take absolute values
    IntrepidRSTools::absval(worksetDeltaSoln);
    IntrepidRSTools::absval(worksetDeltaSolnGrad);
    // apply cubature weights to differences in values and grads for use in integration
    Intrepid::FieldContainer<double> worksetDeltaSolnWeighted(worksetSize, numCubPointsErr);
    Intrepid::FieldContainer<double> worksetDeltaSolnGradWeighted(worksetSize, numCubPointsErr, spaceDim);
    IntrepidFSTools::scalarMultiplyDataData<double>(worksetDeltaSolnWeighted,
                                                    worksetCubWeightsE, worksetDeltaSoln);
    IntrepidFSTools::scalarMultiplyDataData<double>(worksetDeltaSolnGradWeighted,
                                                    worksetCubWeightsE, worksetDeltaSolnGrad);

    // integrate to get errors on each element
    Intrepid::FieldContainer<double> worksetL2err(worksetSize);
    Intrepid::FieldContainer<double> worksetH1err(worksetSize);
    IntrepidFSTools::integrate<double>(worksetL2err, worksetDeltaSoln,
                                       worksetDeltaSolnWeighted, Intrepid::COMP_BLAS);
    IntrepidFSTools::integrate<double>(worksetH1err, worksetDeltaSolnGrad,
                                       worksetDeltaSolnGradWeighted, Intrepid::COMP_BLAS);

    // loop over cells to get errors for total workset
    cellCounter = 0;
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // loop over cubature points
      for(int nPt = 0; nPt < numCubPointsErr; nPt++){

        Linferr = std::max(Linferr, worksetDeltaSoln(cellCounter,nPt));

      }

      L2err += worksetL2err(cellCounter);
      H1err += worksetH1err(cellCounter);

      cellCounter++;

    } // end cell loop

  } // end loop over worksets

#ifdef HAVE_MPI
  // sum over all processors
  TC_sumAll(Comm,L2err,L2errTot);
  TC_sumAll(Comm,H1err,H1errTot);
  TC_maxAll(Comm,Linferr,LinferrTot);
#else
  L2errTot = L2err;
  H1errTot = H1err;
  LinferrTot = Linferr;
#endif


  if (MyPID == 0) {
    std::cout << "\n" << "L2 Error:  " << sqrt(L2errTot) <<"\n";
    std::cout << "H1 Error:  " << sqrt(H1errTot) <<"\n";
    std::cout << "LInf Error:  " << LinferrTot <<"\n\n";
  }

  const bool alwaysWriteLocal = false;
  const bool writeGlobalStats = true;
  const bool writeZeroTimers  = false;
  const bool ignoreZeroTimers = true;
  const std::string filter    = "";
  TimeMonitor::summarize(Comm.ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                         writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);


  // Cleanup
  for(long long b = 0; b < numElemBlk; b++){
    delete [] elmt_node_linkage[b];
    delete [] element_types[b];
  }
  delete [] block_ids;
  delete [] nodes_per_element;
  delete [] element_attributes;
  delete [] element_types;
  delete [] elmt_node_linkage;
  delete [] elements;
  delete [] P1_globalNodeIds;
  delete [] P1_nodeIsOwned;
  if(num_node_comm_maps > 0){
    delete [] node_comm_proc_ids;
    delete [] node_cmap_node_cnts;
    delete [] node_cmap_ids;
    for(long long i=0;i<num_node_comm_maps;i++){
      delete [] comm_node_ids[i];
      delete [] comm_node_proc_ids[i];
    }

    delete [] comm_node_ids;
    delete [] comm_node_proc_ids;
  }

  delete [] nodeCoordx;
  delete [] nodeCoordy;

  // delete mesh
  Delete_Pamgen_Mesh();

  return 0;

}
/**********************************************************************************/
/********************************* END MAIN ***************************************/
/**********************************************************************************/

/**********************************************************************************/
/************ USER DEFINED FUNCTIONS FOR EXACT SOLUTION ***************************/
/**********************************************************************************/

template<typename Scalar>
const Scalar exactSolution(const Scalar& x, const Scalar& y) {

  // Patch test: bi-linear function is in the FE space and should be recovered
  //  return 1. + x + y + x*y;

  // Analytic solution with homogeneous Dirichlet boundary data
  //return sin(M_PI*x)*sin(M_PI*y)**exp(x+y);

  // Analytic solution with inhomogeneous Dirichlet boundary data
  return exp(x + y )/(1. + x*y);
}


template<typename Scalar>
void materialTensor(Scalar material[][2], const Scalar& x, const Scalar& y) {

  material[0][0] = 1.;
  material[0][1] = 0.;
  //
  material[1][0] = 0.;
  material[1][1] = 1.;
}

/**********************************************************************************/
/************** AUXILIARY FUNCTIONS FROM EXACT SOLUTION ***************************/
/**********************************************************************************/

/************ Grad of Exact Solution ****************/
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[2], const Scalar& x, const Scalar& y) {

  // To enable derivatives of the gradient (i.e., 2nd derivatives of the exact solution) need 2 levels of fad types
  Sacado::Fad::SFad<Scalar,2> fad_x = x;
  Sacado::Fad::SFad<Scalar,2> fad_y = y;
  Sacado::Fad::SFad<Scalar,2> u;

  // Indicate the independent variables
  fad_x.diff(0,2);
  fad_y.diff(1,2);

  u = exactSolution(fad_x, fad_y);

  gradExact[0] = u.dx(0);
  gradExact[1] = u.dx(1);
}

/************ Source Term (RHS) ****************/
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y){

  Scalar u;
  Scalar grad_u[2];
  Scalar flux[2] = {0.0, 0.0};
  Scalar material[2][2];
  Scalar f = 0.;

  // Indicate the independent variables
  x.diff(0,2);
  y.diff(1,2);

  // Get exact solution and its gradient
  u = exactSolution(x, y);
  exactSolutionGrad(grad_u, x, y);

  // Get material tensor
  materialTensor<Scalar>(material, x, y);

  // Compute total flux = (A.grad u)
  for(int i = 0; i < 2; i++){

    // Add diffusive flux
    for(int j = 0; j < 2; j++){
      flux[i] += material[i][j]*grad_u[j];
    }
  }

  // Compute source term (right hand side): f = -div(A.grad u)
  f = -(flux[0].dx(0) + flux[1].dx(1));

  return f;
}

/**********************************************************************************/
/*************************** EVALUATION METHODS ***********************************/
/**********************************************************************************/

/************ Material Tensor ****************/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        matTensorValues,
                            const ArrayIn &   evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints        = evaluationPoints.dimension(1);
  int spaceDim         = evaluationPoints.dimension(2);

  double material[2][2];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);

      materialTensor<double>(material, x, y);

      for(int row = 0; row < spaceDim; row++){
        for(int col = 0; col < spaceDim; col++){
          matTensorValues(cell, pt, row, col) = material[row][col];
        }
      }
    }
  }
}
/************ Source Term (RHS) ****************/
template<class ArrayOut, class ArrayIn>
void evaluateSourceTerm(ArrayOut &       sourceTermValues,
                        const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      Sacado::Fad::SFad<double,2> x = evaluationPoints(cell, pt, 0);
      Sacado::Fad::SFad<double,2> y = evaluationPoints(cell, pt, 1);

      sourceTermValues(cell, pt) = sourceTerm<Sacado::Fad::SFad<double,2> >(x, y).val();
    }
  }
}

/************ Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolution(ArrayOut &       exactSolutionValues,
                           const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);

      exactSolutionValues(cell, pt) = exactSolution<double>(x, y);
    }
  }
}
/************ Grad of Exact Solution ****************/
template<class ArrayOut, class ArrayIn>
void evaluateExactSolutionGrad(ArrayOut &       exactSolutionGradValues,
                               const ArrayIn &  evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints = evaluationPoints.dimension(1);
  int spaceDim  = evaluationPoints.dimension(2);

  double gradient[2];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);

      exactSolutionGrad<double>(gradient, x, y);

      for(int row = 0; row < spaceDim; row++){
        exactSolutionGradValues(cell, pt, row) = gradient[row];
      }
    }
  }
}

/**********************************************************************************/
/******************************* TEST MueLu ***************************************/
/**********************************************************************************/

#include "ml_LevelWrap.h"

// Test MueLu
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
                                 ParameterList   & amgList,
                                 RCP<crs_matrix_type> const &A0,
                                 RCP<crs_matrix_type> const &P0,
                                 RCP<crs_matrix_type> const &R0,
                                 RCP<multivector_type> const & xexact,
                                 RCP<multivector_type> & b,
                                 RCP<multivector_type> & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol,
				 std::string &amgType,
				 std::string &solveType)
{
  int mypid = A0->getComm()->getRank();
  int maxIts = 100;
  double tol =1e-10;
  RCP<multivector_type> x = rcp(new multivector_type(xexact->getMap(),1));
  if (amgList.isParameter("solve Ae=0")) {
    b->scale(0.0);
    x->randomize();
  }

  linear_problem_type Problem(linear_problem_type(A0, x, b));
  RCP<multivector_type> lhs = Problem.getLHS();
  RCP<const multivector_type> rhs = Problem.getRHS();

  RCP<TimeMonitor>
  tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Total solve time")));

  // ==================== //
  // call MueLu and Belos //
  // ==================== //

  if (amgType != "ML" && amgType != "MueLu") {
    std::ostringstream errStr;
    errStr << "Error: \"amgType\" has value \"" << amgType << "\". Valid values are \"ML\" or \"MueLu\".";
    throw std::runtime_error(errStr.str());
  }

  int numIterations=0;
  if (amgType == "ML") {
    throw std::runtime_error("Error: ML does not support Tpetra objects");
  } else if (amgType == "MueLu") {
    // Turns a Tpetra::CrsMatrix into a MueLu::Matrix
    RCP<Xpetra::Matrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>> mueluA = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>(A0);
    // Multigrid Hierarchy
    crs_matrix_type * A1;
    RCP<Xpetra::Matrix <scalar_type, local_ordinal_type, global_ordinal_type> > xA1;
    bool userCoarseA = false;
    if (amgList.isParameter("user coarse matrix")) {
      A1=amgList.get<crs_matrix_type*>("user coarse matrix");
      xA1 = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>(rcpFromRef(*A1));
      amgList.remove("user coarse matrix");
      userCoarseA = true;
    }
    MueLu::ParameterListInterpreter<scalar_type,local_ordinal_type,global_ordinal_type> mueLuFactory(amgList);
    RCP<MueLu::Hierarchy<scalar_type,local_ordinal_type,global_ordinal_type> > H = mueLuFactory.CreateHierarchy();
    H->setVerbLevel(Teuchos::VERB_HIGH);
    H->GetLevel(0)->Set("A", mueluA);
    MueLu::FactoryManager<scalar_type,local_ordinal_type,global_ordinal_type,NO> M1, M2;
    if (P0 != Teuchos::null) {
      H->AddNewLevel();
      RCP<Xpetra::Matrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>> xP0 = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>(P0);
      H->GetLevel(1)->AddKeepFlag("P",MueLu::NoFactory::get(),MueLu::UserData);
      H->GetLevel(1)->Set("P", xP0);
      RCP<Xpetra::Matrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>> xR0 = MueLu::TpetraCrs_To_XpetraMatrix<scalar_type,local_ordinal_type,global_ordinal_type,NO>(R0);
      H->GetLevel(1)->AddKeepFlag("R",MueLu::NoFactory::get(),MueLu::UserData);
      H->GetLevel(1)->Set("R", xR0);
      if (mypid==0) std::cout << ">>>>>>>>>>>>>>>>>>>>> "
         << "adding P0 to AMG hierarchy" << std::endl;
    }
    RCP<multivector_type> nullspace;
    if (userCoarseA) {
      H->GetLevel(1)->AddKeepFlag("A",MueLu::NoFactory::get(),MueLu::UserData);
      H->GetLevel(1)->AddKeepFlag("Nullspace",MueLu::NoFactory::get(),MueLu::UserData);
      H->GetLevel(1)->Set("A", xA1);
      if (mypid==0) std::cout << ">>>>>>>>>>>>>>>>>>>>> "
         << "adding coarse auxiliary matrix to AMG hierarchy" << std::endl;
      // create nullspace for Level 1 (first coarse AMG level)
      nullspace = rcp(new multivector_type(xexact->getMap(),1));
      ArrayRCP<scalar_type> data = nullspace->getDataNonConst(0);
      for (int i=0; i<data.size(); ++i)
        data[i] = 1.0;
      data = Teuchos::null;
      RCP<Xpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,NO>> xnullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector<scalar_type,local_ordinal_type,global_ordinal_type,NO>(nullspace);
      H->GetLevel(1)->Set("Nullspace", xnullspace);
    }
    // Multigrid setup phase
    mueLuFactory.SetupHierarchy(*H);
    //copied directly from MueLu's blockReader.cpp

    H->IsPreconditioner(true);
    RCP<muelu_tpetra_operator> M = rcp(new muelu_tpetra_operator(H));

    Problem.setRightPrec(M);

    bool set = Problem.setProblem();
    if (set == false) {
      std::cout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // Belos parameter list
    ParameterList belosList;
    belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
    belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
    belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosList.set("Output Frequency",      1);
    belosList.set("Output Style",          Belos::Brief);
    //if (!scaleResidualHist)
    belosList.set("Implicit Residual Scaling", "None");

    // Create an iterative solver manager
    RCP< Belos::SolverManager<scalar_type, multivector_type, operator_type> > solver;
    if(solveType == "cg")
      solver = rcp(new Belos::PseudoBlockCGSolMgr<scalar_type, multivector_type, operator_type>(rcpFromRef(Problem), rcp(&belosList, false)));
    else if (solveType == "gmres") { 
      solver = rcp(new Belos::PseudoBlockGmresSolMgr<scalar_type, multivector_type, operator_type>(rcpFromRef(Problem), rcp(&belosList, false)));
    }
    else if (solveType == "fixed point" || solveType == "fixed-point") {
      solver = rcp(new Belos::FixedPointSolMgr<scalar_type, multivector_type, operator_type>(rcpFromRef(Problem), rcp(&belosList, false)));   
    }
    else {
      std::cout << "\nERROR:  Invalid solver '"<<solveType<<"'" << std::endl;
      return EXIT_FAILURE;
    }

    // Perform solve
    solver->solve();
    //writer_type::writeDenseFile("sol.m", x);
    numIterations = solver->getNumIters();
  }

  uh = lhs;

  ArrayRCP<const scalar_type> lhsdata = lhs->getData(0);
  ArrayRCP<const scalar_type> xexactdata = xexact->getData(0);
  // ======================================================= //
  // compute difference between exact solution and MueLu one //
  // ======================================================= //
  double d = 0.0, d_tot = 0.0 , s =0.0, s_tot=0.0;
  //for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i ) {
  for( size_t i=0 ; i<lhs->getMap()->getLocalNumElements() ; ++i ) {
    d += (lhsdata[i] - xexactdata[i]) * (lhsdata[i] - xexactdata[i]);
    s +=  xexactdata[i]* xexactdata[i];
  }
  lhsdata = Teuchos::null;
  xexactdata = Teuchos::null;

  TC_sumAll(A0->getComm(),d,d_tot);
  TC_sumAll(A0->getComm(),s,s_tot);

  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  Array<typename ScalarTraits::magnitudeType> Norm(1);
  vector_type Ax(rhs->getMap());
  A0->apply(*lhs, Ax);
  Ax.update(1.0, *rhs, -1.0);
  Ax.norm2(Norm);

  string msg = ProblemType;

  int numProc = A0->getComm()->getSize();
  if (A0->getComm()->getRank() == 0) {
    cout << msg << endl << "......Using " << numProc << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm[0] << endl;
    cout << msg << "......||x_exact - x||_2/||x_exact||_2 = " << sqrt(d_tot/s_tot) << endl;
  }

  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm[0];

  return(numIterations);

}

/*********************************************************************************************************/
/*********************************************************************************************************/
/*********************************************************************************************************/


void GenerateEdgeEnumeration(const FieldContainer<int> & elemToNode, const FieldContainer<double> & nodeCoord, FieldContainer<int> & elemToEdge, FieldContainer<int> & elemToEdgeOrient, FieldContainer<double> & edgeCoord) {
  // Not especially efficient, but effective... at least in serial!
  
  int numElems        = elemToNode.dimension(0);
  int numNodesperElem = elemToNode.dimension(1);
  int dim             = nodeCoord.dimension(1);

  // Sanity checks
  if(numNodesperElem !=4) throw std::runtime_error("Error: GenerateEdgeEnumeration only works on Quads!");
  if(elemToEdge.dimension(0)!=numElems || elemToEdge.dimension(1)!=4 || elemToEdge.dimension(0)!=elemToEdgeOrient.dimension(0) || elemToEdge.dimension(1)!=elemToEdgeOrient.dimension(1)) 
    throw std::runtime_error("Error: GenerateEdgeEnumeration array size mismatch");

  int edge_node0_id[4]={0,1,2,3};
  int edge_node1_id[4]={1,2,3,0};
  
  // Run over all the elements and start enumerating edges
  typedef std::map<std::pair<int,int>,int> en_map_type;
  en_map_type e2n;

  int num_edges=0;
  en_map_type edge_map;
  for(int i=0; i<numElems; i++) {

    // Generate edge pairs, based on the global orientation of "edges go from low ID to high ID"
    for(int j=0; j<4; j++) {
      int lo = std::min(elemToNode(i,edge_node0_id[j]),elemToNode(i,edge_node1_id[j]));
      int hi = std::max(elemToNode(i,edge_node0_id[j]),elemToNode(i,edge_node1_id[j]));

      std::pair<int,int> ep(lo,hi);

      int edge_id;
      en_map_type::iterator iter = edge_map.find(ep);      
      if(iter==edge_map.end()) {
        edge_map[ep] = num_edges;
        edge_id = num_edges;
        num_edges++;
      }
      else
        edge_id = (*iter).second;
      
      elemToEdge(i,j) = edge_id;
      elemToEdgeOrient(i,j) = (lo==elemToNode(i,edge_node0_id[j]))?1:-1;
    }
  }      

  // Fill out the edge centers (clobbering data if needed)
  edgeCoord.resize(num_edges,dim);
  for(int i=0; i<numElems; i++) {
    for(int j=0; j<4; j++) {      
      int n0 = elemToNode(i,edge_node0_id[j]);
      int n1 = elemToNode(i,edge_node1_id[j]);
      for(int k=0; k<dim; k++)
        edgeCoord(elemToEdge(i,j),k) = (nodeCoord(n0,k)+nodeCoord(n1,k))/2.0;
    }
  }
  
}

/*********************************************************************************************************/

void PromoteMesh(const int degree,
                 const FieldContainer<int>    & P1_elemToNode,
                 const FieldContainer<double> & P1_nodeCoord,
                 const FieldContainer<double> & P1_edgeCoord,
                 const FieldContainer<int>    & P1_elemToEdge,
                 const FieldContainer<int>    & P1_elemToEdgeOrient,
                 const FieldContainer<int>    & P1_nodeOnBoundary,
                 FieldContainer<int>          & P2_elemToNode,
                 FieldContainer<double>       & P2_nodeCoord,
                 FieldContainer<int>          & P2_nodeOnBoundary,
                 std::vector<int> & P2_edgeNodes,
                 std::vector<int> & P2_cellNodes) {

  int numElems           = P1_elemToNode.dimension(0);
  int P1_numNodesperElem = P1_elemToNode.dimension(1);
  int P1_numEdgesperElem = P1_elemToEdge.dimension(1);
  int P2_numNodesperElem = P2_elemToNode.dimension(1);
  int P1_numNodes        = P1_nodeCoord.dimension(0);
  int P1_numEdges        = P1_edgeCoord.dimension(0);
  int dim                = P1_nodeCoord.dimension(1);
  

  // Sanity checks
  if(P1_numNodesperElem !=4 || P2_numNodesperElem !=9 ) throw std::runtime_error("Error: GenerateEdgeEnumeration only works on Quads!");
  if(P1_elemToEdge.dimension(0)!=numElems || P1_elemToEdge.dimension(1)!=4 || P1_elemToEdge.dimension(0)!=P1_elemToEdgeOrient.dimension(0) || P1_elemToEdge.dimension(1)!=P1_elemToEdgeOrient.dimension(1) ||
     P2_elemToNode.dimension(0)!=numElems || P2_nodeCoord.dimension(0) != P1_numNodes+P1_numEdges+numElems)
    throw std::runtime_error("Error: GenerateEdgeEnumeration array size mismatch");

  /*Quad-9 Layout:   
    inode3 -- inode6 -- inode2    
    |                   |
    inode7    inode8    inode5
    |                   |       
    inode0 -- inode4 -- inode1
  */
  // Make the new el2node array
  for(int i=0; i<numElems; i++)  {    
    // P1 nodes
    for(int j=0; j<P1_numNodesperElem; j++) 
      P2_elemToNode(i,j) = P1_elemToNode(i,j);
    
    // P1 edges
    for(int j=0; j<P1_numEdgesperElem; j++)  {
      P2_elemToNode(i,P1_numNodesperElem+j) = P1_numNodes+P1_elemToEdge(i,j);
      P2_edgeNodes.push_back(P1_numNodes+P1_elemToEdge(i,j));
    }

    // P1 cells
    P2_elemToNode(i,P1_numNodesperElem+P1_numEdgesperElem) = P1_numNodes+P1_numEdges+i;
    P2_cellNodes.push_back(P1_numNodes+P1_numEdges+i);

  }

  // Make the new coordinates
  // P1 nodes
  for(int i=0; i<P1_numNodes; i++) 
    for(int j=0; j<dim; j++)
      P2_nodeCoord(i,j) = P1_nodeCoord(i,j);

  // P1 edges
  for(int i=0; i<P1_numEdges; i++) 
    for(int j=0; j<dim; j++)
      P2_nodeCoord(P1_numNodes+i,j) = P1_edgeCoord(i,j);

  // P1 cells
  for(int i=0; i<numElems; i++) 
    for(int j=0; j<dim; j++) {
      P2_nodeCoord(P1_numNodes+P1_numEdges+i,j) =0;
      for(int k=0; k<P1_numNodesperElem; k++)
        P2_nodeCoord(P1_numNodes+P1_numEdges+i,j) += P1_nodeCoord(P1_elemToNode(i,k),j) / P1_numNodesperElem;
    }
  


  // Update the boundary conditions
  int edge_node0_id[4]={0,1,2,3};
  int edge_node1_id[4]={1,2,3,0};

  // P1 nodes
  for(int i=0; i<P1_numNodes; i++) 
    P2_nodeOnBoundary(i) = P1_nodeOnBoundary(i);

  
  // P1 edges
  // Not all that efficient
  for(int i=0; i<numElems; i++) {
    for(int j=0; j<P1_numEdgesperElem; j++) {
      int P1_edge = P1_elemToEdge(i,j);
      int n0 = P1_elemToNode(i,edge_node0_id[j]);
      int n1 = P1_elemToNode(i,edge_node1_id[j]);

      if(P1_nodeOnBoundary(n0) && P2_nodeOnBoundary(n1)) 
        P2_nodeOnBoundary(P1_numNodes+P1_edge) = 1;
    }
  }
}

/*********************************************************************************************************/

void CreateP1MeshFromP2Mesh(
                            FieldContainer<int> const    & P2_elemToNode,
                            FieldContainer<int>          & P1_elemToNode)
{

  /*
    Main idea:
    For each P2 element
    Create four P1 elements
    For each of those P1 elements
    Create element to node map
  */

  int P2_numElems        = P2_elemToNode.dimension(0);

  int p1ElemCtr=0;
  for (int i=0; i<P2_numElems; ++i) {

    /*
      How "P1" elements are traversed in a P2 element

      inode3 -- inode6 -- inode2    

      |   3    |     2    |

      inode7 -- inode8 -- inode5

      |   0    |     1    |       

      inode0 -- inode4 -- inode1
    */

    P1_elemToNode(p1ElemCtr,0) = P2_elemToNode(i,0);
    P1_elemToNode(p1ElemCtr,1) = P2_elemToNode(i,4);
    P1_elemToNode(p1ElemCtr,2) = P2_elemToNode(i,8);
    P1_elemToNode(p1ElemCtr++,3) = P2_elemToNode(i,7);

    P1_elemToNode(p1ElemCtr,0) = P2_elemToNode(i,4);
    P1_elemToNode(p1ElemCtr,1) = P2_elemToNode(i,1);
    P1_elemToNode(p1ElemCtr,2) = P2_elemToNode(i,5);
    P1_elemToNode(p1ElemCtr++,3) = P2_elemToNode(i,8);

    P1_elemToNode(p1ElemCtr,0) = P2_elemToNode(i,8);
    P1_elemToNode(p1ElemCtr,1) = P2_elemToNode(i,5);
    P1_elemToNode(p1ElemCtr,2) = P2_elemToNode(i,2);
    P1_elemToNode(p1ElemCtr++,3) = P2_elemToNode(i,6);

    P1_elemToNode(p1ElemCtr,0) = P2_elemToNode(i,7);
    P1_elemToNode(p1ElemCtr,1) = P2_elemToNode(i,8);
    P1_elemToNode(p1ElemCtr,2) = P2_elemToNode(i,6);
    P1_elemToNode(p1ElemCtr++,3) = P2_elemToNode(i,3);

  }
} //CreateP1MeshFromP2Mesh

/*********************************************************************************************************/

void CreateLinearSystem(int numWorksets,
                        int desiredWorksetSize,
                        FieldContainer<int>    const &elemToNode,
                        FieldContainer<double> const &nodeCoord,
                        FieldContainer<double> const &cubPoints,
                        FieldContainer<double> const &cubWeights,
                        FieldContainer<double> const &HGBGrads,
                        FieldContainer<double> const &HGBValues,
                        std::vector<global_ordinal_type>       const &globalNodeIds,
                        shards::CellTopology const &cellType,
                        crs_matrix_type &StiffMatrix,
                        RCP<multivector_type> &rhsVector,
                        std::string &msg
                        )
{
  int numCubPoints = cubPoints.dimension(0);
  int cubDim = cubPoints.dimension(1);
  int spaceDim = cellType.getDimension();
  int numFieldsG = HGBGrads.dimension(0);
  long long numElems = elemToNode.dimension(0);
  int numNodesPerElem = cellType.getNodeCount();

  
    std::cout << "CreateLinearSystem:" << std::endl;
    std::cout << "     numCubPoints = " << numCubPoints << std::endl;
    std::cout << "     cubDim = " << cubDim << std::endl;
    std::cout << "     spaceDim = " << spaceDim << std::endl;
    std::cout << "     numFieldsG = " << numFieldsG << std::endl;
    std::cout << "     numElems = " << numElems << std::endl;
    std::cout << "     numNodesPerElem = " << numNodesPerElem << std::endl;
    std::cout << "     length(globalNodeIds) = " << globalNodeIds.size() << std::endl;
    std::cout << "     length(nodeCoord) = " << nodeCoord.dimension(0) << std::endl;
  

  if(nodeCoord.dimension(0) != Teuchos::as<int>(globalNodeIds.size())) {
    std::ostringstream errStr;
    errStr << "Error: CreateLinearSystem: length of coordinates != #nodes ("
           << nodeCoord.dimension(0) << "!=" << globalNodeIds.size() << ")";
    throw std::runtime_error(errStr.str());
  }

  for(int workset = 0; workset < numWorksets; workset++){

    // Compute cell numbers where the workset starts and ends
    int worksetSize  = 0;
    int worksetBegin = (workset + 0)*desiredWorksetSize;
    int worksetEnd   = (workset + 1)*desiredWorksetSize;

    // When numElems is not divisible by desiredWorksetSize, the last workset ends at numElems
    worksetEnd   = (worksetEnd <= numElems) ? worksetEnd : numElems;

    // Now we know the actual workset size and can allocate the array for the cell nodes
    worksetSize  = worksetEnd - worksetBegin;
    //std::cout << "     cellWorkset: begin=" << worksetBegin << ", end=" << worksetEnd
    //          << ", size=" << worksetSize << std::endl;
    FieldContainer<double> cellWorkset(worksetSize, numNodesPerElem, spaceDim);

    // Copy coordinates into cell workset
    int cellCounter = 0;
    for(int cell = worksetBegin; cell < worksetEnd; cell++){
      for (int node = 0; node < numNodesPerElem; node++) {
        cellWorkset(cellCounter, node, 0) = nodeCoord( elemToNode(cell, node), 0);
        cellWorkset(cellCounter, node, 1) = nodeCoord( elemToNode(cell, node), 1);
      }
      cellCounter++;
    }

    /**********************************************************************************/
    /*                                Allocate arrays                                 */
    /**********************************************************************************/

    // Containers for Jacobians, integration measure & cubature points in workset cells
    FieldContainer<double> worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobInv  (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetJacobDet  (worksetSize, numCubPoints);
    FieldContainer<double> worksetCubWeights(worksetSize, numCubPoints);
    FieldContainer<double> worksetCubPoints (worksetSize, numCubPoints, cubDim);

    // Containers for basis values transformed to workset cells and them multiplied by cubature weights
    FieldContainer<double> worksetHGBValues        (worksetSize, numFieldsG, numCubPoints);
    FieldContainer<double> worksetHGBValuesWeighted(worksetSize, numFieldsG, numCubPoints);
    FieldContainer<double> worksetHGBGrads         (worksetSize, numFieldsG, numCubPoints, spaceDim);
    FieldContainer<double> worksetHGBGradsWeighted (worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for diffusive & advective fluxes & non-conservative adv. term and reactive terms
    FieldContainer<double> worksetDiffusiveFlux(worksetSize, numFieldsG, numCubPoints, spaceDim);

    // Containers for material values and source term. Require user-defined functions
    FieldContainer<double> worksetMaterialVals (worksetSize, numCubPoints, spaceDim, spaceDim);
    FieldContainer<double> worksetSourceTerm   (worksetSize, numCubPoints);

    // Containers for workset contributions to the discretization matrix and the right hand side
    FieldContainer<double> worksetStiffMatrix (worksetSize, numFieldsG, numFieldsG);
    FieldContainer<double> worksetRHS         (worksetSize, numFieldsG);


    /**********************************************************************************/
    /*                                Calculate Jacobians                             */
    /**********************************************************************************/

    IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
    IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
    IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

    /**********************************************************************************/
    /*          Cubature Points to Physical Frame and Compute Data                    */
    /**********************************************************************************/

    // map cubature points to physical frame
    IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);

    // get A at cubature points
    evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

    // get source term at cubature points
    evaluateSourceTerm (worksetSourceTerm, worksetCubPoints);


    /**********************************************************************************/
    /*                         Compute Stiffness Matrix                               */
    /**********************************************************************************/

    // Transform basis gradients to physical frame:
    IntrepidFSTools::HGRADtransformGRAD<double>(worksetHGBGrads,                // DF^{-T}(grad u)
                                                worksetJacobInv,   HGBGrads);

    // Compute integration measure for workset cells:
    IntrepidFSTools::computeCellMeasure<double>(worksetCubWeights,              // Det(DF)*w = J*w
                                                worksetJacobDet, cubWeights);


    // Multiply transformed (workset) gradients with weighted measure
    IntrepidFSTools::multiplyMeasure<double>(worksetHGBGradsWeighted,           // DF^{-T}(grad u)*J*w
                                             worksetCubWeights, worksetHGBGrads);


    // Compute the diffusive flux:
    IntrepidFSTools::tensorMultiplyDataField<double>(worksetDiffusiveFlux,      //  A*(DF^{-T}(grad u)
                                                     worksetMaterialVals,
                                                     worksetHGBGrads);

    // Integrate to compute workset diffusion contribution to global matrix:
    IntrepidFSTools::integrate<double>(worksetStiffMatrix,                    //(DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
                                       worksetHGBGradsWeighted,
                                       worksetDiffusiveFlux, COMP_BLAS);

    /**********************************************************************************/
    /*                                   Compute RHS                                  */
    /**********************************************************************************/

    // Transform basis values to physical frame:
    IntrepidFSTools::HGRADtransformVALUE<double>(worksetHGBValues,              // clones basis values (u)
                                                 HGBValues);

    // Multiply transformed (workset) values with weighted measure
    IntrepidFSTools::multiplyMeasure<double>(worksetHGBValuesWeighted,          // (u)*J*w
                                             worksetCubWeights, worksetHGBValues);

    // Integrate worksetSourceTerm against weighted basis function set
    IntrepidFSTools::integrate<double>(worksetRHS,                             // f.(u)*J*w
                                       worksetSourceTerm,
                                       worksetHGBValuesWeighted,  COMP_BLAS);

    /**********************************************************************************/
    /*                         Assemble into Global Matrix and RHS                    */
    /**********************************************************************************/

    //"WORKSET CELL" loop: local cell ordinal is relative to numElems
    
    Array<scalar_type> vals1(1);
    Array<global_ordinal_type> cols1(1);
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // Compute cell ordinal relative to the current workset
      int worksetCellOrdinal = cell - worksetBegin;

      // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
      for (int cellRow = 0; cellRow < numFieldsG; cellRow++){

        int localRow  = elemToNode(cell, cellRow);
        int globalRow = globalNodeIds[localRow];
        double sourceTermContribution =  worksetRHS(worksetCellOrdinal, cellRow);

        //rhsVector.sumIntoGlobalValue(globalRow, 0, sourceTermContribution); 
        rhsVector->sumIntoGlobalValue(globalRow, 0, sourceTermContribution); 

        // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
        for (int cellCol = 0; cellCol < numFieldsG; cellCol++){

          int localCol  = elemToNode(cell, cellCol);
          int globalCol = globalNodeIds[localCol];
          double operatorMatrixContribution = worksetStiffMatrix(worksetCellOrdinal, cellRow, cellCol);
          cols1[0] = globalCol;
          vals1[0] = operatorMatrixContribution;
          StiffMatrix.insertGlobalValues(globalRow, cols1(), vals1());

        }// *** cell col loop ***
      }// *** cell row loop ***
    }// *** workset cell loop **

  }// *** workset loop ***

} //CreateLinearSystem

/*********************************************************************************************************/
void GenerateLinearCoarsening_p2_to_p1(const FieldContainer<int> & P2_elemToNode, RCP<driver_map_type> &P1_map, RCP<driver_map_type> & P2_map,RCP<crs_matrix_type>& P) {

  // Generate a P matrix that uses the linear coarsening from p2 to p1 on the base mesh.
  // This presumes that the P2 element has all of the P1 nodes numbered first
  // Resulting matrix is #P2nodes x #P1nodes
  double one     = 1.0;
  double half    = 0.5;    
  double quarter = 0.25;
  int edge_node0_id[4]={0,1,2,3};
  int edge_node1_id[4]={1,2,3,0};
  
  int Nelem=P2_elemToNode.dimension(0);
  if(P2_elemToNode.dimension(1) != 9) throw std::runtime_error("Unidentified element type");
  
  P = rcp(new crs_matrix_type(P2_map,0));

  Array<scalar_type> vals1(1);
  Array<global_ordinal_type> cols1(1);
  Array<scalar_type> vals2(2);
  Array<global_ordinal_type> cols2(2);
  Array<scalar_type> vals4(4);
  Array<global_ordinal_type> cols4(4);

  for(int i=0; i<Nelem; i++)
    for(int j=0; j<P2_elemToNode.dimension(1); j++) {
      int row = P2_elemToNode(i,j);

      if(j<4) {
        vals1[0] = one; cols1[0] = row;
        P->insertGlobalValues(row,cols1(),vals1());
      }
      else if (j>=4 && j<8){
        int col0 = P2_elemToNode(i,edge_node0_id[j-4]);
        int col1 = P2_elemToNode(i,edge_node1_id[j-4]);
        vals2[0] = half; vals2[1] = half;
        cols2[0] = col0; cols2[1] = col1;
        P->insertGlobalValues(row,cols2(),vals2());
      }
      else {
        cols4[0] = P2_elemToNode(i,0); cols4[1] = P2_elemToNode(i,1);
        cols4[2] = P2_elemToNode(i,2); cols4[3] = P2_elemToNode(i,3);
        vals4[0] = quarter; vals4[1] = quarter;
        vals4[2] = quarter; vals4[3] = quarter;
        P->insertGlobalValues(row,cols4(),vals4());
      }
    }
  P->fillComplete(P1_map,P2_map);
}

/*********************************************************************************************************/
void GenerateIdentityCoarsening_p2_to_p1(const FieldContainer<int> & P2_elemToNode,
                      const FieldContainer<int> & P1_elemToNode,
                      RCP<const driver_map_type> const & P1_map_aux, RCP<const driver_map_type> const &P2_map,
                      RCP<crs_matrix_type> & P,
                      RCP<crs_matrix_type> & R) {

  // Generate prolongator matrix that will be used to transfer from P1 on auxiliary mesh to  P2 on base mesh.
  // It's just the identity matrix.
  // By construction, the node numbering on the P1 auxiliary mesh is the same as on the P2 base mesh
  // (see CreateP1MeshFromP2Mesh).  The P2 map is the range map, the P1 auxiliary map is the domain map.

  double one = 1.0;
  P = rcp(new crs_matrix_type(P2_map,1));

  //We must keep track of the nodes already encountered.  Inserting more than once will cause
  //the values to be summed.  Using a hashtable would work -- we abuse std::map for this purpose.
  std::map<int,int> hashTable;
  int Nelem=P2_elemToNode.dimension(0);
  Array<scalar_type> vals1(1);
  vals1[0] = one;
  Array<global_ordinal_type> cols1(1);
  for(int i=0; i<Nelem; i++) {
    for(int j=0; j<P2_elemToNode.dimension(1); j++) {
      int row = P2_elemToNode(i,j);
      if (hashTable.find(row) == hashTable.end()) {
        //not found
        cols1[0] = row;
        P->insertGlobalValues(row,cols1(),vals1());
        hashTable[row] = 1;
      }
    }
  }
  P->fillComplete(P1_map_aux,P2_map);

  hashTable.clear();
  R = rcp(new crs_matrix_type(P1_map_aux,1));
  Nelem = P1_elemToNode.dimension(0);
  int nodesPerElem = P1_elemToNode.dimension(1);
  for(int i=0; i<Nelem; ++i) {
    for(int j=0; j<nodesPerElem; ++j) {
      int row = P1_elemToNode(i,j);
      if (hashTable.find(row) == hashTable.end()) {
        //not found
        cols1[0] = row;
        R->insertGlobalValues(row,cols1(),vals1());
        hashTable[row] = 1;
      }
    }
  }
  R->fillComplete(P2_map,P1_map_aux);
}

/*********************************************************************************************************/
void Apply_Dirichlet_BCs(std::vector<int> &BCNodes, crs_matrix_type & A, multivector_type & x, multivector_type & b,
                         const multivector_type & soln) {
  int N=(int)BCNodes.size();
  ArrayRCP<scalar_type> xdata = x.getDataNonConst(0);
  ArrayRCP<scalar_type> bdata = b.getDataNonConst(0);
  ArrayRCP<const scalar_type> solndata = soln.getData(0);

  A.resumeFill();

  for(int i=0; i<N; i++) {
    local_ordinal_type lrid = BCNodes[i];

    xdata[lrid]=bdata[lrid] = solndata[lrid];

    size_t numEntriesInRow = A.getNumEntriesInLocalRow(lrid);    typename crs_matrix_type::nonconst_local_inds_host_view_type cols("cols", numEntriesInRow);
    typename crs_matrix_type::nonconst_values_host_view_type vals("vals", numEntriesInRow);
    A.getLocalRowCopy(lrid, cols, vals, numEntriesInRow);

    for(size_t j=0; j<vals.extent(0); j++)
      vals(j) = (cols(j) == lrid) ? 1.0 : 0.0;

    A.replaceLocalValues(lrid, cols, vals);
  }

  A.fillComplete();

}
