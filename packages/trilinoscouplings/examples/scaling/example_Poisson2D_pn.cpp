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
    \li     Epetra to handle the global matrix and vector.
    \li     Isorropia to partition the matrix. (Optional)
    \li     ML to solve the linear system.


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
//#define DUMP_DATA_COORD
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
#include "Intrepid_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Epetra includes
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Import.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "pamgen_im_exodusII_l.h"
#include "pamgen_im_ne_nemesisI_l.h"
#include "pamgen_extras.h"

// AztecOO includes
#include "AztecOO.h"

// ML Includes
#include "ml_MultiLevelPreconditioner.h"

#ifdef HAVE_INTREPID_KOKKOS
#include "Sacado.hpp"
#else
// Sacado includes
#include "Sacado_No_Kokkos.hpp"
#endif

#include <sstream>


#if defined(HAVE_TRILINOSCOUPLINGS_MUELU)
#include "MueLu_IntrepidPCoarsenFactory.hpp"
#endif

using namespace std;
using namespace Intrepid;


// for debugging
#define zwrap(x) (std::abs(x)<1e-10 ? 0.0 : x)

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef Sacado::Fad::SFad<double,2>      Fad2; //# ind. vars fixed at 2
typedef Intrepid::FunctionSpaceTools     IntrepidFSTools;
typedef Intrepid::RealSpaceTools<double> IntrepidRSTools;
typedef Intrepid::CellTools<double>      IntrepidCTools;


struct fecomp{
  bool operator () ( topo_entity* x,  topo_entity*  y )const
  {
    if(x->sorted_local_node_ids < y->sorted_local_node_ids)return true;
    return false;
  }
};


int global_MyPID;//DEBUG

#if 0
 long long  *nodes_per_element   = new long long [numElemBlk];
  long long  *element_attributes  = new long long [numElemBlk];
  long long  *elements            = new long long [numElemBlk];
  char      **element_types       = new char * [numElemBlk];
  long long **elmt_node_linkage   = new long long * [numElemBlk];
#endif


//elements,nodes_per_element,elmt_node_linkage,numElemBlk,
struct PamgenMesh{ 
  int dim;

  /* Mesh connectivity info */
  long long numElemBlk;
  std::vector<long long> nodes_per_element; 
  std::vector<long long> element_attributes;
  std::vector<long long> elements;            
  std::vector<std::vector<char> > element_types;       
  std::vector<std::vector<long long> >elmt_node_linkage;   

  /* Communicator info */
  std::vector<long long> node_comm_proc_ids;
  std::vector<long long> node_cmap_node_cnts;
  std::vector<long long> node_cmap_ids;
  std::vector<long long*>  comm_node_ids;
  std::vector<std::vector<long long> > comm_node_proc_ids;
  long long num_node_comm_maps;

  PamgenMesh(int mydim):dim(mydim){}
  ~PamgenMesh() {
    for(int i=0;i<(int)comm_node_ids.size(); i++)
      delete [] comm_node_ids[i];
  }   


  void allocateConnectivity(long long my_numElemBlk) {
    numElemBlk = my_numElemBlk;
    nodes_per_element.resize(numElemBlk);
    element_attributes.resize(numElemBlk);
    elements.resize(numElemBlk);
    element_types.resize(numElemBlk);
    elmt_node_linkage.resize(numElemBlk);
  }



};



// forward declarations
void PromoteMesh_Pn_Kirby(const int degree, const EPointType & pointType,long long P1_globalNumNodes,long long P1_globalNumEdges,long long P1_globalNumElem,
                 const FieldContainer<int>    & P1_elemToNode,
                 const FieldContainer<double> & P1_nodeCoord,
                 const FieldContainer<double> & P1_edgeCoord,
                 const FieldContainer<int>    & P1_elemToEdge,
                 const FieldContainer<int>    & P1_elemToEdgeOrient,
                 const FieldContainer<int>    & P1_nodeOnBoundary,
		 const bool*                    P1_nodeIsOwned,
		 const std::vector<bool>      & P1_edgeIsOwned,
		 const std::vector<long long> & P1_globalNodeIds,
		 const std::vector<long long> & P1_globalEdgeIds,
		 const std::vector<long long> & P1_globalElementIds,
		 long long                    & Pn_globalNumNodes,
                 FieldContainer<int>          & Pn_elemToNode,
                 FieldContainer<double>       & Pn_nodeCoord,
		 FieldContainer<int>          & Pn_nodeOnBoundary,
		 std::vector<bool>            & Pn_nodeIsOwned,
		 std::vector<long long>       & Pn_globalNodeIds);


void PamgenEnumerateEdges(int numNodesPerElem, int numEdgesPerElem, int numNodesPerEdge,
			  FieldContainer<int> refEdgeToNode,const FieldContainer<double> & nodeCoord,
			  std::vector<long long> globalNodeIds,
			  PamgenMesh & mesh,
			  Epetra_Comm & Comm,
			  /*Output args */
			  std::vector<long long> & globalEdgeIds,
			  std::vector<bool> & edgeIsOwned,
			  std::vector<int> & ownedEdgeIds,
			  FieldContainer<int> &elemToEdge, FieldContainer<int> & elemToEdgeOrient, FieldContainer<int> &edgeToNode,FieldContainer<double> & edgeCoord,
			  long long & numEdgesGlobal);

void EnumerateElements(Epetra_Comm & Comm, int numMyElements, std::vector<long long> & globalElementIds);

void CreateP1MeshFromP2Mesh(const FieldContainer<int> & P2_elemToNode, FieldContainer<int> &aux_P1_elemToNode);

void CreateLinearSystem(int numWorkSets,
                        int desiredWorksetSize,
                        FieldContainer<int> const &elemToNode,
                        FieldContainer<double> const &nodeCoord,
                        FieldContainer<double> const &cubPoints,
                        FieldContainer<double> const &cubWeights,
                        Teuchos::RCP<Basis<double,FieldContainer<double> > > &myBasis_rcp,
                        FieldContainer<double> const &HGBGrads,
                        FieldContainer<double> const &HGBValues,
                        std::vector<long long> const &globalNodeIds,
                        Epetra_FECrsMatrix &StiffMatrix,
                        Epetra_FEVector &rhsVector,
                        std::string &msg,
                        Epetra_Time &Time
                        );


void GenerateLinearCoarsening_pn_kirby_to_p1(const int degree,const FieldContainer<int> & Pn_elemToNode, const std::vector<bool> & Pn_nodeIsOwned, Teuchos::RCP<Basis_HGRAD_QUAD_Cn_FEM<double,FieldContainer<double> > > &PnBasis_rcp,Teuchos::RCP<Basis<double,FieldContainer<double> > > &P1Basis_rcp,Epetra_Map & P1_map, Epetra_Map & Pn_map,Teuchos::RCP<Epetra_CrsMatrix>& P);

void GenerateIdentityCoarsening_pn_to_p1(const FieldContainer<int> & Pn_elemToNode,
                      Epetra_Map const & P1_map_aux, Epetra_Map const &Pn_map,
                      Teuchos::RCP<Epetra_CrsMatrix> & I);

void Apply_Dirichlet_BCs(std::vector<int> BCNodes, Epetra_FECrsMatrix & A, Epetra_MultiVector & x, Epetra_MultiVector & b, const Epetra_MultiVector & soln);

/**********************************************************************************/
/***************** FUNCTION DECLARATION FOR ML PRECONDITIONER *********************/
/**********************************************************************************/

/** \brief  ML Preconditioner

    \param  ProblemType        [in]    problem type
    \param  MLList             [in]    ML parameter list
    \param  A                  [in]    discrete operator matrix
    \param  xexact             [in]    exact solution
    \param  b                  [in]    right-hand-side vector
    \param  uh                 [out]   solution vector
    \param  TotalErrorResidual [out]   error residual
    \param  TotalErrorExactSol [out]   error in uh

*/
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
                                 Teuchos::ParameterList   & MLList,
                                 Teuchos::RCP<Epetra_CrsMatrix>   const & A,
                                 Teuchos::RCP<Epetra_CrsMatrix>   const & P,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol);

void CalculateError(const Epetra_FEVector & femCoefficients, const Epetra_Map &overlapMap,Epetra_Time & Time,  shards::CellTopology &Pn_cellType, Teuchos::RCP<Basis<double,FieldContainer<double> > > myHGradBasis_rcp, const FieldContainer<int> & elemToNode, const FieldContainer<double> & nodeCoord, int degree);

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

#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  rank=mpiSession.getRank();
  numProcs=mpiSession.getNProc();
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  //  if(numProcs!=1) {printf("Error: This test only currently works in serial\n");return 1;}
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  global_MyPID = MyPID;//DEBUG
  Epetra_Time Time(Comm);

  //Check number of arguments
  if (argc > 3) {
    std::cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
    std::cout <<"Usage:\n\n";
    std::cout <<"  ./TrilinosCouplings_examples_scaling_Example_Poisson.exe [meshfile.xml] [solver.xml]\n\n";
    std::cout <<"   meshfile.xml(optional) - xml file with description of Pamgen mesh\n\n";
    std::cout <<"   solver.xml(optional) - xml file with ML solver options\n\n";
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
  Teuchos::ParameterList inputMeshList;
  Teuchos::ParameterList inputSolverList;

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
  std::cout<<std::endl;

  if(xmlSolverInFileName.length()) {
    if (MyPID == 0)
      std::cout << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n\n";
    Teuchos::updateParametersFromXmlFile(xmlSolverInFileName, Teuchos::inoutArg(inputSolverList));
  } else if (MyPID == 0) std::cout << "Using default solver values ..." << std::endl;

  // Get pamgen mesh definition
  std::string meshInput = Teuchos::getParameter<std::string>(inputMeshList,"meshInput");
  int degree = inputMeshList.get("degree",2);

  // Get Isorropia and Zoltan parameters.
  Teuchos::ParameterList iso_paramlist = inputMeshList.sublist
    ("Isorropia Input") ;
  if (MyPID == 0) {
    std::cout << "Isorropia/Zoltan parameters" << std::endl;
    iso_paramlist.print(std::cout,2,true,true);
  }


  /**********************************************************************************/
  /***************************** GET CELL TOPOLOGY **********************************/
  /**********************************************************************************/

  // Get cell topology for base quad
  shards::CellTopology P1_cellType(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  shards::CellTopology Pn_cellType(shards::getCellTopologyData<shards::Quadrilateral<> >() );
  assert(P1_cellType.getDimension() == Pn_cellType.getDimension());

  // Get dimensions
  int P1_numNodesPerElem = P1_cellType.getNodeCount();
  int P1_numEdgesPerElem = P1_cellType.getEdgeCount();
  int spaceDim = P1_cellType.getDimension();
  int dim = 2;
  int P1_numNodesPerEdge=2;

  // Build reference element edge to node map
  FieldContainer<int> P1_refEdgeToNode(P1_numEdgesPerElem,P1_numNodesPerEdge);
  for (int i=0; i<P1_numEdgesPerElem; i++){
    P1_refEdgeToNode(i,0)=P1_cellType.getNodeMap(1, i, 0);
    P1_refEdgeToNode(i,1)=P1_cellType.getNodeMap(1, i, 1);
  }

  /**********************************************************************************/
  /******************************* GENERATE MESH ************************************/
  /**********************************************************************************/

  if (MyPID == 0) {
    std::cout << "Generating mesh ... \n\n";
  }

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  long long cr_result = Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);
  TrilinosCouplings::pamgen_error_check(std::cout,cr_result);
  PamgenMesh P1_mesh(dim);


  string msg("Poisson: ");
  if(MyPID == 0) {cout << msg << "Pamgen Setup     = " << Time.ElapsedTime() << endl; Time.ResetStartTime();}

  // Get mesh size info
  char title[100];
  long long numDim;
  long long numNodes;
  long long numElems;
  long long numNodeSets;
  long long numSideSets;
  int id = 0;
  long long my_numElemBlk;
  
  im_ex_get_init_l(id, title, &numDim, &numNodes,
                   &numElems, &my_numElemBlk, &numNodeSets,
                   &numSideSets);
  P1_mesh.allocateConnectivity(my_numElemBlk);


  long long P1_globalNumNodes;
  long long P1_globalNumElems;
  long long numElemBlkGlobal;
  long long numNodeSetsGlobal;
  long long numSideSetsGlobal;

  im_ne_get_init_global_l(id, &P1_globalNumNodes, &P1_globalNumElems,
                          &numElemBlkGlobal, &numNodeSetsGlobal,
                          &numSideSetsGlobal);

  // Print mesh information
  if (MyPID == 0){
    std::cout << " Number of Global Elements: " << P1_globalNumElems << " \n";
    std::cout << " Number of Global Nodes: "    << P1_globalNumNodes << " \n\n";
  }

  long long * block_ids = new long long [P1_mesh.numElemBlk];
  error += im_ex_get_elem_blk_ids_l(id, block_ids);

 

  for(long long i = 0; i < P1_mesh.numElemBlk; i ++){
    P1_mesh.element_types[i].resize(MAX_STR_LENGTH + 1);
    error += im_ex_get_elem_block_l(id,
                                    block_ids[i],
                                    P1_mesh.element_types[i].data(),
                                    &(P1_mesh.elements[i]),
                                    &(P1_mesh.nodes_per_element[i]),
                                    &(P1_mesh.element_attributes[i]));
  }

  /*connectivity*/
  for(long long b = 0; b < P1_mesh.numElemBlk; b++){
    P1_mesh.elmt_node_linkage[b].resize(P1_mesh.nodes_per_element[b]* P1_mesh.elements[b]);
    error += im_ex_get_elem_conn_l(id,block_ids[b],P1_mesh.elmt_node_linkage[b].data());
  }

  // Get node-element connectivity
  int telct = 0;
  FieldContainer<int> P1_elemToNode(numElems,P1_numNodesPerElem);
  for(long long b = 0; b <P1_mesh.numElemBlk; b++){
    for(long long el = 0; el < P1_mesh.elements[b]; el++){
      for (int j=0; j<P1_numNodesPerElem; j++) {
        P1_elemToNode(telct,j) = P1_mesh.elmt_node_linkage[b][el*P1_numNodesPerElem + j]-1;
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
  long long num_elem_comm_maps;
  im_ne_get_loadbal_param_l( id,
                             &num_internal_nodes,
                             &num_border_nodes,
                             &num_external_nodes,
                             &num_internal_elems,
                             &num_border_elems,
                             &P1_mesh.num_node_comm_maps,
                             &num_elem_comm_maps,
                             0/*unused*/ );

  if(P1_mesh.num_node_comm_maps > 0){
    P1_mesh.node_comm_proc_ids.resize(P1_mesh.num_node_comm_maps);
    P1_mesh.node_cmap_node_cnts.resize(P1_mesh.num_node_comm_maps);
    P1_mesh.node_cmap_ids.resize(P1_mesh.num_node_comm_maps);
    P1_mesh.comm_node_ids.resize(P1_mesh.num_node_comm_maps);
    P1_mesh.comm_node_proc_ids.resize(P1_mesh.num_node_comm_maps);

    long long *  elem_cmap_ids        = new long long [num_elem_comm_maps];
    long long *  elem_cmap_elem_cnts  = new long long [num_elem_comm_maps];


    if ( im_ne_get_cmap_params_l( id,
                                  P1_mesh.node_cmap_ids.data(),
                                  (long long*)P1_mesh.node_cmap_node_cnts.data(),
                                  elem_cmap_ids,
                                  (long long*)elem_cmap_elem_cnts,
                                  0/*not used proc_id*/ ) < 0 )++error;

    for(long long j = 0; j < P1_mesh.num_node_comm_maps; j++) {
      P1_mesh.comm_node_ids[j] = new long long [P1_mesh.node_cmap_node_cnts[j]];
      P1_mesh.comm_node_proc_ids[j].resize(P1_mesh.node_cmap_node_cnts[j]);
      if ( im_ne_get_node_cmap_l( id,
				  P1_mesh.node_cmap_ids[j],
                                  P1_mesh.comm_node_ids[j],
                                  P1_mesh.comm_node_proc_ids[j].data(),
                                  0/*not used proc_id*/ ) < 0 )++error;
      P1_mesh.node_comm_proc_ids[j] = P1_mesh.comm_node_proc_ids[j][0];
    }

    delete [] elem_cmap_ids;
    delete [] elem_cmap_elem_cnts;
  }

  if(!Comm.MyPID()) {cout << msg << "Mesh Queries     = " << Time.ElapsedTime() << endl; Time.ResetStartTime();}

  //Calculate global node ids
  std::vector<long long> P1_globalNodeIds(numNodes);
  bool*      P1_nodeIsOwned = new bool[numNodes];

  calc_global_node_ids(P1_globalNodeIds.data(),
                       P1_nodeIsOwned,
                       numNodes,
                       P1_mesh.num_node_comm_maps,
                       P1_mesh.node_cmap_node_cnts.data(),
                       P1_mesh.node_comm_proc_ids.data(),
                       P1_mesh.comm_node_ids.data(),
                       rank);


  if(MyPID==0) {cout << msg << "Global Node Nums = " << Time.ElapsedTime() << endl; Time.ResetStartTime();}

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

  if(MyPID ==0) {cout << msg << "Boundary Conds   = " << Time.ElapsedTime() << endl; Time.ResetStartTime();}


  /******************************************/
  // Enumerate Edges
  if(MyPID==0) printf("Using new edge enumeration\n");
  std::vector<long long> P1_globalEdgeIds;
  std::vector<bool> P1_edgeIsOwned;
  std::vector<int> P1_ownedEdgeIds;
  FieldContainer<int> P1_elemToEdge(numElems,P1_numEdgesPerElem);
  FieldContainer<int> P1_elemToEdgeOrient(numElems,P1_numEdgesPerElem);
  FieldContainer<int> P1_edgeToNode;
  FieldContainer<double> P1_edgeCoord;
  long long P1_globalNumEdges;

  PamgenEnumerateEdges(P1_numNodesPerElem,P1_numEdgesPerElem,P1_numNodesPerEdge,P1_refEdgeToNode,P1_nodeCoord,
		       P1_globalNodeIds,P1_mesh,
		       Comm, P1_globalEdgeIds,P1_edgeIsOwned,P1_ownedEdgeIds,P1_elemToEdge,P1_elemToEdgeOrient,P1_edgeToNode,P1_edgeCoord,P1_globalNumEdges);

#if 0
  {
    ostringstream ss;
    ss<<"***** Edge Enumeration ["<<Comm.MyPID()<<"] *****"<<endl;
    for (int i=0;i<numElems; i++) {
      for (int j=0; j<P1_numEdgesPerElem; j++)
	ss<<P1_elemToEdge(i,j)<<"("<<PG_globalEdgeIds[P1_elemToEdge(i,j)]<<")["<<P1_elemToEdgeOrient(i,j)<<"] ";
      ss<<endl;
    }
    ss<<"***** Edge2Node ["<<Comm.MyPID()<<"] *****"<<endl;
    for (int i=0; i<PG_edgeToNode.dimension(0); i++) {
      for (int j=0; j<PG_edgeToNode.dimension(1); j++)
	ss<<PG_edgeToNode(i,j)<<" ";
      ss<<endl;
    }
    printf("%s",ss.str().c_str());
    fflush(stdout);
  }
#endif

  /******************************************/
  // Enumerate Elements
  std::vector<long long> P1_globalElementIds;
  EnumerateElements(Comm,numElems,P1_globalElementIds);


  /******************************************/
  // Generate higher order mesh
  // The number of *my* Pn nodes will be related to the number of *my* nodes, edges (owned and ghosted) and elements (owned)
  int Pn_numNodes = numNodes + (degree-1)*P1_edgeCoord.dimension(0) + (degree-1)*(degree-1)*numElems;  // local
  int Pn_numNodesperElem = (degree+1)*(degree+1); // Quads
  FieldContainer<int>    elemToNode(numElems,Pn_numNodesperElem); 
  FieldContainer<double> nodeCoord(Pn_numNodes,dim);
  FieldContainer<int>   nodeOnBoundary(Pn_numNodes);

  //  printf("[%d] P1_numNodes = %d Pn_numNodes = %d Pn_numNodesperElem = %d degree = %d\n",Comm.MyPID(),(int)numNodes,Pn_numNodes,Pn_numNodesperElem,degree);

  std::vector<bool>Pn_nodeIsOwned(Pn_numNodes,false);
  std::vector<long long>Pn_globalNodeIds(Pn_numNodes);
  long long Pn_globalNumNodes=0;
  long long Pn_globalNumElems=P1_globalNumElems;

  PromoteMesh_Pn_Kirby(degree,POINTTYPE_EQUISPACED,P1_globalNumNodes,P1_globalNumEdges,P1_globalNumElems,
		       P1_elemToNode,P1_nodeCoord,P1_edgeCoord,P1_elemToEdge,P1_elemToEdgeOrient,P1_nodeOnBoundary,P1_nodeIsOwned,P1_edgeIsOwned,P1_globalNodeIds,P1_globalEdgeIds,P1_globalElementIds,
                       Pn_globalNumNodes,elemToNode, nodeCoord, nodeOnBoundary,Pn_nodeIsOwned,Pn_globalNodeIds);


  if(!MyPID) printf("P1 global (nodes,edges) = (%lld,%lld) Pn global nodes = %lld\n",P1_globalNumNodes,P1_globalNumEdges,Pn_globalNumNodes);

#if 0
  {
    ostringstream ss;
    ss<<"***** Node ids & ownership ["<<Comm.MyPID()<<"] *****"<<endl;
    for (int i=0; i<Pn_numNodes; i++)
      ss<<i<<"["<<Pn_globalNodeIds[i]<<","<<Pn_nodeIsOwned[i]<<"] ";
    ss<<endl;

    printf("%s",ss.str().c_str());
    fflush(stdout);
  }
#endif


  /******************************************/
  // Generate Auxillary Mesh
  int numElems_aux = numElems*4;  //4 P1 elements per Pn element in auxiliary mesh
  long long globalNumElems_aux = P1_globalNumElems*4;
  FieldContainer<int> aux_P1_elemToNode(numElems_aux,P1_numNodesPerElem); //4 P1 elements per Pn element
  CreateP1MeshFromP2Mesh(elemToNode, aux_P1_elemToNode);

  // Print mesh information
  if (MyPID == 0){
    std::cout << " Number of Pn Global Elements: " << Pn_globalNumElems << " \n";
    std::cout << " Number of Pn Global Nodes: " << Pn_globalNumNodes << " \n";
    std::cout << " Number of faux P1 Global Elements: " << globalNumElems_aux << " \n\n";
  }

  // Coordinates in single vector form
  std::vector<double> Pn_nodeCoordx(Pn_numNodes);
  std::vector<double> Pn_nodeCoordy(Pn_numNodes);
  for (int i=0; i<Pn_numNodes; i++) {
    Pn_nodeCoordx[i] = nodeCoord(i,0);
    Pn_nodeCoordy[i] = nodeCoord(i,1);
  }

  // Reset constants
  int P1_numNodes =numNodes;
  numNodes = Pn_numNodes;

  /**********************************************************************************/
  /********************************* GET CUBATURE ***********************************/
  /**********************************************************************************/

  // Get numerical integration points and weights
  DefaultCubatureFactory<double>  cubFactory;
  int cubDegree = 2*degree;
  Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(Pn_cellType, cubDegree);

  int cubDim       = myCub->getDimension();
  int numCubPoints = myCub->getNumPoints();

  FieldContainer<double> cubPoints(numCubPoints, cubDim);
  FieldContainer<double> cubWeights(numCubPoints);

  myCub->getCubature(cubPoints, cubWeights);

  if(MyPID==0) {std::cout << "Getting cubature                            "
                          << Time.ElapsedTime() << " sec \n"  ; Time.ResetStartTime();}

  /**********************************************************************************/
  /*********************************** GET BASIS ************************************/
  /**********************************************************************************/

  // Define basis
  Basis_HGRAD_QUAD_Cn_FEM<double, FieldContainer<double> > myHGradBasis(degree,POINTTYPE_EQUISPACED);
  Teuchos::RCP<Basis<double,FieldContainer<double> > > myHGradBasis_rcp(&myHGradBasis, false);
  Teuchos::RCP<Basis_HGRAD_QUAD_Cn_FEM<double,FieldContainer<double> > > myHGradBasisWithDofCoords_rcp(&myHGradBasis, false);

  // Auxillary p=1 basis
  Basis_HGRAD_QUAD_C1_FEM<double, FieldContainer<double> > myHGradBasis_aux;
  Teuchos::RCP<Basis<double,FieldContainer<double> > > myHGradBasis_aux_rcp(&myHGradBasis_aux, false);

  int numFieldsG = myHGradBasis.getCardinality();

  FieldContainer<double> HGBValues(numFieldsG, numCubPoints);
  FieldContainer<double> HGBGrads(numFieldsG, numCubPoints, spaceDim);

  // Evaluate basis values and gradients at cubature points
  myHGradBasis.getValues(HGBValues, cubPoints, OPERATOR_VALUE);
  myHGradBasis.getValues(HGBGrads, cubPoints, OPERATOR_GRAD);


  //#define OUTPUT_REFERENCE_FUNCTIONS
#ifdef OUTPUT_REFERENCE_FUNCTIONS
  // Evaluate basis values and gradients at DOF points 
  FieldContainer<double> DofCoords(numFieldsG,spaceDim);
  myHGradBasis.getDofCoords(DofCoords);
  FieldContainer<double> HGBValues_at_Dofs(numFieldsG,numFieldsG);
  myHGradBasis.getValues(HGBValues_at_Dofs, DofCoords, OPERATOR_VALUE);


  printf("*** Dof Coords ***\n");
  for(int j=0; j<numFieldsG; j++)
    printf("[%d] %10.2e %10.2e\n",j,DofCoords(j,0),DofCoords(j,1));
  
  printf("*** CubPoints & Weights ***\n");
  for(int j=0; j<numCubPoints; j++)
    printf("[%d] %10.2e %10.2e (%10.2e)\n",j,cubPoints(j,0),cubPoints(j,1),cubWeights(j));
  
  printf("*** HGBValues @ Dofs ***\n");
  for(int i=0; i<numFieldsG; i++) {
    printf("[%d] ",i);
    for(int j=0; j<numFieldsG; j++) {
      printf("%10.2e ",zwrap(HGBValues_at_Dofs(i,j)));      
    }
    printf("\n");
  }
  
  printf("*** HGBValues ***\n");
  for(int i=0; i<numFieldsG; i++) {
    printf("[%d] ",i);
    for(int j=0; j<numCubPoints; j++) {
      printf("%10.2e ",zwrap(HGBValues(i,j)));
      
    }
    printf("\n");
  }
  
  printf("*** HGBGrad ***\n");
  for(int i=0; i<numFieldsG; i++) {
    printf("[%d] ",i);
    for(int j=0; j<numCubPoints; j++) {
      printf("(%10.2e %10.2e)",zwrap(HGBGrads(i,j,0)),zwrap(HGBGrads(i,j,1)));
      printf("\n");
    }
  }    

  printf("*** Pn_nodeIsOwned ***\n");
  for(int i=0; i<Pn_numNodes; i++)
    printf("%d ",(int)Pn_nodeIsOwned[i]);
  printf("\n");

  printf("*** Pn_nodeOnBoundary ***\n");
  for(int i=0; i<Pn_numNodes; i++)
    printf("%d ",(int)nodeOnBoundary[i]);
  printf("\n");

#endif

  
  if(MyPID==0) {std::cout << "Getting basis                               "
                          << Time.ElapsedTime() << " sec \n"  ; Time.ResetStartTime();}

  /**********************************************************************************/
  /********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
  /**********************************************************************************/
  // Count owned nodes (Pn)
  int Pn_ownedNodes=0;
  for(int i=0;i<numNodes;i++)
    if(Pn_nodeIsOwned[i]) Pn_ownedNodes++;

  // Build a list of the OWNED global ids...
  // NOTE: We cast all of this down to ints for ease of use w/ epetra
  std::vector<int> Pn_ownedGIDs(Pn_ownedNodes);
  int oidx=0;
  for(int i=0;i<numNodes;i++)
    if(Pn_nodeIsOwned[i]){
      Pn_ownedGIDs[oidx]=(int)Pn_globalNodeIds[i];
      oidx++;
    }

  // Count owned nodes (P1)
  int P1_ownedNodes=0;
  for(int i=0;i<P1_numNodes;i++)
    if(P1_nodeIsOwned[i]) P1_ownedNodes++;
  std::vector<int> P1_ownedGIDs(P1_ownedNodes);
  oidx=0;
  for(int i=0;i<P1_numNodes;i++)
    if(P1_nodeIsOwned[i]){
      P1_ownedGIDs[oidx]=(int)P1_globalNodeIds[i];
      oidx++;
    }
       
  // Generate epetra map for nodes
  Epetra_Map globalMapG(-1,Pn_ownedNodes,&Pn_ownedGIDs[0],0,Comm);
    
  // Generate p1 map
  Epetra_Map P1_globalMap(-1,P1_ownedNodes,&P1_ownedGIDs[0],0,Comm);

  // Genetrate Pn-to-P1 coarsening.
  if (inputSolverList.isParameter("aux P1") && inputSolverList.isParameter("linear P1"))
    throw std::runtime_error("Can only specify \"aux P1\" or \"linear P1\", not both.");
  Teuchos::RCP<Epetra_CrsMatrix> P_linear;
  if (inputSolverList.isParameter("linear P1")) {
    printf("Generating Linear Pn-to-P1 coarsening...\n");
    GenerateLinearCoarsening_pn_kirby_to_p1(degree,elemToNode,Pn_nodeIsOwned,myHGradBasisWithDofCoords_rcp, myHGradBasis_aux_rcp,P1_globalMap,globalMapG,P_linear);
    inputSolverList.remove("linear P1"); //even though LevelWrap happily accepts this parameter
  }

  // Global arrays in Epetra format
  Epetra_FECrsMatrix StiffMatrix(Copy, globalMapG, 20*numFieldsG);
  Epetra_FEVector rhsVector(globalMapG);
  Epetra_FEVector femCoefficients(globalMapG);


  if(MyPID==0) {std::cout << msg << "Build global maps                           "
                          << Time.ElapsedTime() << " sec \n";  Time.ResetStartTime();}



#ifdef DUMP_DATA_COORD
  /**********************************************************************************/
  /**** PUT COORDINATES AND NODAL VALUES IN ARRAYS FOR OUTPUT (FOR PLOTTING ONLY) ***/
  /**********************************************************************************/

  // Put coordinates in multivector for output
  Epetra_MultiVector nCoord(globalMapG,dim);
  Epetra_MultiVector nBound(globalMapG,1);

  int indOwned = 0;
  for (int inode=0; inode<numNodes; inode++) {
    if (Pn_nodeIsOwned[inode]) {
      nCoord[0][indOwned]=nodeCoord(inode,0);
      nCoord[1][indOwned]=nodeCoord(inode,1);
      indOwned++;
    }
  }
  EpetraExt::MultiVectorToMatrixMarketFile("coords.dat",nCoord,0,0,false);

  if(MyPID==0) {Time.ResetStartTime();}

#endif

  /**********************************************************************************/
  /************************** DIRICHLET BC SETUP ************************************/
  /**********************************************************************************/

  int numBCNodes = 0;
  for (int inode = 0; inode < numNodes; inode++){
    if (nodeOnBoundary(inode) && Pn_nodeIsOwned[inode]){
      numBCNodes++;
    }
  }

  // Vector for use in applying BCs
  Epetra_MultiVector v(globalMapG,true);
  v.PutScalar(0.0);

  // Set v to boundary values on Dirichlet nodes
  std::vector<int> BCNodes(numBCNodes);
  int indbc=0;
  int iOwned=0;
  for (int inode=0; inode<numNodes; inode++){
    if (Pn_nodeIsOwned[inode]){
      if (nodeOnBoundary(inode)){
        BCNodes[indbc]=iOwned;
        indbc++;
        double x  = nodeCoord(inode, 0);
        double y  = nodeCoord(inode, 1);
        v[0][iOwned]=exactSolution(x, y);
      }
      iOwned++;
    }
  }
    
  if(MyPID==0) {std::cout << msg << "Get Dirichlet boundary values               "
                          << Time.ElapsedTime() << " sec \n\n"; Time.ResetStartTime();}
  
  /**********************************************************************************/
  /******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
  /**********************************************************************************/

  // Define desired workset size and count how many worksets there are on this processor's mesh block
  int desiredWorksetSize = numElems;                      // change to desired workset size!
  //int desiredWorksetSize = 100;                      // change to desired workset size!
  int numWorksets        = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksets*desiredWorksetSize < numElems) numWorksets += 1;

  if (MyPID == 0) {
    std::cout << "Building discretization matrix and right hand side... \n\n";
    std::cout << "\tDesired workset size:                 " << desiredWorksetSize <<"\n";
    std::cout << "\tNumber of worksets (per processor):   " << numWorksets <<"\n\n";
    Time.ResetStartTime();
  }

  // Create Pn matrix and RHS
  CreateLinearSystem(numWorksets,
                     desiredWorksetSize,
                     elemToNode,
                     nodeCoord,
                     cubPoints,
                     cubWeights,
                     myHGradBasis_rcp,
                     HGBGrads,
                     HGBValues,
                     Pn_globalNodeIds,
                     StiffMatrix,
                     rhsVector,
                     msg,
                     Time
                     );

  /**********************************************************************************/
  /********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
  /**********************************************************************************/

  StiffMatrix.GlobalAssemble();
  StiffMatrix.FillComplete();
  rhsVector.GlobalAssemble();

  if(MyPID==0) {std::cout << msg << "Global assembly                             "
                          << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

  /////////////////////////////////////////////////////////////////////
  // Create P1 matrix and RHS to be used in preconditioner
  /////////////////////////////////////////////////////////////////////
  //Cubature
  // Get numerical integration points and weights
  int cubDegree_aux = 2; //TODO  This was 3 for P=2.  I think this should be 2 now .... Ask Chris.
  Teuchos::RCP<Cubature<double> > myCub_aux = cubFactory.create(P1_cellType, cubDegree_aux);

  int cubDim_aux       = myCub_aux->getDimension();
  int numCubPoints_aux = myCub_aux->getNumPoints();

  FieldContainer<double> cubPoints_aux(numCubPoints_aux, cubDim_aux);
  FieldContainer<double> cubWeights_aux(numCubPoints_aux);
  myCub_aux->getCubature(cubPoints_aux, cubWeights_aux);

  if(MyPID==0) {std::cout << "Getting cubature for auxiliary P1 mesh      "
                          << Time.ElapsedTime() << " sec \n"  ; Time.ResetStartTime();}

  //Basis
  // Define basis
  int numFieldsG_aux = myHGradBasis_aux.getCardinality();
  FieldContainer<double> HGBValues_aux(numFieldsG_aux, numCubPoints_aux);
  FieldContainer<double> HGBGrads_aux(numFieldsG_aux, numCubPoints_aux, spaceDim);

  // Evaluate basis values and gradients at cubature points
  myHGradBasis_aux.getValues(HGBValues_aux, cubPoints_aux, OPERATOR_VALUE);
  myHGradBasis_aux.getValues(HGBGrads_aux, cubPoints_aux, OPERATOR_GRAD);

  if(MyPID==0) {std::cout << "Getting basis for auxiliary P1 mesh         "
                          << Time.ElapsedTime() << " sec \n"  ; Time.ResetStartTime();}

  Epetra_FECrsMatrix StiffMatrix_aux(Copy, globalMapG, 20*numFieldsG_aux);
  Epetra_FEVector rhsVector_aux(globalMapG);

  desiredWorksetSize = numElems_aux;
  numWorksets        = numElems_aux/desiredWorksetSize;
  if(numWorksets*desiredWorksetSize < numElems_aux) numWorksets += 1;
  CreateLinearSystem(numWorksets,
                     desiredWorksetSize,
                     aux_P1_elemToNode,
                     nodeCoord,
                     cubPoints_aux,
                     cubWeights_aux,
                     myHGradBasis_aux_rcp,
                     HGBGrads_aux,
                     HGBValues_aux,
                     Pn_globalNodeIds,
                     StiffMatrix_aux,
                     rhsVector_aux,
                     msg,
                     Time
                     );
  /**********************************************************************************/
  /********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
  /**********************************************************************************/

  StiffMatrix_aux.GlobalAssemble();
  StiffMatrix_aux.FillComplete();
  rhsVector_aux.GlobalAssemble();

  if(MyPID==0) {std::cout << msg << "Global assembly (auxiliary system)          "
                          << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}



  // Generate Pn-to-P1 identity coarsening (base mesh to auxiliary mesh).
  Teuchos::RCP<Epetra_CrsMatrix> P_identity;
  if (inputSolverList.isParameter("aux P1")) {
    printf("Generating Identity Pn-to-P1 coarsening...\n");
    GenerateIdentityCoarsening_pn_to_p1(elemToNode, StiffMatrix_aux.DomainMap(), StiffMatrix.RangeMap(), P_identity);
    inputSolverList.remove("aux P1"); //even though LevelWrap happily accepts this parameter
  }

  
/**********************************************************************************/
/******************************* ADJUST MATRIX DUE TO BC **************************/
/**********************************************************************************/

  // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
  //  and add one to diagonal.
  Apply_Dirichlet_BCs(BCNodes,StiffMatrix,femCoefficients,rhsVector,v);

  ///////////////////////////////////////////
  // Apply BCs to auxiliary P1 matrix
  ///////////////////////////////////////////
  // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
  //  and add one to diagonal.
  Apply_Dirichlet_BCs(BCNodes,StiffMatrix_aux,rhsVector_aux,rhsVector_aux,v);

  if(MyPID==0) {std::cout << msg << "Adjust global matrix and rhs due to BCs     " << Time.ElapsedTime()
                          << " sec \n"; Time.ResetStartTime();}


#ifdef DUMP_DATA
  // Dump matrices to disk
  EpetraExt::RowMatrixToMatlabFile("stiff_matrix.dat",StiffMatrix);
  EpetraExt::MultiVectorToMatrixMarketFile("rhs_vector.dat",rhsVector,0,0,false);
  EpetraExt::RowMatrixToMatlabFile("stiff_matrix_aux.dat",StiffMatrix_aux);
  EpetraExt::MultiVectorToMatrixMarketFile("rhs_vector_aux.dat",rhsVector_aux,0,0,false);
  //#undef DUMP_DATA
#endif

  /**********************************************************************************/
  /*********************************** SOLVE ****************************************/
  /**********************************************************************************/
  // Run the solver
  Teuchos::ParameterList MLList = inputSolverList;
  ML_Epetra::SetDefaults("SA", MLList, 0, 0, false);
  MLList.set("repartition: Zoltan dimensions",2);
  MLList.set("x-coordinates",Pn_nodeCoordx.data());
  MLList.set("y-coordinates",Pn_nodeCoordy.data());

  Epetra_FEVector exactNodalVals(globalMapG);
  double TotalErrorResidual = 0.0;
  double TotalErrorExactSol = 0.0;

  // Get exact solution at nodes
  for (int i = 0; i<numNodes; i++) {
    if (Pn_nodeIsOwned[i]){
      double x = nodeCoord(i,0);
      double y = nodeCoord(i,1);
      double exactu = exactSolution(x, y);

      int rowindex=Pn_globalNodeIds[i];
      exactNodalVals.SumIntoGlobalValues(1, &rowindex, &exactu);
    }
  }
  exactNodalVals.GlobalAssemble();

  char probType[10] = "laplace";

  Teuchos::RCP<Epetra_CrsMatrix> interpolationMatrix;
  if (P_identity != Teuchos::null) {
    MLList.set("user coarse matrix",(Epetra_CrsMatrix*)&StiffMatrix_aux);
    interpolationMatrix = P_identity;
  }
  if (P_linear != Teuchos::null) {
    interpolationMatrix = P_linear;
  }

  TestMultiLevelPreconditionerLaplace(probType, MLList,
                                      Teuchos::rcpFromRef(StiffMatrix), interpolationMatrix, exactNodalVals,
                                      rhsVector,            femCoefficients,
                                      TotalErrorResidual,   TotalErrorExactSol);


#ifdef DUMP_DATA
  // Dump matrices to disk
  EpetraExt::MultiVectorToMatrixMarketFile("lhs_vector.dat",femCoefficients,0,0,false);
  EpetraExt::MultiVectorToMatrixMarketFile("lhs_vector_exact.dat",exactNodalVals,0,0,false);
#endif

  // Build an overlap map for error calculation.  This is kind of like StiffMatrix's ColMap, but in a different order
  std::vector<int>Pn_globalNodeIds_int(Pn_numNodes);
  for(int i=0;i<Pn_numNodes; i++)
    Pn_globalNodeIds_int[i] = (int)Pn_globalNodeIds[i];
  Epetra_Map OverlapMap(-1,Pn_numNodes,Pn_globalNodeIds_int.data(),0,Comm);

  // Calculate Error
  CalculateError(femCoefficients,OverlapMap,Time,Pn_cellType,myHGradBasis_rcp,elemToNode,nodeCoord,degree);

  // Cleanup
  delete [] block_ids;

  delete [] nodeCoordx;
  delete [] nodeCoordy;
  delete [] P1_nodeIsOwned;

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
/******************************* TEST ML ******************************************/
/**********************************************************************************/

#include "ml_LevelWrap.h"

// Test ML
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
                                 Teuchos::ParameterList   & MLList,
                                 Teuchos::RCP<Epetra_CrsMatrix> const &A0,
                                 Teuchos::RCP<Epetra_CrsMatrix> const &P0,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol)
{
  Epetra_MultiVector x(xexact);
  x = uh;

  Epetra_LinearProblem Problem(&*A0,&x,&b);
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();

  Epetra_Time Time(A0->Comm());

  // =================== //
  // call ML and AztecOO //
  // =================== //

  AztecOO solver(Problem);
  Epetra_Operator *MLPrec;
  if (P0 == Teuchos::null)
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A0, MLList, true);
  else
    MLPrec = new ML_Epetra::LevelWrap(A0, P0, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);

  solver.Iterate(200, 1e-12);

  delete MLPrec;

  uh = *lhs;


#ifdef DUMP_DATA
    EpetraExt::MultiVectorToMatrixMarketFile("lhs_vector.dat",*lhs,0,0,false);
#undef DUMP_DATA
#endif

  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //
  double d = 0.0, d_tot = 0.0 , s =0.0, s_tot=0.0;
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i ) {
    d += ((*lhs)[0][i] - xexact[0][i]) * ((*lhs)[0][i] - xexact[0][i]);
    s +=  xexact[0][i]* xexact[0][i];
  }

  A0->Comm().SumAll(&d,&d_tot,1);
  A0->Comm().SumAll(&s,&s_tot,1);

  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A0->Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);

  string msg = ProblemType;

  if (A0->Comm().MyPID() == 0) {
    cout << msg << "......Using " << A0->Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2/||x_exact||_2 = " << sqrt(d_tot/s_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }

  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;

  return( solver.NumIters() );

}

/*********************************************************************************************************/
/*********************************************************************************************************/


/*********************************************************************************************************/
void PromoteMesh_Pn_Kirby(const int degree, const EPointType & pointType,long long P1_globalNumNodes, long long P1_globalNumEdges,long long P1_globalNumElem,
                 const FieldContainer<int>    & P1_elemToNode,
                 const FieldContainer<double> & P1_nodeCoord,
                 const FieldContainer<double> & P1_edgeCoord,
                 const FieldContainer<int>    & P1_elemToEdge,
                 const FieldContainer<int>    & P1_elemToEdgeOrient,
                 const FieldContainer<int>    & P1_nodeOnBoundary,
		 const bool * P1_nodeIsOwned,
		 const std::vector<bool>      & P1_edgeIsOwned,
		 const std::vector<long long> & P1_globalNodeIds,
		 const std::vector<long long> & P1_globalEdgeIds,
		 const std::vector<long long> & P1_globalElementIds,
		 long long                    & Pn_globalNumNodes,
                 FieldContainer<int>          & Pn_elemToNode,
                 FieldContainer<double>       & Pn_nodeCoord,
		 FieldContainer<int>          & Pn_nodeOnBoundary,
		 std::vector<bool>            & Pn_nodeIsOwned,
		 std::vector<long long>       & Pn_globalNodeIds
) {
//#define DEBUG_PROMOTE_MESH
  int numElems           = P1_elemToNode.dimension(0);
  int P1_numNodesperElem = P1_elemToNode.dimension(1);
  int P1_numEdgesperElem = P1_elemToEdge.dimension(1);
  int Pn_numNodesperElem = Pn_elemToNode.dimension(1);
  int P1_numNodes        = P1_nodeCoord.dimension(0);
  int P1_numEdges        = P1_edgeCoord.dimension(0);
  int dim                = P1_nodeCoord.dimension(1);


#ifdef DEBUG_PROMOTE_MESH
  int Pn_numNodes        = Pn_nodeCoord.dimension(0);
#endif

  int Pn_ExpectedNodesperElem = P1_numNodesperElem + (degree-1)*P1_numEdgesperElem + (degree-1)*(degree-1);
  int Pn_ExpectedNumNodes     = P1_numNodes + (degree-1)*P1_numEdges + (degree-1)*(degree-1)*numElems;  

  // Sanity checks
  if(P1_numNodesperElem !=4 || Pn_numNodesperElem !=Pn_ExpectedNodesperElem ) throw std::runtime_error("Error: PromoteMesh_Pn_Kirby only works on Quads!");
  if(P1_elemToEdge.dimension(0)!=numElems || P1_elemToEdge.dimension(1)!=4 || P1_elemToEdge.dimension(0)!=P1_elemToEdgeOrient.dimension(0) || P1_elemToEdge.dimension(1)!=P1_elemToEdgeOrient.dimension(1) ||
     Pn_elemToNode.dimension(0)!=numElems || Pn_nodeCoord.dimension(0) != Pn_ExpectedNumNodes)
    throw std::runtime_error("Error: PromoteMesh_Pn_Kirby array size mismatch");

  const CellTopologyData &cellTopoData = *shards::getCellTopologyData<shards::Quadrilateral<4> >();                                                  
  shards::CellTopology cellTopo(&cellTopoData);

  // Kirby elements are ordered strictly in lex ordering from the bottom left to the top right
  /*
    Kirby Quad-9 Layout (p=2)
    inode6 -- inode7 -- inode8    
    |                   |
    inode3    inode4    inode5
    |                   |       
    inode0 -- inode1 -- inode2


    Kirby Quad-16 Layout (p=3)
    inode12-- inode13-- inode14-- inode15
    |                             |       
    inode8    inode9    inode10   inode11
    |                             |       
    inode4    inode5    inode6    inode7
    |                             |       
    inode0 -- inode1 -- inode2 -- inode3
  */

  // We still number the global nodes in the exodus-style (node,edge,face) ordering, *by global orientation* but
  // the elem2node map has to account for Kirby-style elements
  int p1_node_in_pn[4] = {0,degree, (degree+1)*(degree+1)-1, degree*(degree+1)};
  int edge_node0_id[4]={0,1,2,3};
  int edge_node1_id[4]={1,2,3,0};

  // As you advance along each edge, we'll skip by the following amount to get to the next Kirby dof
  // edge_skip is in *local* orientation, so we'll multiply it by the elemToEdgeOrient value
  int edge_skip[4] = {1,degree+1,-1,-(degree+1)};
  int center_root = degree+2;

  //  printf("[%d] Global Edge Offset = %d Global Element Offset = %d\n",global_MyPID, P1_globalNumNodes, P1_globalNumNodes+P1_globalNumEdges*(degree-1));

  // Make the new el2node array
  for(int i=0; i<numElems; i++)  {    
    // P1 nodes
    for(int j=0; j<P1_numNodesperElem; j++) {
      int lid = P1_elemToNode(i,j);
      Pn_elemToNode(i,p1_node_in_pn[j]) = lid;
      Pn_nodeIsOwned[lid]   = P1_nodeIsOwned[lid];
      Pn_globalNodeIds[lid] = P1_globalNodeIds[lid];
    }  

    // P1 edges
    for(int j=0; j<P1_numEdgesperElem; j++){
      int orient   = P1_elemToEdgeOrient(i,j);
      int base_id  = (orient==1) ? p1_node_in_pn[edge_node0_id[j]] : p1_node_in_pn[edge_node1_id[j]];
      int skip     = orient*edge_skip[j];
      bool is_owned= P1_edgeIsOwned[P1_elemToEdge(i,j)];// new node is owned if the edge was
      for(int k=0; k<degree-1; k++) {
	int lid       = P1_numNodes+P1_elemToEdge(i,j)*(degree-1)+k;
	long long gid = P1_globalNumNodes + P1_globalEdgeIds[P1_elemToEdge(i,j)]*(degree-1)+k;
        Pn_elemToNode(i,base_id+(k+1)*skip) = lid;
	Pn_nodeIsOwned[lid]   = is_owned;
	Pn_globalNodeIds[lid] = gid;
      }
    }
 
    //    printf("GID Base = %d degree = %d\n", P1_globalNumNodes+P1_globalNumEdges*(degree-1) + P1_globalElementIds[i]*(degree-1)*(degree-1),degree);
    // P1 cells
    for(int j=0; j<degree-1; j++) 
      for(int k=0; k<degree-1; k++) {
	int offset    = j*(degree-1)+k;
	int lid       = P1_numNodes+P1_numEdges*(degree-1)+i*(degree-1)*(degree-1)+offset;
	long long gid = P1_globalNumNodes+P1_globalNumEdges*(degree-1) + P1_globalElementIds[i]*(degree-1)*(degree-1) + offset;
	//	printf("offset = %d lid = %d gid = %d P1_globalElementIds = %d\n",offset,lid,gid, P1_globalElementIds[i]);
	Pn_elemToNode(i,center_root+j*(degree+1)+k) = lid;
	Pn_nodeIsOwned[lid] = true; //elements are always owned
	Pn_globalNodeIds[lid] = gid;
      }

  }  

  Pn_globalNumNodes = P1_globalNumNodes+P1_globalNumEdges*(degree-1) + P1_globalNumElem*(degree-1)*(degree-1);

  // Get the p=n node locations
  FieldContainer<double> RefNodeCoords(Pn_numNodesperElem,dim);  
  FieldContainer<double> RefNodeCoords2(1,Pn_numNodesperElem,dim);  
  FieldContainer<double> PhysNodeCoords(1,Pn_numNodesperElem,dim);  
  Basis_HGRAD_QUAD_Cn_FEM<double, FieldContainer<double> > BasisPn(degree,pointType);
  BasisPn.getDofCoords(RefNodeCoords);

#ifdef DEBUG_PROMOTE_MESH
  printf("**** P=%d reference nodal coordinates ***\n",degree);
  for(int i=0; i<Pn_numNodesperElem; i++)
    printf("[%2d] %10.2f %10.2f\n",i,RefNodeCoords(i,0),RefNodeCoords(i,1));
#endif

  for(int i=0; i<Pn_numNodesperElem; i++)
    for(int k=0; k<dim; k++)
      RefNodeCoords2(0,i,k) = RefNodeCoords(i,k);

  // Make the new coordinates (inefficient, but correct)
 for(int i=0; i<numElems; i++)  {    
   // Get Element coords
   FieldContainer<double> my_p1_nodes(1,P1_numNodesperElem,dim);
   for(int j=0; j<P1_numNodesperElem; j++) {
     int jj = P1_elemToNode(i,j);
     for(int k=0; k<dim; k++)
       my_p1_nodes(0,j,k) = P1_nodeCoord(jj,k);
   }

   // Map the reference node location to physical   
   Intrepid::CellTools<double>::mapToPhysicalFrame(PhysNodeCoords,RefNodeCoords2,my_p1_nodes,cellTopo);

#ifdef DEBUG_PROMOTE_MESH
   printf("[%2d] PhysNodes  : ",i);
   for(int j=0; j<Pn_numNodesperElem; j++)
     printf("(%10.2f %10.2f) ",PhysNodeCoords(0,j,0),PhysNodeCoords(0,j,1));
   printf("\n");
#endif

   // Copy the physical node locations to the Pn_nodeCoord array
   for(int j=0; j<Pn_numNodesperElem; j++) {
     int jj = Pn_elemToNode(i,j);
     for(int k=0; k<dim; k++)
       Pn_nodeCoord(jj,k) = PhysNodeCoords(0,j,k);
   }

 }

#ifdef DEBUG_PROMOTE_MESH
 for(int i=0; i<numElems; i++)  { 
   printf("[%2d] Effective  : ",i);
   for(int j=0; j<Pn_numNodesperElem; j++)
     printf("(%10.2f %10.2f) ",Pn_nodeCoord(Pn_elemToNode(i,j),0),Pn_nodeCoord(Pn_elemToNode(i,j),1));
   printf("\n");
 }
#endif

  // Update the boundary conditions

  // P1 nodes
  for(int i=0; i<P1_numNodes; i++) 
    Pn_nodeOnBoundary(i) = P1_nodeOnBoundary(i);
  
  // P1 edges
  // Not all that efficient
  for(int i=0; i<numElems; i++) {
    for(int j=0; j<P1_numEdgesperElem; j++) {
      int n0 = P1_elemToNode(i,edge_node0_id[j]);
      int n1 = P1_elemToNode(i,edge_node1_id[j]);
      if(P1_nodeOnBoundary(n0) && P1_nodeOnBoundary(n1)) {
        for(int k=0; k<degree-1; k++) 
          Pn_nodeOnBoundary(P1_numNodes+P1_elemToEdge(i,j)*(degree-1)+k) =1;
      }
    }
  }

#ifdef DEBUG_PROMOTE_MESH
  // debug
  printf("**** P=%d elem2node  ***\n",degree);
  for(int i=0; i<numElems; i++) {
    printf("[%2d] ",i);
    for(int j=0; j<Pn_numNodesperElem; j++)
      printf("%2d ",Pn_elemToNode(i,j));
    printf("\n");
  }


  printf("**** P=%d nodes  ***\n",degree);
  for(int i=0; i<Pn_numNodes; i++) {
    printf("[%2d] %10.2e %10.2e (%d)\n",i,Pn_nodeCoord(i,0),Pn_nodeCoord(i,1),(int)Pn_nodeOnBoundary(i));   
  }
  
  printf("**** P=1 nodes  ***\n");
  for(int i=0; i<P1_numNodes; i++) {
    printf("[%2d] %10.2e %10.2e (%d)\n",i,P1_nodeCoord(i,0),P1_nodeCoord(i,1),(int)P1_nodeOnBoundary(i));   
  }
#endif

}


/*********************************************************************************************************/

void CreateP1MeshFromP2Mesh(FieldContainer<int> const    & P2_elemToNode,                           
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
                        Teuchos::RCP<Basis<double,FieldContainer<double> > >&myBasis_rcp,
                        FieldContainer<double> const &HGBGrads,
                        FieldContainer<double> const &HGBValues,
                        std::vector<long long> const &globalNodeIds,
                        Epetra_FECrsMatrix &StiffMatrix,
                        Epetra_FEVector &rhsVector,
                        std::string &msg0,
                        Epetra_Time &Time
                        )
{
  int MyPID = StiffMatrix.Comm().MyPID();
  int numCubPoints = cubPoints.dimension(0);
  int cubDim = cubPoints.dimension(1);
  int spaceDim = nodeCoord.dimension(1);
  int numFieldsG = HGBGrads.dimension(0);
  long long numElems = elemToNode.dimension(0);
  int numNodesPerElem = elemToNode.dimension(1);


  
  if(!global_MyPID) {
    std::cout << "CreateLinearSystem:" << std::endl;
    std::cout << "     numCubPoints = " << numCubPoints << std::endl;
    std::cout << "     cubDim = " << cubDim << std::endl;
    std::cout << "     spaceDim = " << spaceDim << std::endl;
    std::cout << "     numFieldsG = " << numFieldsG << std::endl;
    //    std::cout << "     numElems = " << numElems << std::endl;
    std::cout << "     numNodesPerElem = " << numNodesPerElem << std::endl;
    //    std::cout << "     length(globalNodeIds) = " << globalNodeIds.size() << std::endl;
    //    std::cout << "     length(nodeCoord) = " << nodeCoord.dimension(0) << std::endl;
  }
  std::string msg = "     " + msg0;

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

    if(MyPID==0) {std::cout << msg << "Allocate arrays                             "
                            << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


    /**********************************************************************************/
    /*                                Calculate Jacobians                             */
    /**********************************************************************************/
    // Due to the way setJacobian works, we need to use the basis form here
    //    IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
    IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, myBasis_rcp);
    IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
    IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

    //#define DEBUG_OUTPUT
#ifdef DEBUG_OUTPUT
      printf("*** cubPoints ***\n");
      for(int j=0; j<numCubPoints; j++)
        printf("(%d) [%10.2e %10.2e]\n",j,cubPoints(j,0),cubPoints(j,1));

      printf("*** cellWorkset ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numNodesPerElem; j++)
          printf("(%d,%d) [%10.2e %10.2e]\n",i,j,cellWorkset(i,j,0),cellWorkset(i,j,1));
      
      printf("*** Jacobian ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numCubPoints; j++)
          printf("(%d,%d) [%10.2e %10.2e; %10.2e %10.2e]\n",i,j,zwrap(worksetJacobian(i,j,0,0)),zwrap(worksetJacobian(i,j,0,1)),zwrap(worksetJacobian(i,j,1,0)),zwrap(worksetJacobian(i,j,1,1)));

      printf("*** det J ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numCubPoints; j++)
          printf("(%d,%d) %10.2e\n",i,j,zwrap(worksetJacobDet(i,j)));

#endif
    

    if(MyPID==0) {std::cout << msg << "Calculate Jacobians                         "
                            << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

    /**********************************************************************************/
    /*          Cubature Points to Physical Frame and Compute Data                    */
    /**********************************************************************************/

    // map cubature points to physical frame
       IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, myBasis_rcp);

    // get A at cubature points
    evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

    // get source term at cubature points
    evaluateSourceTerm (worksetSourceTerm, worksetCubPoints);

    if(MyPID==0) {std::cout << msg << "Map to physical frame and get source term   "
                            << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


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

    if(MyPID==0) {std::cout << msg << "Compute stiffness matrix                    "
                            << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}


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
#ifdef DEBUG_OUTPUT      
      printf("*** worksetSourceTerm ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numCubPoints; j++)
          printf("(%d,%d) %10.2e\n",i,j,worksetSourceTerm(i,j));


      printf("*** worksetCubWeights ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numCubPoints; j++)
          printf("(%d,%d) %10.2e\n",i,j,worksetCubWeights(i,j));


      printf("*** worksetHGBValues ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numFieldsG; j++) {
          printf("(%d,%d) ",i,j);
          for(int k=0; k<numCubPoints; k++) {
            printf("%10.2e ",zwrap(worksetHGBValues(i,j,k)));       
          }
          printf("\n");
        }
      
      printf("*** worksetHGBValuesWeighted ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numFieldsG; j++) {
          printf("(%d,%d) ",i,j);
          for(int k=0; k<numCubPoints; k++) {
            printf("%10.2e ",zwrap(worksetHGBValues(i,j,k)));       
          }
          printf("\n");
        }     

      printf("*** worksetRHS ***\n");
      for(int i=0; i<worksetSize; i++)
        for(int j=0; j<numFieldsG; j++)
          printf("(%d,%d) %10.2e\n",i,j,worksetRHS(i,j));

      printf("*** worksetStiffMatrix ***\n");
      for(int i=0; i<worksetSize; i++) {
        printf("(%d) [ ",i);
        for(int j=0; j<numFieldsG; j++) {
          for(int k=0; k<numFieldsG; k++) {
            printf("%10.2e ",zwrap(worksetStiffMatrix(i,j,k)));     
          }
          if(j!=numFieldsG-1) printf(";");
        }
        printf("]\n");
      }     



#endif

    if(MyPID==0) {std::cout << msg << "Compute right-hand side                     "
                            << Time.ElapsedTime() << " sec \n"; Time.ResetStartTime();}

    /**********************************************************************************/
    /*                         Assemble into Global Matrix and RHS                    */
    /**********************************************************************************/

    //"WORKSET CELL" loop: local cell ordinal is relative to numElems
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // Compute cell ordinal relative to the current workset
      int worksetCellOrdinal = cell - worksetBegin;

      // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
      for (int cellRow = 0; cellRow < numFieldsG; cellRow++){

        int localRow  = elemToNode(cell, cellRow);
        int globalRow = (int) globalNodeIds[localRow];
        double sourceTermContribution =  worksetRHS(worksetCellOrdinal, cellRow);

        rhsVector.SumIntoGlobalValues(1,&globalRow,&sourceTermContribution);

        // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
        for (int cellCol = 0; cellCol < numFieldsG; cellCol++){

          int localCol  = elemToNode(cell, cellCol);
          int globalCol = (int) globalNodeIds[localCol];
          double operatorMatrixContribution = worksetStiffMatrix(worksetCellOrdinal, cellRow, cellCol);
          StiffMatrix.InsertGlobalValues(1, &globalRow, 1, &globalCol, &operatorMatrixContribution);

        }// *** cell col loop ***
      }// *** cell row loop ***
    }// *** workset cell loop **

  }// *** workset loop ***

} //CreateLinearSystem

/*********************************************************************************************************/
void GenerateLinearCoarsening_pn_kirby_to_p1(const int degree,const FieldContainer<int> & Pn_elemToNode, const std::vector<bool> & Pn_nodeIsOwned,Teuchos::RCP<Basis_HGRAD_QUAD_Cn_FEM<double,FieldContainer<double> > > &PnBasis_rcp,Teuchos::RCP<Basis<double,FieldContainer<double> > > &P1Basis_rcp, Epetra_Map & P1_map, Epetra_Map & Pn_map,Teuchos::RCP<Epetra_CrsMatrix>& P) {

  // Sanity checks
  assert(Pn_elemToNode.dimension(1) == PnBasis_rcp->getCardinality());
  assert(P1Basis_rcp->getCardinality() == 4);

  // Generate a P matrix that uses the linear coarsening from pn to p1 on the base mesh.
  // This presumes that the Pn element is number according to the Kirby convention (aka straight across, bottom to top)
  // Resulting matrix is #Pnnodes x #P1nodes
  int p1_node_in_pn[4] = {0,degree, (degree+1)*(degree+1)-1, degree*(degree+1)};

  // Get the reference coordinates for the Pn element  
  int numFieldsPn = PnBasis_rcp->getCardinality();
  int spaceDim    = PnBasis_rcp->getBaseCellTopology().getDimension();
  FieldContainer<double> PnDofCoords(numFieldsPn,spaceDim);
  PnBasis_rcp->getDofCoords(PnDofCoords);

  // Evaluate the linear basis functions at the Pn nodes
  int numFieldsP1 = P1Basis_rcp->getCardinality();
  FieldContainer<double> P1Values_at_PnDofs(numFieldsP1,numFieldsPn);
  P1Basis_rcp->getValues(P1Values_at_PnDofs, PnDofCoords, OPERATOR_VALUE);

  // Generate P
  int Nelem=Pn_elemToNode.dimension(0);  
  P = Teuchos::rcp(new Epetra_CrsMatrix(Copy,Pn_map,0));

  // Assemble
  std::vector<int> touched(P1_map.NumMyElements(),0);

  for(int i=0; i<Nelem; i++) {
    for(int j=0; j<numFieldsPn; j++) {
      int row_lid = Pn_elemToNode(i,j);
      int row_gid = Pn_map.GID(row_lid);
      for(int k=0; k<numFieldsP1; k++) {
        int col_lid = Pn_elemToNode(i,p1_node_in_pn[k]);
	int col_gid = P1_map.GID(col_lid);
        double val = P1Values_at_PnDofs(k,j);
	if(Pn_nodeIsOwned[row_lid] && !touched[row_lid]) {
	  P->InsertGlobalValues(row_gid,1,&val,&col_gid);
	  touched[row_lid]=1;
	}
      }
    }
  }
  P->FillComplete(P1_map,Pn_map);
}

/*********************************************************************************************************/
void GenerateIdentityCoarsening_pn_to_p1(const FieldContainer<int> & Pn_elemToNode,
                      Epetra_Map const & P1_map_aux, Epetra_Map const &Pn_map,
                      Teuchos::RCP<Epetra_CrsMatrix> & P) {

  // Generate prolongator matrix that will be used to transfer from P1 on auxiliary mesh to Pn on base mesh.
  // It's just the identity matrix.
  // By construction, the node numbering on the P1 auxiliary mesh is the same as on the Pn base mesh
  // (see CreateP1MeshFromPnMesh).  The P2 map is the range map, the P1 auxiliary map is the domain map.
  
  double one = 1.0;
  P = Teuchos::rcp(new Epetra_CrsMatrix(Copy,Pn_map,0));
 
  // Identity matrices are easy - exploit the fact that the maps are identical
  for(int i=0; i<Pn_map.NumMyElements(); i++) {
    int row_gid = Pn_map.GID(i);
    P->InsertGlobalValues(row_gid,1,&one,&row_gid);
  }
  P->FillComplete(P1_map_aux,Pn_map);
}





/*********************************************************************************************************/
void Apply_Dirichlet_BCs(std::vector<int> BCNodes, Epetra_FECrsMatrix & A, Epetra_MultiVector & x, Epetra_MultiVector & b, const Epetra_MultiVector & soln) {
  int N=(int)BCNodes.size();
  for(int i=0; i<N; i++) {
    int lrid = BCNodes[i];
    
    int NumEntries, *Indices;
    double *Values;

    x[0][lrid]=b[0][lrid] = soln[0][lrid];

    A.ExtractMyRowView(lrid,NumEntries,Values,Indices);
    
    for(int j=0; j<NumEntries; j++)
      Values[j] = (Indices[j] == lrid) ? 1.0 : 0.0;      
  }
}


/*********************************************************************************************************/
void PamgenEnumerateEdges(int numNodesPerElem, int numEdgesPerElem, int numNodesPerEdge,
			  FieldContainer<int> refEdgeToNode,const FieldContainer<double> & nodeCoord,
			  std::vector<long long> globalNodeIds,
			  PamgenMesh & mesh,
			  Epetra_Comm & Comm,
			  /*Output args */
			  std::vector<long long> & globalEdgeIds,
			  std::vector<bool> & edgeIsOwned,
			  std::vector<int> & ownedEdgeIds,
			  FieldContainer<int> &elemToEdge, FieldContainer<int> & elemToEdgeOrient, FieldContainer<int> &edgeToNode,FieldContainer<double> & edgeCoord,
			  long long & numEdgesGlobal) {
  std::vector < topo_entity * > edge_vector;
  std::set < topo_entity * , fecomp > edge_set;
  std::vector < int > edge_comm_procs;
  int rank = Comm.MyPID();

  /***** Hensinger Stuff *****/
  // Calculate edge and ids
  int elct = 0;
  for(long long b = 0; b < mesh.numElemBlk; b++){
    //loop over all elements and push their edges onto a set if they are not there already
    for(long long el = 0; el < mesh.elements[b]; el++){
      std::set< topo_entity *, fecomp > ::iterator fit;
      for (int i=0; i < numEdgesPerElem; i++){
	topo_entity * teof = new topo_entity;
	for(int j = 0; j < numNodesPerEdge;j++)
	  teof->add_node(mesh.elmt_node_linkage[b][el*numNodesPerElem + refEdgeToNode(i,j)],globalNodeIds.data());
	teof->sort();
	fit = edge_set.find(teof);
	if(fit == edge_set.end()){
	  teof->local_id = edge_vector.size();
	  edge_set.insert(teof);
	  //	    printf("[%d] Adding edge %d to element %d/%d\n",Comm.MyPID(),edge_vector.size(),elct,elemToEdge.dimension(0));fflush(stdout);
	  elemToEdge(elct,i)= edge_vector.size();
	  edge_vector.push_back(teof);
	}
	else{
	  elemToEdge(elct,i) = (*fit)->local_id;
	  delete teof;
	}
      }        
      elct ++;
    }
  }

  // Edge to Node connectivity 
  assert(numNodesPerEdge==2);
  edgeToNode.resize(edge_vector.size(), numNodesPerEdge);
  for(unsigned ect = 0; ect != edge_vector.size(); ect++){
    int n[2];
    std::list<long long>::iterator elit=edge_vector[ect]->local_node_ids.begin();
    n[0] = *elit-1;
    elit++;
    n[1] = *elit-1;
    long long nid0 = globalNodeIds[n[0]];
    long long nid1 = globalNodeIds[n[1]];
    long long lo = std::min(nid0,nid1);
    edgeToNode(ect,0)= (lo==nid0)? n[0] : n[1];
    edgeToNode(ect,1)= (lo==nid0)? n[1] : n[0];
  }


  int numEdges = edge_vector.size();
   
  // Calculate global edge and face numbering
  std::string doing_type;
  doing_type = "EDGES";
  calc_global_ids(edge_vector,
		  mesh.comm_node_ids.data(),
		  mesh.node_comm_proc_ids.data(),
		  mesh.node_cmap_node_cnts.data(),
		  mesh.num_node_comm_maps,
		  rank,
		  doing_type);
   
  // Build list of owned global edge ids
  globalEdgeIds.resize(numEdges);
  edgeIsOwned.resize(numEdges);
  int numOwnedEdges=0;
  for (int i=0; i<numEdges; i++) {
    edgeIsOwned[i] = edge_vector[i]->owned;
    globalEdgeIds[i] = edge_vector[i]->global_id;
    if (edgeIsOwned[i]){
      numOwnedEdges++;
    }
  }
  ownedEdgeIds.resize(numOwnedEdges);
  int nedge=0;
  for (int i=0; i<numEdges; i++) {
    if (edgeIsOwned[i]){
      ownedEdgeIds[nedge]=(int)globalEdgeIds[i];
      nedge++;
    }
  }

  // Calculate number of global edges
#ifdef HAVE_MPI
  long long numOwnedEdges_ll = numOwnedEdges;
  Comm.SumAll(&numOwnedEdges_ll,&numEdgesGlobal,1);
#else
  numEdgesGlobal = numEdges;
#endif

  assert(numNodesPerEdge==2);

  /***** Non-Hensinger Stuff ****/
  // Fill out the edge centers (clobbering data if needed)
  edgeCoord.resize(numEdges,mesh.dim);
  for (int i=0; i<numEdges; i++) {
    for(int j=0; j<mesh.dim; j++) {
      edgeCoord(i,j)=0;
      for(int k=0; k<numNodesPerEdge; k++)
	edgeCoord(i,j)+=nodeCoord(edgeToNode(i,k),j)/numNodesPerEdge;
    }   
  }
  
  //#define DEBUG_EDGE_ENUMERATION
#ifdef DEBUG_EDGE_ENUMERATION
  printf("**** New Edge coordinates ***\n");
  for(int i=0; i<edgeCoord.dimension(0); i++)
    printf("[%2d] %10.2f %10.2f\n",i,edgeCoord(i,0),edgeCoord(i,1));
#endif

  // Build the element edge orientation list
  int elid=0;
  for(long long b = 0; b < mesh.numElemBlk; b++){
    for(long long el = 0; el < mesh.elements[b]; el++){
      for (int i=0; i < numEdgesPerElem; i++){	
	int edge = elemToEdge(elid,i);
	int n0 = mesh.elmt_node_linkage[b][el*numNodesPerElem + refEdgeToNode(i,0)]-1;
	elemToEdgeOrient(elid,i) = (n0==edgeToNode(edge,0))?  1 : -1;
      }
      elid++;
    }
  }
 
}



/*********************************************************************************************************/
void EnumerateElements(Epetra_Comm & Comm, int numMyElements, std::vector<long long> & globalElementIds) { 
  long long numMyElements_ll = numMyElements;
  long long myGlobalElementBase = 0;

  Comm.ScanSum(&numMyElements_ll,&myGlobalElementBase,1);
  myGlobalElementBase -=numMyElements;
  //  printf("[%d] MyGlobalElementBase = %lld\n",Comm.MyPID(),myGlobalElementBase);

  globalElementIds.resize(numMyElements);
  for(int i=0; i<numMyElements; i++)
    globalElementIds[i] = myGlobalElementBase + i;
}




/**********************************************************************************/
/**************************** CALCULATE ERROR *************************************/
/**********************************************************************************/
void CalculateError(const Epetra_FEVector & femCoefficients, const Epetra_Map &overlapMap,Epetra_Time & Time,  shards::CellTopology &Pn_cellType, Teuchos::RCP<Basis<double,FieldContainer<double> > > myHGradBasis_rcp, const FieldContainer<int> & elemToNode, const FieldContainer<double> & nodeCoord, int degree) {

  const Epetra_BlockMap & globalMapG = femCoefficients.Map();
  const Epetra_Comm & Comm = globalMapG.Comm();
  int MyPID = Comm.MyPID();
  int numFieldsG = myHGradBasis_rcp->getCardinality();
  int spaceDim = Pn_cellType.getDimension();
  int numElems = elemToNode.dimension(0);
  string msg("Poisson: ");

  if (MyPID == 0) {Time.ResetStartTime();}

  double L2err = 0.0;
  double L2errTot = 0.0;
  double H1err = 0.0;
  double H1errTot = 0.0;
  double Linferr = 0.0;
  double LinferrTot = 0.0;

#ifdef HAVE_MPI
  // Get ghost information from solution to compute error
  Epetra_Import  solnImporter(overlapMap, globalMapG);
  Epetra_Vector  uCoeff(overlapMap);
  uCoeff.Import(femCoefficients, solnImporter, Insert);
#endif

  // Define desired workset size
  int desiredWorksetSize = numElems;
  int numWorksetsErr    = numElems/desiredWorksetSize;

  // When numElems is not divisible by desiredWorksetSize, increase workset count by 1
  if(numWorksetsErr*desiredWorksetSize < numElems) numWorksetsErr += 1;

  // Get cubature points and weights for error calc (may be different from previous)
  Intrepid::DefaultCubatureFactory<double>  cubFactoryErr;
  int cubDegErr = 3*degree;
  Teuchos::RCP<Intrepid::Cubature<double> > cellCubatureErr = cubFactoryErr.create(Pn_cellType, cubDegErr);
  int cubDimErr       = cellCubatureErr->getDimension();
  int numCubPointsErr = cellCubatureErr->getNumPoints();
  Intrepid::FieldContainer<double> cubPointsErr(numCubPointsErr, cubDimErr);
  Intrepid::FieldContainer<double> cubWeightsErr(numCubPointsErr);
  cellCubatureErr->getCubature(cubPointsErr, cubWeightsErr);

  // Evaluate basis values and gradients at cubature points
  Intrepid::FieldContainer<double> uhGVals(numFieldsG, numCubPointsErr);
  Intrepid::FieldContainer<double> uhGrads(numFieldsG, numCubPointsErr, spaceDim);
  myHGradBasis_rcp->getValues(uhGVals, cubPointsErr, Intrepid::OPERATOR_VALUE);
  myHGradBasis_rcp->getValues(uhGrads, cubPointsErr, Intrepid::OPERATOR_GRAD);

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
    Intrepid::FieldContainer<double> cellWorksetEr(worksetSize, numFieldsG, spaceDim);
    Intrepid::FieldContainer<double> worksetApproxSolnCoef(worksetSize, numFieldsG);

    // loop over cells to fill arrays with coordinates and discrete solution coefficient
    int cellCounter = 0;
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      for (int node = 0; node < numFieldsG; node++) {
        cellWorksetEr(cellCounter, node, 0) = nodeCoord( elemToNode(cell, node), 0);
        cellWorksetEr(cellCounter, node, 1) = nodeCoord( elemToNode(cell, node), 1);

        int rowIndex  = elemToNode(cell, node);
#ifdef HAVE_MPI
        worksetApproxSolnCoef(cellCounter, node) = uCoeff.Values()[rowIndex];
#else
        worksetApproxSolnCoef(cellCounter, node) = femCoefficients.Values()[rowIndex];
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
    IntrepidCTools::setJacobian(worksetJacobianE, cubPointsErr, cellWorksetEr,  myHGradBasis_rcp);
    IntrepidCTools::setJacobianInv(worksetJacobInvE, worksetJacobianE );
    IntrepidCTools::setJacobianDet(worksetJacobDetE, worksetJacobianE );

    // map cubature points to physical frame
    Intrepid::FieldContainer<double> worksetCubPoints(worksetSize, numCubPointsErr, cubDimErr);
    IntrepidCTools::mapToPhysicalFrame(worksetCubPoints, cubPointsErr, cellWorksetEr, myHGradBasis_rcp);

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
  Comm.SumAll(&L2err,&L2errTot,1);
  Comm.SumAll(&H1err,&H1errTot,1);
  Comm.MaxAll(&Linferr,&LinferrTot,1);
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

  if(MyPID==0) {std::cout << msg << "Calculate error                             "
                          << Time.ElapsedTime() << " s \n"; Time.ResetStartTime();}


}
