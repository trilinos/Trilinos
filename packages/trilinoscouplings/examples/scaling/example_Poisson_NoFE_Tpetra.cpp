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

/** \file   example_Poisson_NoFE_Tpetra.cpp
    \brief  Example solution of a Poisson equation on a hexahedral mesh using
            nodal (Hgrad) elements.  The system is assembled but not solved.

           This example uses the following Trilinos packages:
    \li     Pamgen to generate a Hexahedral mesh.
    \li     Sacado to form the source term from user-specified manufactured solution.
    \li     Intrepid to build the discretization matrix and right-hand side.
    \li     Tpetra to handle the global matrix and vector.


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
            Converted to Tpetra by I. Kalashnikova (ikalash@sandia.gov)

    \remark Usage:
    \code   ./TrilinosCouplings_examples_scaling_example_Poisson.exe \endcode

    \remark Example driver requires input file named Poisson.xml with Pamgen
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

// Intrepid includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

//Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Operator.hpp"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// Pamgen includes
#include "create_inline_mesh.h"
#include "im_exodusII_l.h"
#include "im_ne_nemesisI_l.h"
#include "pamgen_extras.h"

// Belos includes
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosTpetraAdapter.hpp"

// Sacado includes
#include "Sacado.hpp"

using namespace std;
using namespace Intrepid;
using Tpetra::global_size_t;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TimeMonitor;

/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef Sacado::Fad::SFad<double,3>                             Fad3; //# ind. vars fixed at 3
typedef Intrepid::FunctionSpaceTools                            IntrepidFSTools;
typedef Intrepid::RealSpaceTools<double>                        IntrepidRSTools;
typedef Intrepid::CellTools<double>                             IntrepidCTools;
//Tpetra typedefs
typedef Tpetra::DefaultPlatform::DefaultPlatformType            Platform;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
typedef double                                                  ST;
typedef int                                                     Ordinal;
typedef Tpetra::Map<Ordinal,Ordinal,Node>                       Map;
typedef Tpetra::Export<Ordinal,Ordinal,Node>                    export_type;
typedef Tpetra::Import<Ordinal,Ordinal,Node>                    import_type;
typedef Teuchos::ArrayView<Ordinal>::size_type                  size_type;
typedef Tpetra::CrsMatrix<ST, Ordinal, Ordinal, Node>           sparse_matrix_type;
typedef Tpetra::CrsGraph<Ordinal, Ordinal, Node>                sparse_graph_type;
typedef Tpetra::MultiVector<ST, Ordinal, Ordinal, Node>         MV;
typedef Tpetra::Vector<ST, Ordinal, Ordinal, Node>              vector_type;
typedef Tpetra::Operator<ST, Ordinal, Ordinal, Node>            OP;
typedef Belos::MultiVecTraits<ST, MV>                           MVT;
typedef Belos::OperatorTraits<ST, MV, OP>                       OPT;



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
const Scalar exactSolution(const Scalar& x, const Scalar& y, const Scalar& z);

/** \brief  User-defined material tensor.

    \param  material    [out]   3 x 3 material tensor evaluated at (x,y,z)
    \param  x           [in]    x-coordinate of the evaluation point
    \param  y           [in]    y-coordinate of the evaluation point
    \param  z           [in]    z-coordinate of the evaluation point

    \warning Symmetric and positive definite tensor is required for every (x,y,z).
*/
template<typename Scalar>
void materialTensor(Scalar material[][3], const Scalar&  x, const Scalar&  y, const Scalar&  z);

/** \brief  Computes gradient of the exact solution. Requires user-defined exact solution.

    \param  gradExact  [out]   gradient of the exact solution evaluated at (x,y,z)
    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point
 */
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z);


/** \brief Computes source term: f = -div(A.grad u).  Requires user-defined exact solution
           and material tensor.

    \param  x          [in]    x-coordinate of the evaluation point
    \param  y          [in]    y-coordinate of the evaluation point
    \param  z          [in]    z-coordinate of the evaluation point

    \return Source term corresponding to the user-defined exact solution evaluated at (x,y,z)
 */
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y, Scalar& z);


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

  Teuchos::FancyOStream fos(Teuchos::rcpFromRef(cout));

  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  rank=mpiSession.getRank();
  numProcs=mpiSession.getNProc();


//Get the default communicator and node for Tpetra
//rewrite using with IFDEF for MPI/no MPI??
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Teuchos::Comm<int> > CommT = platform.getComm();
  RCP<Node> node = platform.getNode();
  int MyPID = CommT->getRank();


  //Check number of arguments
  if (argc > 3) {
      cout <<"\n>>> ERROR: Invalid number of arguments.\n\n";
      cout <<"Usage:\n\n";
      cout <<"  ./TrilinosCouplings_examples_scaling_Example_Poisson.exe [meshfile.xml] [solver.xml]\n\n";
      cout <<"   meshfile.xml(optional) - xml file with description of Pamgen mesh\n\n";
      cout <<"   solver.xml(optional) - xml file with ML solver options\n\n";
      exit(1);
   }

 if (MyPID == 0){
  cout \
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
    cout << "PARALLEL executable \n";
  }
#else
  if (MyPID == 0) {
    cout << "SERIAL executable \n";
  }
#endif

/**********************************************************************************/
/********************************** GET XML INPUTS ********************************/
/**********************************************************************************/

  // Command line for xml file, otherwise use default
    std::string   xmlMeshInFileName, xmlSolverInFileName;
    if(argc>=2) xmlMeshInFileName=string(argv[1]);
    else xmlMeshInFileName="Poisson.xml";
    if(argc>=3) xmlSolverInFileName=string(argv[2]);

  // Read xml file into parameter list
    ParameterList inputMeshList;
    ParameterList inputSolverList;

   if(xmlMeshInFileName.length()) {
     if (MyPID == 0) {
      cout << "\nReading parameter list from the XML file \""<<xmlMeshInFileName<<"\" ...\n\n";
     }
     Teuchos::updateParametersFromXmlFile(xmlMeshInFileName, Teuchos::inoutArg(inputMeshList));
     if (MyPID == 0) {
      inputMeshList.print(cout,2,true,true);
      cout << "\n";
     }
    }
    else
    {
      cout << "Cannot read input file: " << xmlMeshInFileName << "\n";
      return 0;
    }

   if(xmlSolverInFileName.length()) {
     if (MyPID == 0)
        cout << "\nReading parameter list from the XML file \""<<xmlSolverInFileName<<"\" ...\n\n";
     Teuchos::updateParametersFromXmlFile(xmlSolverInFileName, Teuchos::inoutArg(inputSolverList));
   } else if (MyPID == 0) cout << "Using default solver values ..." << endl;

   // Get pamgen mesh definition
    std::string meshInput = Teuchos::getParameter<std::string>(inputMeshList,"meshInput");

   // Get Isorropia and Zoltan parameters.
    ParameterList iso_paramlist = inputMeshList.sublist
                                                    ("Isorropia Input") ;
    if (MyPID == 0) {
      cout << "Isorropia/Zoltan parameters" << endl;
        iso_paramlist.print(cout,2,true,true);
    }


/**********************************************************************************/
/***************************** GET CELL TOPOLOGY **********************************/
/**********************************************************************************/

   // Get cell topology for base hexahedron
    shards::CellTopology cellType(shards::getCellTopologyData<shards::Hexahedron<8> >() );

   // Get dimensions
    int numNodesPerElem = cellType.getNodeCount();
    int spaceDim = cellType.getDimension();
    int dim = 3;

/**********************************************************************************/
/******************************* GENERATE MESH ************************************/
/**********************************************************************************/

  if (MyPID == 0) {
    cout << "Generating mesh ... \n\n";
  }

  long long *  node_comm_proc_ids   = NULL;
  long long *  node_cmap_node_cnts  = NULL;
  long long *  node_cmap_ids        = NULL;
  long long ** comm_node_ids        = NULL;
  long long ** comm_node_proc_ids   = NULL;

   // Generate mesh with Pamgen
    long long maxInt = 9223372036854775807LL;
    Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);

    string msg("Poisson: ");


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
       cout << " Number of Global Elements: " << numElemsGlobal << " \n";
       cout << "    Number of Global Nodes: " << numNodesGlobal << " \n\n";
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
    FieldContainer<int> elemToNode(numElems,numNodesPerElem);
    for(long long b = 0; b < numElemBlk; b++){
      for(long long el = 0; el < elements[b]; el++){
        for (int j=0; j<numNodesPerElem; j++) {
          elemToNode(telct,j) = elmt_node_linkage[b][el*numNodesPerElem + j]-1;
        }
        telct ++;
      }
    }


   // Read node coordinates and place in field container
    FieldContainer<double> nodeCoord(numNodes,dim);
    double * nodeCoordx = new double [numNodes];
    double * nodeCoordy = new double [numNodes];
    double * nodeCoordz = new double [numNodes];
    im_ex_get_coord_l(id,nodeCoordx,nodeCoordy,nodeCoordz);
    for (int i=0; i<numNodes; i++) {
      nodeCoord(i,0)=nodeCoordx[i];
      nodeCoord(i,1)=nodeCoordy[i];
      nodeCoord(i,2)=nodeCoordz[i];
    }
    delete [] nodeCoordx;
    delete [] nodeCoordy;
    delete [] nodeCoordz;


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



    //Calculate global node ids
    long long * globalNodeIds = new long long[numNodes];
    bool * nodeIsOwned = new bool[numNodes];

    calc_global_node_ids(globalNodeIds,
                         nodeIsOwned,
                         numNodes,
                         num_node_comm_maps,
                         node_cmap_node_cnts,
                         node_comm_proc_ids,
                         comm_node_ids,
                         rank);



   // Container indicating whether a node is on the boundary (1-yes 0-no)
    FieldContainer<int> nodeOnBoundary(numNodes);

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

   // Get numerical integration points and weights
    DefaultCubatureFactory<double>  cubFactory;
    int cubDegree = 2;
    RCP<Cubature<double> > hexCub = cubFactory.create(cellType, cubDegree);

    int cubDim       = hexCub->getDimension();
    int numCubPoints = hexCub->getNumPoints();

    FieldContainer<double> cubPoints(numCubPoints, cubDim);
    FieldContainer<double> cubWeights(numCubPoints);

    hexCub->getCubature(cubPoints, cubWeights);



/**********************************************************************************/
/*********************************** GET BASIS ************************************/
/**********************************************************************************/

   // Define basis
     Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> > hexHGradBasis;
     int numFieldsG = hexHGradBasis.getCardinality();
     FieldContainer<double> HGBValues(numFieldsG, numCubPoints);
     FieldContainer<double> HGBGrads(numFieldsG, numCubPoints, spaceDim);

  // Evaluate basis values and gradients at cubature points
     hexHGradBasis.getValues(HGBValues, cubPoints, OPERATOR_VALUE);
     hexHGradBasis.getValues(HGBGrads, cubPoints, OPERATOR_GRAD);




/**********************************************************************************/
/********************* BUILD MAPS FOR GLOBAL SOLUTION *****************************/
/**********************************************************************************/
    //timer for building global maps
    RCP<Teuchos::Time> timerBuildGlobalMaps = TimeMonitor::getNewTimer("Build global maps: Total Time");
    Teuchos::Array<int> ownedGIDs;
    RCP<const Map > globalMapGT;
    {
    TimeMonitor timerBuildGlobalMapsL(*timerBuildGlobalMaps);
    // Count owned nodes
    int ownedNodes=0;
    for(int i=0;i<numNodes;i++){
      if(nodeIsOwned[i]) ownedNodes++;
    }


    // Build a list of the OWNED global ids...
    // NTS: will need to switch back to long long
    ownedGIDs.resize(ownedNodes);
    int oidx=0;
    for(int i=0;i<numNodes;i++)
      if(nodeIsOwned[i]){
        ownedGIDs[oidx]=(int)globalNodeIds[i];
        oidx++;
      }
    RCP<Teuchos::Time> timerBuildGlobalMaps1 = TimeMonitor::getNewTimer("Build global maps: globalMapGT");
    {
    TimeMonitor timerBuildGlobalMaps1L(*timerBuildGlobalMaps1);
    //Generate Tpetra map for nodes
    globalMapGT = rcp(new Map(-1, ownedGIDs(), 0, CommT));
    }
    }


/**********************************************************************************/
/********************* BUILD MAPS FOR OVERLAPPED SOLUTION *************************/
/**********************************************************************************/
    RCP<Teuchos::Time> timerBuildOverlapMaps = TimeMonitor::getNewTimer("Build overlapped maps: Total Time");
    Teuchos::Array<int> overlappedGIDs;
    RCP<const Map > overlappedMapGT;
    RCP<const export_type> exporterT;
    {
    TimeMonitor timerBuildOverlapMapsL(*timerBuildOverlapMaps);
    // Count owned nodes
    int overlappedNodes=numNodes;

    // Build a list of the OVERLAPPED global ids...
    overlappedGIDs.resize(overlappedNodes);
    for(int i=0;i<numNodes;i++)
        overlappedGIDs[i]=(int)globalNodeIds[i];

    RCP<Teuchos::Time> timerBuildOverlapMaps1 = TimeMonitor::getNewTimer("Build overlapped maps: overlappedMapGT");
    {
    TimeMonitor timerBuildOverlapMaps1L(*timerBuildOverlapMaps1);
    //Generate Tpetra map for nodes
    overlappedMapGT = rcp(new Map(-1, overlappedGIDs(), 0, CommT));
    }
    //build Tpetra Export/Import
    RCP<Teuchos::Time> timerBuildOverlapMaps2 = TimeMonitor::getNewTimer("Build overlapped maps: exporterT");
    {
    TimeMonitor timerBuildOverlapMaps2L(*timerBuildOverlapMaps2);
    exporterT = rcp(new export_type(overlappedMapGT, globalMapGT));
    }
    RCP<Teuchos::Time> timerBuildOverlapMaps3 = TimeMonitor::getNewTimer("Build overlapped maps: importerT");
    {
    TimeMonitor timerBuildOverlapMaps3L(*timerBuildOverlapMaps3);
    RCP<const import_type> importerT = rcp(new import_type(overlappedMapGT, globalMapGT));
    }
    }

/**********************************************************************************/
/********************* BUILD GRAPH FOR OVERLAPPED SOLUTION *************************/
/**********************************************************************************/

    RCP<Teuchos::Time> timerBuildOverlapGraph = TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: Total Time");
    RCP<sparse_matrix_type> gl_StiffMatrixT;
    RCP<vector_type> gl_rhsVectorT;
    RCP<sparse_matrix_type> StiffMatrixT;
    RCP<vector_type> rhsVectorT;
    {
    TimeMonitor timerBuildOverlapGraphL(*timerBuildOverlapGraph);
    RCP<sparse_graph_type> overlappedGraphT;
    RCP<sparse_graph_type> ownedGraphT;
    RCP<Teuchos::Time> timerBuildOverlapGraph2 = TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: create overlappedGraphT & ownedGraphT");
    {
    TimeMonitor timerBuildOverlapGraph2L(*timerBuildOverlapGraph2);
    //construct Tpetra CrsGraphs
    overlappedGraphT= rcp(new sparse_graph_type(overlappedMapGT, 0));
    ownedGraphT= rcp(new sparse_graph_type(globalMapGT, 0));
    }

    {

    RCP<Teuchos::Time> timerBuildOverlapGraph3 = TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: insertGlobalIndices overlappedGraphT");
    {
    TimeMonitor timerBuildOverlapGraph3L(*timerBuildOverlapGraph3);
      // Define desired workset size and count how many worksets there are on this processor's mesh block
      int desiredWorksetSize = numElems;                      // change to desired workset size!
      //int desiredWorksetSize = 100;                      // change to desired workset size!
      int numWorksets        = numElems/desiredWorksetSize;
      for(int workset = 0; workset < numWorksets; workset++){

        // Compute cell numbers where the workset starts and ends
        int worksetSize  = 0;
        int worksetBegin = (workset + 0)*desiredWorksetSize;
        int worksetEnd   = (workset + 1)*desiredWorksetSize;

        // When numElems is not divisible by desiredWorksetSize, the last workset ends at numElems
        worksetEnd   = (worksetEnd <= numElems) ? worksetEnd : numElems;

        // Now we know the actual workset size and can allocate the array for the cell nodes
        worksetSize  = worksetEnd - worksetBegin;

        //"WORKSET CELL" loop: local cell ordinal is relative to numElems
        for(int cell = worksetBegin; cell < worksetEnd; cell++){

          // Compute cell ordinal relative to the current workset
          int worksetCellOrdinal = cell - worksetBegin;

          // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
          for (int cellRow = 0; cellRow < numFieldsG; cellRow++){

            int localRow  = elemToNode(cell, cellRow);
            //globalRow for Tpetra Graph
            global_size_t globalRowT = globalNodeIds[localRow];

            // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
            for (int cellCol = 0; cellCol < numFieldsG; cellCol++){

              int localCol  = elemToNode(cell, cellCol);
              int globalCol = globalNodeIds[localCol];
              //create ArrayView globalCol object for Tpetra
              Teuchos::ArrayView<int> globalColAV = Teuchos::arrayView(&globalCol, 1);

              //Update Tpetra overlap Graph
              overlappedGraphT->insertGlobalIndices(globalRowT, globalColAV);

            }// *** cell col loop ***
          }// *** cell row loop ***
        }// *** workset cell loop **
      }// *** workset loop ***
      }
    }

    RCP<Teuchos::Time> timerBuildOverlapGraph3 = TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: fillcomplete overlappedGraphT");
    {
    TimeMonitor timerBuildOverlapGraph3L(*timerBuildOverlapGraph3);
    //Fillcomplete Tpetra overlap Graph
    overlappedGraphT->fillComplete();
    }

    RCP<Teuchos::Time> timerBuildOverlapGraph4 = TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: export and fillcomplete ownedGraphT");
    {
    TimeMonitor timerBuildOverlapGraph4L(*timerBuildOverlapGraph4);
    //build global map - Tpetra
    ownedGraphT->doExport(*overlappedGraphT, *exporterT, Tpetra::INSERT);
    ownedGraphT->fillComplete();
    }

    RCP<Teuchos::Time> timerBuildOverlapGraph5 = TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: create and zero global/local matrices/vectors");
    {
    TimeMonitor timerBuildOverlapGraph5L(*timerBuildOverlapGraph5);
    //build Tpetra analogs of Stiffness matrix and rhsVector
    gl_StiffMatrixT = rcp(new sparse_matrix_type(ownedGraphT.getConst()));
    gl_rhsVectorT = rcp(new vector_type(globalMapGT));

    StiffMatrixT = rcp(new sparse_matrix_type(overlappedGraphT.getConst()));
    rhsVectorT = rcp(new vector_type(overlappedMapGT));
    StiffMatrixT->setAllToScalar(0.0);
    rhsVectorT->putScalar(0.0);
    }
}

#ifdef DUMP_DATA
/**********************************************************************************/
/**** PUT COORDINATES AND NODAL VALUES IN ARRAYS FOR OUTPUT (FOR PLOTTING ONLY) ***/
/**********************************************************************************/


   //Put coordinate in Tpetra multivector for ourput
   RCP<MV> nCoordT = rcp(new MV(globalMapGT, 3));
   RCP<MV> nBoundT = rcp(new MV(globalMapGT, 1));

     int indOwned = 0;
     for (int inode=0; inode<numNodes; inode++) {
       if (nodeIsOwned[inode]) {
          nCoordT->replaceLocalValue(indOwned, 0, nodeCoord(inode,0));
          nCoordT->replaceLocalValue(indOwned, 1, nodeCoord(inode,1));
          nCoordT->replaceLocalValue(indOwned, 2, nodeCoord(inode,2));
          nBoundT->replaceLocalValue(indOwned, 0, nodeOnBoundary(inode));
          indOwned++;
       }
     }
     //Write Tpetra multivectors to MatrixMarket files
      Tpetra::MatrixMarket::Writer<sparse_matrix_type >::writeDenseFile("coordsT.dat",nCoordT);
      Tpetra::MatrixMarket::Writer<sparse_matrix_type >::writeDenseFile("nodesOnBoundT.dat",nBoundT);


     //Put element to node mapping in Tpetra multivector for output
     RCP<const Map> globalMapElemT = rcp(new Map(numElemsGlobal, numElems, 0, CommT));
     RCP<MV>elem2nodeT = rcp(new MV(globalMapElemT, numNodesPerElem));
     for (int ielem=0; ielem<numElems; ielem++) {
        for (int inode=0; inode<numNodesPerElem; inode++) {
          elem2nodeT->replaceLocalValue(ielem, inode, globalNodeIds[elemToNode(ielem,inode)]);
        }
      }

     //Write Tpetra multivector to MatrixMarket file
      Tpetra::MatrixMarket::Writer<sparse_matrix_type >::writeDenseFile("elem2nodeT.dat",elem2nodeT);

#endif

/**********************************************************************************/
/************************** DIRICHLET BC SETUP ************************************/
/**********************************************************************************/

  //timer for Dirichlet BC setup
  RCP<Teuchos::Time> timerDirichletBC = TimeMonitor::getNewTimer("Get Dirichlet boundary values: Total Time");
  int numBCNodes = 0;
  RCP<MV> vT;
  Teuchos::Array<int> BCNodes;
  {
  TimeMonitor timerDirichletBCL(*timerDirichletBC);
  for (int inode = 0; inode < numNodes; inode++){
     if (nodeOnBoundary(inode) && nodeIsOwned[inode]){
        numBCNodes++;
     }
  }


   //Vector for use in applying BCs - Tpetra
   RCP<Teuchos::Time> timerDirichletBC1 = TimeMonitor::getNewTimer("Get Dirichlet boundary values: create and zero vT");
   {
   TimeMonitor timerDirichletBC1L(*timerDirichletBC1);
   vT = rcp(new MV(globalMapGT, true));
   vT->putScalar(0.0);
   }
   RCP<Teuchos::Time> timerDirichletBC2 = TimeMonitor::getNewTimer("Get Dirichlet boundary values: set vT for Dirichlet nodes");
   {
   TimeMonitor timerDirichletBC2L(*timerDirichletBC2);
   // Set v to boundary values on Dirichlet nodes
    //BCNodes = new int [numBCNodes];
    BCNodes.resize(numBCNodes);
    int indbc=0;
    int iOwned=0;
    for (int inode=0; inode<numNodes; inode++){
      if (nodeIsOwned[inode]){
        if (nodeOnBoundary(inode)){
           BCNodes[indbc]=iOwned;
           indbc++;
           double x  = nodeCoord(inode, 0);
           double y  = nodeCoord(inode, 1);
           double z  = nodeCoord(inode, 2);
           vT->replaceLocalValue(iOwned, 0, exactSolution(x,y,z));
        }
         iOwned++;
      }
    }
   }
   }



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
    cout << "Building discretization matrix and right hand side... \n\n";
    cout << "\tDesired workset size:                 " << desiredWorksetSize <<"\n";
    cout << "\tNumber of worksets (per processor):   " << numWorksets <<"\n\n";
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
     FieldContainer<double> cellWorkset(worksetSize, numNodesPerElem, spaceDim);

    // Copy coordinates into cell workset
    int cellCounter = 0;
    for(int cell = worksetBegin; cell < worksetEnd; cell++){
      for (int node = 0; node < numNodesPerElem; node++) {
        cellWorkset(cellCounter, node, 0) = nodeCoord( elemToNode(cell, node), 0);
        cellWorkset(cellCounter, node, 1) = nodeCoord( elemToNode(cell, node), 1);
        cellWorkset(cellCounter, node, 2) = nodeCoord( elemToNode(cell, node), 2);
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
    IntrepidFSTools::integrate<double>(worksetStiffMatrix,                      // (DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
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
 /*                         Assemble into Global Matrix                            */
 /**********************************************************************************/

  RCP<Teuchos::Time> timerAssembleGlobalMatrix = TimeMonitor::getNewTimer("Assemble global matrix: Total Time");
    {
    TimeMonitor timerAssembleGlobalMatrixL(*timerAssembleGlobalMatrix);
    //"WORKSET CELL" loop: local cell ordinal is relative to numElems
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // Compute cell ordinal relative to the current workset
      int worksetCellOrdinal = cell - worksetBegin;

      // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
      for (int cellRow = 0; cellRow < numFieldsG; cellRow++){

        int localRow  = elemToNode(cell, cellRow);
        int globalRow = globalNodeIds[localRow];
        double sourceTermContribution =  worksetRHS(worksetCellOrdinal, cellRow);
        Teuchos::ArrayView<double> sourceTermContributionAV = Teuchos::arrayView(&sourceTermContribution, 1);

        rhsVectorT->sumIntoGlobalValue(globalRow, sourceTermContribution);

        // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
        for (int cellCol = 0; cellCol < numFieldsG; cellCol++){

          int localCol  = elemToNode(cell, cellCol);
          int globalCol = globalNodeIds[localCol];
          Teuchos::ArrayView<int> globalColAV = Teuchos::arrayView(&globalCol, 1);
          double operatorMatrixContribution = worksetStiffMatrix(worksetCellOrdinal, cellRow, cellCol);
          Teuchos::ArrayView<double> operatorMatrixContributionAV = Teuchos::arrayView(&operatorMatrixContribution, 1);

          StiffMatrixT->sumIntoGlobalValues(globalRow, globalColAV, operatorMatrixContributionAV);

        }// *** cell col loop ***
      }// *** cell row loop ***
    }// *** workset cell loop **
   } //***timer stop ***
  }// *** workset loop ***


/**********************************************************************************/
/********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
/**********************************************************************************/


   RCP<Teuchos::Time> timerAssembMultProc = TimeMonitor::getNewTimer("Global assembly over multiple processors: Total Time");
   {
   TimeMonitor timerAssembMultProcL(*timerAssembMultProc);
   gl_StiffMatrixT->setAllToScalar(0.0);
   gl_StiffMatrixT->doExport(*StiffMatrixT, *exporterT, Tpetra::ADD);
   //if target of export has static graph, no need to do setAllToScalar(0.0); replace mode export will clobber values
   gl_StiffMatrixT->fillComplete();

   gl_rhsVectorT->putScalar(0.0);
   gl_rhsVectorT->doExport(*rhsVectorT, *exporterT, Tpetra::ADD);
   }


/**********************************************************************************/
/******************************* ADJUST MATRIX DUE TO BC **************************/
/**********************************************************************************/

   RCP<Teuchos::Time> timerAdjustMatrixBC = TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: Total Time");
   {
   TimeMonitor timerAdjustMatrixBCL(*timerAdjustMatrixBC);

   RCP<MV> rhsDirT;
   RCP<Teuchos::Time> timerAdjustMatrixBC1 = TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: apply stiffmatrix to vT");
   {
   TimeMonitor timerAdjustMatrixBC1L(*timerAdjustMatrixBC1);
   //Apply stiffness matrix to v - Tpetra
   rhsDirT = rcp(new MV(globalMapGT, true));
   gl_StiffMatrixT->apply(*vT.getConst(), *rhsDirT);
   gl_StiffMatrixT->resumeFill(); //call to start messing w/ fillCompleted Matrix again
   }

   RCP<Teuchos::Time> timerAdjustMatrixBC2 = TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: update rhs");
   {
   TimeMonitor timerAdjustMatrixBC2L(*timerAdjustMatrixBC2);
   //Update right-hand side - Tpetra
   gl_rhsVectorT->update(-1.0, *rhsDirT, 1.0);
   }

    Teuchos::ArrayRCP<const double> vTArrayRCP = vT->getData(0);
   RCP<Teuchos::Time> timerAdjustMatrixBC3 = TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: adjust rhs to Dirichlet BCs");
   {
   TimeMonitor timerAdjustMatrixBC3L(*timerAdjustMatrixBC3);
    // Adjust rhs due to Dirichlet boundary conditions
   for (int inode=0; inode<numNodes; inode++){
      if (nodeIsOwned[inode]){
        if (nodeOnBoundary(inode)){
           // get global node number
           int gni = globalNodeIds[inode];
           int lidT = globalMapGT->getLocalElement(gni);
           double v_valT = vTArrayRCP[lidT];
           gl_rhsVectorT->replaceGlobalValue(gni, v_valT);
        }
      }
    }
    }

   RCP<Teuchos::Time> timerAdjustMatrixBC4 = TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: zero out rows/cols of stiffness for Dirichlet edges");
   {
   TimeMonitor timerAdjustMatrixBC4L(*timerAdjustMatrixBC4);
    // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
   //  and add one to diagonal.
   // The following is the Tpetra analog of Apply_OAZToMatrix function

   /* Find the local column numbers to nuke */
   RCP<const Map> ColMapT = gl_StiffMatrixT->getColMap();
   RCP<const Map> globalMapT = rcp(new Map(gl_StiffMatrixT->getGlobalNumCols(), 0, CommT));
   // create the exporter from this proc's column map to global 1-1 column map
   RCP<const export_type > ExporterT = rcp(new export_type(ColMapT, globalMapT));
   // create a vector of global column indices that we will export to
   RCP<Tpetra::Vector<int, Ordinal, Ordinal, Node> > globColsToZeroT = rcp(new Tpetra::Vector<int, Ordinal, Ordinal, Node>(globalMapT));
   // create a vector of local column indices that we will export from
   RCP<Tpetra::Vector<int, Ordinal, Ordinal, Node> > myColsToZeroT = rcp(new Tpetra::Vector<int, Ordinal, Ordinal, Node>(ColMapT));
   myColsToZeroT->putScalar(0);
   // flag all local columns corresponding to the local rows specified
   for (int i=0; i < numBCNodes; i++){
       int globalRow = gl_StiffMatrixT->getRowMap()->getGlobalElement(BCNodes[i]);
       int localCol = gl_StiffMatrixT->getColMap()->getLocalElement(globalRow);
       myColsToZeroT->replaceLocalValue(localCol, 1);
   }
   // export to the global column map
    globColsToZeroT->doExport(*myColsToZeroT,*ExporterT,Tpetra::ADD);
    // now import from the global column map to the local column map
    myColsToZeroT->doImport(*globColsToZeroT,*ExporterT,Tpetra::INSERT);

    Teuchos::Array<double> values;
    Teuchos::Array<int> indices;
    Teuchos::ArrayRCP<const int> myColsToZeroArrayRCP = myColsToZeroT->getData(0);
    size_t NumEntries = 0;
    /* Zero the columns */
    for (int i=0; i < gl_StiffMatrixT->getNodeNumRows(); i++) {
       NumEntries = gl_StiffMatrixT->getNumEntriesInLocalRow(i);
       values.resize(NumEntries);
       indices.resize(NumEntries);
       gl_StiffMatrixT->getLocalRowCopy(i, indices(), values(), NumEntries);
       //Matrix.ExtractMyRowView(i,numEntries,vals,cols);
       for (int j=0; j < NumEntries; j++){
          //Teuchos::ArrayRCP<const int> myColsToZeroj = myColsToZeroT->getData();
          if (myColsToZeroArrayRCP[indices[j]] == 1)
              values[j] = 0.0;
       }
       gl_StiffMatrixT->replaceLocalValues(i, indices(), values());
    }/*end for*/

    /* Zero the rows, add ones to diagonal */
   for (int i = 0; i<numBCNodes; i++) {
      NumEntries = gl_StiffMatrixT->getNumEntriesInLocalRow(BCNodes[i]);
      indices.resize(NumEntries);
      values.resize(NumEntries);
      gl_StiffMatrixT->getLocalRowCopy(BCNodes[i], indices(), values(), NumEntries);
      int globalRow = gl_StiffMatrixT->getRowMap()->getGlobalElement(BCNodes[i]);
      int localCol = gl_StiffMatrixT->getColMap()->getLocalElement(globalRow);
      for (int j = 0; j<NumEntries; j++){
         values[j] = 0.0;
         if (indices[j] == localCol)
            values[j] = 1.0;
      }
      gl_StiffMatrixT->replaceLocalValues(BCNodes[i], indices(), values());
   }
   }
   }
   gl_StiffMatrixT->fillComplete();

   //save global stiffness matrix and rhs vector to matrix market file
   Tpetra::MatrixMarket::Writer<sparse_matrix_type >::writeSparseFile("gl_StiffMatrixT.dat",gl_StiffMatrixT);
   Tpetra::MatrixMarket::Writer<sparse_matrix_type >::writeDenseFile("gl_rhsVectorT.dat",gl_rhsVectorT);

   /**********************************************************************************/
  /*******************SOLVE GLOBAL SYSTEM USING BELOS + CG **************************/
  /**********************************************************************************/

   if (MyPID == 0) {
      cout << "Solving global linear system using Belos with unpreconditioned CG... \n\n";
    }
    RCP<Teuchos::Time> timerBelosSolve = TimeMonitor::getNewTimer("Solve global system using Belos + CG: Total Time");
    RCP<vector_type> gl_solVectorT;
    RCP<Belos::SolverManager<ST, MV, OP > >solver;
    double tol = 1e-10;
    Belos::ReturnType ret;
    {
    TimeMonitor timerBelosSolveL(*timerBelosSolve);
    //Create vector to store solution & set initial guess
    gl_solVectorT = rcp(new vector_type(globalMapGT));
    gl_solVectorT->putScalar(1.0);

    //create parameter list for block CG solver manager
    ParameterList belosList;
    belosList.set("Block Size", 1);
    belosList.set("Maximum Iterations", 200);
    belosList.set("Convergence Tolerance", tol);

   //construct unpreconditioned linear problem
    RCP<Belos::LinearProblem<ST, MV, OP > > problem = rcp(new Belos::LinearProblem<ST, MV, OP >(gl_StiffMatrixT, gl_solVectorT, gl_rhsVectorT));
   //set problem
   bool set = problem->setProblem ();
   if (set == false) {
      cout << endl << "ERROR: Belos::LinearProblem failed to set up correctly!" << endl;
      return -1;
   }
   //create an iterative solver manager
   solver = rcp(new Belos::PseudoBlockCGSolMgr<ST, MV, OP >(problem, rcp(&belosList, false)));

   //Perform solve
   ret = solver->solve();
   }

   //Get # iterations for this solve
   int numIters = solver->getNumIters();
   if (MyPID == 0) cout << "Number of iterations performed for the linear solve: "<< numIters << endl;

   //compute actual residuals
   bool badRes = false;
   vector<double> actual_resids(1);
   vector<double> rhs_norm(1);
   RCP<MV> residT = rcp(new MV(globalMapGT, 1));
   OPT::Apply(*gl_StiffMatrixT, *gl_solVectorT, *residT);
   MVT::MvAddMv(-1.0, *residT, 1.0, *gl_rhsVectorT, *residT);
   MVT::MvNorm(*residT, actual_resids);
   MVT::MvNorm(*gl_rhsVectorT, rhs_norm);
   double actRes = actual_resids[0]/rhs_norm[0];
   if (MyPID == 0) cout << "Actual residual (normalized): " << actRes << endl;
   if (actRes > tol) badRes = true;

   if (ret!=Belos::Converged || badRes) {
      if (MyPID == 0) cout << "ERROR: Belos failed to converge!" << endl;
      return -1;
   }
   if (MyPID == 0) cout << "Belos converged!" << endl;


   //write gl_solVector to MatrixMarket file
   Tpetra::MatrixMarket::Writer<sparse_matrix_type >::writeDenseFile("gl_solVectorT.dat", gl_solVectorT);


   //summarize timings
   TimeMonitor::summarize( cout );

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
   delete [] globalNodeIds;
   delete [] nodeIsOwned;
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
const Scalar exactSolution(const Scalar& x, const Scalar& y, const Scalar& z) {

  // Patch test: tri-linear function is in the FE space and should be recovered
   return 1. + x + y + z + x*y + x*z + y*z + x*y*z;

  // Analytic solution with homogeneous Dirichlet boundary data
  // return sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z)*exp(x+y+z);

  // Analytic solution with inhomogeneous Dirichlet boundary data
  // return exp(x + y + z)/(1. + x*y + y*z + x*y*z);
}


template<typename Scalar>
void materialTensor(Scalar material[][3], const Scalar& x, const Scalar& y, const Scalar& z) {

  material[0][0] = 1.;
  material[0][1] = 0.;
  material[0][2] = 0.;
  //
  material[1][0] = 0.;
  material[1][1] = 1.;
  material[1][2] = 0.;
  //
  material[2][0] = 0.;
  material[2][1] = 0.;
  material[2][2] = 1.;
}

/**********************************************************************************/
/************** AUXILIARY FUNCTIONS FROM EXACT SOLUTION ***************************/
/**********************************************************************************/

/************ Grad of Exact Solution ****************/
template<typename Scalar>
void exactSolutionGrad(Scalar gradExact[3], const Scalar& x, const Scalar& y, const Scalar& z) {

  // To enable derivatives of the gradient (i.e., 2nd derivatives of the exact solution) need 2 levels of fad types
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

/************ Source Term (RHS) ****************/
template<typename Scalar>
const Scalar sourceTerm(Scalar& x, Scalar& y, Scalar& z){

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
  u = exactSolution(x, y, z);
  exactSolutionGrad(grad_u, x, y, z);

  // Get material tensor
  materialTensor<Scalar>(material, x, y, z);

  // Compute total flux = (A.grad u)
  for(int i = 0; i < 3; i++){

    // Add diffusive flux
    for(int j = 0; j < 3; j++){
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

/************ Material Tensor ****************/
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor(ArrayOut &        matTensorValues,
                             const ArrayIn &   evaluationPoints){

  int numWorksetCells  = evaluationPoints.dimension(0);
  int numPoints        = evaluationPoints.dimension(1);
  int spaceDim         = evaluationPoints.dimension(2);

  double material[3][3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      materialTensor<double>(material, x, y, z);

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

      Sacado::Fad::SFad<double,3> x = evaluationPoints(cell, pt, 0);
      Sacado::Fad::SFad<double,3> y = evaluationPoints(cell, pt, 1);
      Sacado::Fad::SFad<double,3> z = evaluationPoints(cell, pt, 2);

      sourceTermValues(cell, pt) = sourceTerm<Sacado::Fad::SFad<double,3> >(x, y, z).val();
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
      double z = evaluationPoints(cell, pt, 2);

      exactSolutionValues(cell, pt) = exactSolution<double>(x, y, z);
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

  double gradient[3];

  for(int cell = 0; cell < numWorksetCells; cell++){
    for(int pt = 0; pt < numPoints; pt++){

      double x = evaluationPoints(cell, pt, 0);
      double y = evaluationPoints(cell, pt, 1);
      double z = evaluationPoints(cell, pt, 2);

      exactSolutionGrad<double>(gradient, x, y, z);

      for(int row = 0; row < spaceDim; row++){
        exactSolutionGradValues(cell, pt, row) = gradient[row];
      }
    }
  }
}



