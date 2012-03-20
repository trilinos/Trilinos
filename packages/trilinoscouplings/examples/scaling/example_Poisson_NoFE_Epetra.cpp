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

/** \file   example_Poisson_NoFE_Epetra.cpp
    \brief  Example solution of a Poisson equation on a hexahedral mesh using
            nodal (Hgrad) elements.  The system is assembled but not solved.

           This example uses the following Trilinos packages:
    \li     Pamgen to generate a Hexahedral mesh.
    \li     Sacado to form the source term from user-specified manufactured solution.
    \li     Intrepid to build the discretization matrix and right-hand side.
    \li     Epetra to handle the global matrix and vector.


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
#include "Epetra_Export.h"

// Teuchos includes
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_VectorOut.h"

// Pamgen includes
#include "create_inline_mesh.h"
#include "im_exodusII_l.h"
#include "im_ne_nemesisI_l.h"
#include "pamgen_extras.h"

// AztecOO includes
#include "AztecOO.h"

// ML Includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

// Sacado includes
#include "Sacado.hpp"

using namespace std;
using namespace Intrepid;


/*********************************************************/
/*                     Typedefs                          */
/*********************************************************/
typedef Sacado::Fad::SFad<double,3>      Fad3; //# ind. vars fixed at 3
typedef Intrepid::FunctionSpaceTools     IntrepidFSTools;
typedef Intrepid::RealSpaceTools<double> IntrepidRSTools;
typedef Intrepid::CellTools<double>      IntrepidCTools;



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
                                 Epetra_CrsMatrix   & A,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
                                 double & TotalErrorExactSol);



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

#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv,0);
  rank=mpiSession.getRank();
  numProcs=mpiSession.getNProc();
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
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
    else xmlMeshInFileName="Poisson.xml";
    if(argc>=3) xmlSolverInFileName=string(argv[2]);

  // Read xml file into parameter list
    Teuchos::ParameterList inputMeshList;
    Teuchos::ParameterList inputSolverList;

   if(xmlMeshInFileName.length()) {
     if (MyPID == 0) {
      std::cout << "\nReading parameter list from the XML file \""<<xmlMeshInFileName<<"\" ...\n\n";
     }
      Teuchos::updateParametersFromXmlFile(xmlMeshInFileName,&inputMeshList);
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
     Teuchos::updateParametersFromXmlFile(xmlSolverInFileName,&inputSolverList);
   } else if (MyPID == 0) std::cout << "Using default solver values ..." << std::endl;

   // Get pamgen mesh definition
    std::string meshInput = Teuchos::getParameter<std::string>(inputMeshList,"meshInput");

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
    std::cout << "Generating mesh ... \n\n";
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
       std::cout << " Number of Global Elements: " << numElemsGlobal << " \n";
       std::cout << "    Number of Global Nodes: " << numNodesGlobal << " \n\n";
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


    //cout << "elemToNode, proc #" << Comm.MyPID() <<  endl; 
   // cout << elemToNode << endl;  
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

   // cout << "nodeCoord, proc #" << Comm.MyPID() <<  endl; 
   // cout << nodeCoord << endl;  

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

/*    cout << "nodeIsOwned, proc #" << Comm.MyPID() << endl; 
    for (int i = 0; i < numNodes; i++) 
      cout <<"node " << i << ": " <<  nodeIsOwned[i] << endl; 

    cout << "globalNodeIds, proc #" << Comm.MyPID()<< endl; 
    for (int i = 0; i < numNodes; i++) 
      cout <<"node " << i << ": " <<  globalNodeIds[i] << endl; 
*/

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
    Teuchos::RCP<Cubature<double> > hexCub = cubFactory.create(cellType, cubDegree); 

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

    Teuchos::RCP<Teuchos::Time> timerBuildGlobalMaps = Teuchos::TimeMonitor::getNewTimer("Build global maps: Total Time");
    int *ownedGIDs; 
    Teuchos::RCP<Epetra_Map> globalMapG; 
    {
    Teuchos::TimeMonitor timerBuildGlobalMapsL(*timerBuildGlobalMaps);
    // Count owned nodes
    int ownedNodes=0;
    for(int i=0;i<numNodes;i++)
      if(nodeIsOwned[i]) ownedNodes++;

    // Build a list of the OWNED global ids...
    // NTS: will need to switch back to long long
    ownedGIDs=new int[ownedNodes];    
    int oidx=0;
    for(int i=0;i<numNodes;i++)
      if(nodeIsOwned[i]){
        ownedGIDs[oidx]=(int)globalNodeIds[i];
        oidx++;
      }

    Teuchos::RCP<Teuchos::Time> timerBuildGlobalMaps1 = Teuchos::TimeMonitor::getNewTimer("Build global maps: globalMapG"); 
    {
    Teuchos::TimeMonitor timerBuildGlobalMaps1L(*timerBuildGlobalMaps1);  
    // Generate epetra map for nodes
     globalMapG = Teuchos::rcp(new Epetra_Map(-1,ownedNodes,ownedGIDs,0,Comm));
    }
     }
    

/**********************************************************************************/
/********************* BUILD MAPS FOR OVERLAPPED SOLUTION *************************/
/**********************************************************************************/

    Teuchos::RCP<Teuchos::Time> timerBuildOverlapMaps = Teuchos::TimeMonitor::getNewTimer("Build overlapped maps: Total Time");
    int *overlappedGIDs; 
    Teuchos::RCP<Epetra_Map> overlappedMapG;
    Teuchos::RCP<Epetra_Export> exporter;  
    {
    Teuchos::TimeMonitor timerBuildOverlapMapsL(*timerBuildOverlapMaps); 
    // Count owned nodes
    int overlappedNodes=numNodes;

    // Build a list of the OVERLAPPED global ids...
    overlappedGIDs=new int[overlappedNodes];    
    for(int i=0;i<numNodes;i++)
        overlappedGIDs[i]=(int)globalNodeIds[i];
    Teuchos::RCP<Teuchos::Time> timerBuildOverlapMaps1 = Teuchos::TimeMonitor::getNewTimer("Build overlapped maps: overlappedMapG");
    {
    Teuchos::TimeMonitor timerBuildOverlapMaps1L(*timerBuildOverlapMaps1);
    // Generate epetra map for nodes
    overlappedMapG = Teuchos::rcp(new Epetra_Map(-1,overlappedNodes,overlappedGIDs,0,Comm));
    }
    // build Export/Import
     Teuchos::RCP<Teuchos::Time> timerBuildOverlapMaps2 = Teuchos::TimeMonitor::getNewTimer("Build overlapped maps: exporter");
    {
    Teuchos::TimeMonitor timerBuildOverlapMaps2L(*timerBuildOverlapMaps2);
    exporter = Teuchos::rcp(new Epetra_Export(*overlappedMapG,*globalMapG));
    }
    Teuchos::RCP<Teuchos::Time> timerBuildOverlapMaps3 = Teuchos::TimeMonitor::getNewTimer("Build overlapped maps: importer");
    {
    Teuchos::TimeMonitor timerBuildOverlapMaps3L(*timerBuildOverlapMaps3);
    Epetra_Import importer(*overlappedMapG,*globalMapG);
    }
    }

/**********************************************************************************/
/********************* BUILD GRAPH FOR OVERLAPPED SOLUTION *************************/
/**********************************************************************************/
 

    Teuchos::RCP<Teuchos::Time> timerBuildOverlapGraph = Teuchos::TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: Total Time");
    Teuchos::RCP<Epetra_Vector> rhsVector;
    Teuchos::RCP<Epetra_CrsMatrix> StiffMatrix;
    Teuchos::RCP<Epetra_CrsMatrix> gl_StiffMatrix;
    Teuchos::RCP<Epetra_Vector> gl_rhsVector;
    { 
    Teuchos::TimeMonitor timerBuildOverlapGraphL(*timerBuildOverlapGraph); 
    
    Teuchos::RCP<Epetra_CrsGraph> overlappedGraph;
    Teuchos::RCP<Epetra_CrsGraph> ownedGraph;

    Teuchos::RCP<Teuchos::Time> timerBuildOverlapGraph2 = Teuchos::TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: create overlappedGraph & ownedGraph");
    {
    Teuchos::TimeMonitor timerBuildOverlapGraph2L(*timerBuildOverlapGraph2);
    overlappedGraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*overlappedMapG,0));
    ownedGraph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*globalMapG,0));
    }
    {
     Teuchos::RCP<Teuchos::Time> timerBuildOverlapGraph3 = Teuchos::TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: insertGlobalIndices overlappedGraph");
    {
    Teuchos::TimeMonitor timerBuildOverlapGraph3L(*timerBuildOverlapGraph3);
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
            int globalRow = globalNodeIds[localRow];

            // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
            for (int cellCol = 0; cellCol < numFieldsG; cellCol++){
    
              int localCol  = elemToNode(cell, cellCol);
              int globalCol = globalNodeIds[localCol];
    
              overlappedGraph->InsertGlobalIndices(globalRow,1, &globalCol);
    
            }// *** cell col loop ***
          }// *** cell row loop ***
        }// *** workset cell loop **
      }// *** workset loop ***
    }
    }
    
    Teuchos::RCP<Teuchos::Time> timerBuildOverlapGraph3 = Teuchos::TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: fillcomplete overlappedGraphT");
    {
    Teuchos::TimeMonitor timerBuildOverlapGraph3L(*timerBuildOverlapGraph3);
    overlappedGraph->FillComplete();
    }
   
    Teuchos::RCP<Teuchos::Time> timerBuildOverlapGraph4 = Teuchos::TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: export and fillcomplete ownedGraphT");
    {
    Teuchos::TimeMonitor timerBuildOverlapGraph4L(*timerBuildOverlapGraph4);
    // build global map
    ownedGraph->Export(*overlappedGraph, *exporter, Insert);
    ownedGraph->FillComplete();
    }
   
    Teuchos::RCP<Teuchos::Time> timerBuildOverlapGraph5 = Teuchos::TimeMonitor::getNewTimer("Build graphs for overlapped and owned solutions: create and zero global/local matrices/vectors");
    {
    Teuchos::TimeMonitor timerBuildOverlapGraph5L(*timerBuildOverlapGraph5);
    gl_StiffMatrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *ownedGraph));
    gl_StiffMatrix->FillComplete();
    gl_rhsVector = Teuchos::rcp(new Epetra_Vector(*globalMapG));

    StiffMatrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *overlappedGraph));
    StiffMatrix->FillComplete();
    rhsVector = Teuchos::rcp(new Epetra_Vector(*overlappedMapG));
    StiffMatrix->PutScalar(0.0);
    rhsVector->PutScalar(0.0);
    }
    }

#ifdef DUMP_DATA
/**********************************************************************************/
/**** PUT COORDINATES AND NODAL VALUES IN ARRAYS FOR OUTPUT (FOR PLOTTING ONLY) ***/
/**********************************************************************************/

  // Put coordinates in multivector for output
    Epetra_MultiVector nCoord(*globalMapG,3);
    Epetra_MultiVector nBound(*globalMapG,1);

     int indOwned = 0;
     for (int inode=0; inode<numNodes; inode++) {
       if (nodeIsOwned[inode]) {
          nCoord[0][indOwned]=nodeCoord(inode,0);
          nCoord[1][indOwned]=nodeCoord(inode,1);
          nCoord[2][indOwned]=nodeCoord(inode,2);
          nBound[0][indOwned]=nodeOnBoundary(inode);
          indOwned++;
       }
     }
     EpetraExt::MultiVectorToMatrixMarketFile("coords.dat",nCoord,0,0,false);
     EpetraExt::MultiVectorToMatrixMarketFile("nodeOnBound.dat",nBound,0,0,false);

    // Put element to node mapping in multivector for output
     Epetra_Map   globalMapElem(numElemsGlobal, numElems, 0, Comm);
     Epetra_MultiVector elem2node(globalMapElem, numNodesPerElem);
     for (int ielem=0; ielem<numElems; ielem++) {
        for (int inode=0; inode<numNodesPerElem; inode++) {
          elem2node[inode][ielem]=globalNodeIds[elemToNode(ielem,inode)];
        }
      }
     EpetraExt::MultiVectorToMatrixMarketFile("elem2node.dat",elem2node,0,0,false);

    if(MyPID==0) {Time.ResetStartTime();}

#endif

/**********************************************************************************/
/************************** DIRICHLET BC SETUP ************************************/
/**********************************************************************************/

  Teuchos::RCP<Teuchos::Time> timerDirichletBCs = Teuchos::TimeMonitor::getNewTimer("Get Dirichlet boundary values: Total Time");
  int numBCNodes = 0;
  int * BCNodes;
  Teuchos::RCP<Epetra_MultiVector> v;
  {
  Teuchos::TimeMonitor timerDirichletBCsL(*timerDirichletBCs);
  for (int inode = 0; inode < numNodes; inode++){
     if (nodeOnBoundary(inode) && nodeIsOwned[inode]){
        numBCNodes++;
     }
  }
   Teuchos::RCP<Teuchos::Time> timerDirichletBC1 = Teuchos::TimeMonitor::getNewTimer("Get Dirichlet boundary values: create and zero v");
   {
   Teuchos::TimeMonitor timerDirichletBC1L(*timerDirichletBC1);
  // Vector for use in applying BCs
   v = Teuchos::rcp(new Epetra_MultiVector(*globalMapG,true));
   v->PutScalar(0.0);
   }
   Teuchos::RCP<Teuchos::Time> timerDirichletBC2 = Teuchos::TimeMonitor::getNewTimer("Get Dirichlet boundary values: set v for Dirichlet nodes");
   {
   Teuchos::TimeMonitor timerDirichletBC2L(*timerDirichletBC2);
   // Set v to boundary values on Dirichlet nodes
    BCNodes = new int [numBCNodes];
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
           (*v)[0][iOwned]=exactSolution(x, y, z);
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
    std::cout << "Building discretization matrix and right hand side... \n\n";
    std::cout << "\tDesired workset size:                 " << desiredWorksetSize <<"\n";
    std::cout << "\tNumber of worksets (per processor):   " << numWorksets <<"\n\n";
    Time.ResetStartTime();
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

    Teuchos::RCP<Teuchos::Time> timerAssembleGlobalMatrix = Teuchos::TimeMonitor::getNewTimer("Assemble global matrix: Total Time");
    {
    Teuchos::TimeMonitor timerAssembleGlobalMatrixL(*timerAssembleGlobalMatrix); 
    
    //"WORKSET CELL" loop: local cell ordinal is relative to numElems
    for(int cell = worksetBegin; cell < worksetEnd; cell++){

      // Compute cell ordinal relative to the current workset
      int worksetCellOrdinal = cell - worksetBegin;

      // "CELL EQUATION" loop for the workset cell: cellRow is relative to the cell DoF numbering
      for (int cellRow = 0; cellRow < numFieldsG; cellRow++){

        int localRow  = elemToNode(cell, cellRow);
        int globalRow = globalNodeIds[localRow];
        double sourceTermContribution =  worksetRHS(worksetCellOrdinal, cellRow);

        rhsVector->SumIntoGlobalValues(1, &sourceTermContribution, &globalRow);

        // "CELL VARIABLE" loop for the workset cell: cellCol is relative to the cell DoF numbering
        for (int cellCol = 0; cellCol < numFieldsG; cellCol++){

          int localCol  = elemToNode(cell, cellCol);
          int globalCol = globalNodeIds[localCol];
          double operatorMatrixContribution = worksetStiffMatrix(worksetCellOrdinal, cellRow, cellCol);

          StiffMatrix->SumIntoGlobalValues(globalRow, 1, &operatorMatrixContribution, &globalCol);

        }// *** cell col loop ***
      }// *** cell row loop ***
    }// *** workset cell loop **
   //timerAssembleGlobalMatrix->stop(); 
   }//timer stop
  }// *** workset loop ***

/**********************************************************************************/
/********************* ASSEMBLE OVER MULTIPLE PROCESSORS **************************/
/**********************************************************************************/

   Teuchos::RCP<Teuchos::Time> timerAssembMultProc = Teuchos::TimeMonitor::getNewTimer("Global assembly over multiple processors: Total Time");
   {
   Teuchos::TimeMonitor timerAssembMultProcL(*timerAssembMultProc);   
   gl_StiffMatrix->PutScalar(0.0);
   gl_StiffMatrix->Export(*StiffMatrix,*exporter,Add);

   gl_rhsVector->PutScalar(0.0);
   gl_rhsVector->Export(*rhsVector,*exporter,Add);
   }

/**********************************************************************************/
/******************************* ADJUST MATRIX DUE TO BC **************************/
/**********************************************************************************/

  Teuchos::RCP<Teuchos::Time> timerAdjustMatrixBC = Teuchos::TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: Total Time");
  //timerAdjustMatrixBC->start(); 
  {
  Teuchos::TimeMonitor timerAdjustMatrixBCL(*timerAdjustMatrixBC); 
  Teuchos::RCP<Epetra_MultiVector> rhsDir;
  Teuchos::RCP<Teuchos::Time> timerAdjustMatrixBC1 = Teuchos::TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: apply stiffmatrix to v");
   {
   Teuchos::TimeMonitor timerAdjustMatrixBC1L(*timerAdjustMatrixBC1);
  // Apply stiffness matrix to v
   rhsDir = Teuchos::rcp(new Epetra_MultiVector(*globalMapG,true));
   gl_StiffMatrix->Apply(*v,*rhsDir);
   }
   
   Teuchos::RCP<Teuchos::Time> timerAdjustMatrixBC2 = Teuchos::TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: update rhs");
   {
   Teuchos::TimeMonitor timerAdjustMatrixBC2L(*timerAdjustMatrixBC2);
  // Update right-hand side
   gl_rhsVector->Update(-1.0,*rhsDir,1.0);
   }
   
   Teuchos::RCP<Teuchos::Time> timerAdjustMatrixBC3 = Teuchos::TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: adjust rhs to Dirichlet BCs");
   {
   Teuchos::TimeMonitor timerAdjustMatrixBC3L(*timerAdjustMatrixBC3);
    // Adjust rhs due to Dirichlet boundary conditions
   for (int inode=0; inode<numNodes; inode++){
      if (nodeIsOwned[inode]){
        if (nodeOnBoundary(inode)){
           // get global node number
           int gni = globalNodeIds[inode];
           int lid = globalMapG->LID(gni);
           double v_val = (*(*v)(0))[lid];
           gl_rhsVector->ReplaceGlobalValues(1,&v_val,&gni);
        }
      }
    }
   }

   Teuchos::RCP<Teuchos::Time> timerAdjustMatrixBC4 = Teuchos::TimeMonitor::getNewTimer("Adjust global matrix and rhs due to BCs: zero out rows/cols of stiffness for Dirichlet edges");
   {
   Teuchos::TimeMonitor timerAdjustMatrixBC4L(*timerAdjustMatrixBC4);
   // Zero out rows and columns of stiffness matrix corresponding to Dirichlet edges
   //  and add one to diagonal.
    ML_Epetra::Apply_OAZToMatrix(BCNodes, numBCNodes, *gl_StiffMatrix);
    }
    }

    delete [] BCNodes;
    
   //output global stiffness matrix and rhs vector to MatrixMarket file
   EpetraExt::RowMatrixToMatrixMarketFile("gl_StiffMatrix.dat",*gl_StiffMatrix);
   EpetraExt::VectorToMatrixMarketFile("gl_rhsVector.dat",*gl_rhsVector);

   //summarize timings
   Teuchos::TimeMonitor::summarize( std::cout );
  
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
   delete [] ownedGIDs;
   delete [] overlappedGIDs;
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

/**********************************************************************************/
/******************************* TEST ML ******************************************/
/**********************************************************************************/

// Test ML
int TestMultiLevelPreconditionerLaplace(char ProblemType[],
				 Teuchos::ParameterList   & MLList,
                                 Epetra_CrsMatrix   & A,
                                 const Epetra_MultiVector & xexact,
                                 Epetra_MultiVector & b,
                                 Epetra_MultiVector & uh,
                                 double & TotalErrorResidual,
				 double & TotalErrorExactSol)
{
  Epetra_MultiVector x(xexact);
  x.PutScalar(0.0);
  
  Epetra_LinearProblem Problem(&A,&x,&b); 
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  
  Epetra_Time Time(A.Comm());
  
  // =================== //
  // call ML and AztecOO //
  // =================== //
  
  AztecOO solver(Problem);  
  ML_Epetra::MultiLevelPreconditioner *MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 1);

  solver.Iterate(200, 1e-10);
  
  delete MLPrec;

  uh = *lhs;
  
  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //  
  double d = 0.0, d_tot = 0.0;  
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - xexact[0][i]) * ((*lhs)[0][i] - xexact[0][i]);
  
  A.Comm().SumAll(&d,&d_tot,1);
  
  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A.Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);
  
  string msg = ProblemType;
  
  if (A.Comm().MyPID() == 0) {
    cout << msg << endl << "......Using " << A.Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
  
  return( solver.NumIters() );
  
}


